#include "sim/nbody_simulation.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace galaxy {

namespace {

using units::kG;

constexpr double kTwoPi = 2.0 * std::numbers::pi;

// Portable uniform [0,1) from raw mt19937_64 bits.
double uniform01(std::mt19937_64& rng) {
  return static_cast<double>(rng() >> 11) * (1.0 / 9007199254740992.0);  // 2^-53
}

// Standard normal via Box-Muller (one value per call; simple beats fast here).
double gaussian(std::mt19937_64& rng) {
  const double u1 = 1.0 - uniform01(rng);  // (0, 1]
  const double u2 = uniform01(rng);
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
}

// Exponential disk: surface density ~ exp(-r/h), so the mass inside radius
// x*h is P(x) = 1 - (1 + x) exp(-x). Invert P(x) = u by bisection.
double exponential_radius_in_scale_lengths(double u, double x_max) {
  double lo = 0.0, hi = x_max;
  for (int it = 0; it < 60; ++it) {
    const double mid = 0.5 * (lo + hi);
    const double p = 1.0 - (1.0 + mid) * std::exp(-mid);
    (p < u ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

double dynamical_time(double r, double gm) { return std::sqrt(r * r * r / gm); }

}  // namespace

// One galaxy's masses and analytic enclosed-mass profile (in the plane). The
// circular speeds of the disk and the Jeans dispersion of the halo are both
// computed from this, so every component feels the same total potential.
struct NBodySimulation::GalaxyModel {
  double m_bh = 0, m_disk = 0, m_halo = 0;
  std::size_t n_disk = 0, n_halo = 0;
  bool counterclockwise = true;
  double radius = 0;  // disk radius

  DiskProfile profile = DiskProfile::Uniform;
  double h = 0, x_max = 0, p_max = 1;  // exponential disk
  double r_in = 0, p_in = 0;           // central hole radius, and its exponential mass fraction

  double halo_core = 0, halo_k = 0, halo_radius = 0;  // halo: Sigma ~ 1/sqrt(R^2 + core^2)

  double disk_fraction_inside(double r) const {
    if (r >= radius) return 1.0;
    if (r <= r_in) return 0.0;
    if (profile == DiskProfile::Uniform) return (r * r - r_in * r_in) / (radius * radius - r_in * r_in);
    const double x = r / h;
    return ((1.0 - (1.0 + x) * std::exp(-x)) - p_in) / (p_max - p_in);
  }
  double halo_fraction_inside(double r) const {
    if (m_halo <= 0.0) return 0.0;
    if (r >= halo_radius) return 1.0;
    if (r <= r_in) return 0.0;
    const double s_in = std::sqrt(r_in * r_in + halo_core * halo_core);
    return (std::sqrt(r * r + halo_core * halo_core) - s_in) / (halo_core + halo_k - s_in);
  }
  double enclosed(double r) const {
    return m_bh + m_disk * disk_fraction_inside(r) + m_halo * halo_fraction_inside(r);
  }

  // Effective enclosed mass M_eff(r) defined by the MEASURED mean inward in-plane
  // acceleration: a_in(r) = G M_eff(r) r / (r^2 + eps^2)^1.5. For flat (2D)
  // mass distributions the true in-plane force differs from the spherical
  // G M(<r)/r^2 (by ~1.5x for a 1/R sheet), and using the spherical value leaves
  // the halo under-supported so the whole galaxy collapses and rebounds.
  std::vector<double> eff_r, eff_m;

  double effective_mass(double r) const {
    if (eff_r.empty()) return enclosed(r);
    if (r <= eff_r.front()) return eff_m.front();
    if (r >= eff_r.back()) return eff_m.back();
    const auto it = std::upper_bound(eff_r.begin(), eff_r.end(), r);
    const auto j = static_cast<std::size_t>(it - eff_r.begin());
    const double f = (r - eff_r[j - 1]) / (eff_r[j] - eff_r[j - 1]);
    return eff_m[j - 1] * (1.0 - f) + eff_m[j] * f;
  }

  double sample_disk_radius(std::mt19937_64& rng) const {
    const double u = uniform01(rng);
    if (profile == DiskProfile::Uniform)
      return std::sqrt(r_in * r_in + u * (radius * radius - r_in * r_in));
    return h * exponential_radius_in_scale_lengths(p_in + u * (p_max - p_in), x_max);
  }
  double sample_halo_radius(std::mt19937_64& rng) const {
    const double s_in = std::sqrt(r_in * r_in + halo_core * halo_core);
    const double s = s_in + uniform01(rng) * (halo_core + halo_k - s_in);
    return std::sqrt(std::max(0.0, s * s - halo_core * halo_core));
  }

  // Build the table from a throwaway realization (more, lighter particles to
  // beat shot noise), pushed through the same tree, theta and softening as the
  // real integration.
  void measure_profile(double eps2, double theta, std::uint64_t seed) {
    constexpr std::size_t kMaxBodies = 150000;
    const std::size_t base = n_disk + n_halo;
    const std::size_t boost = std::clamp<std::size_t>(kMaxBodies / std::max<std::size_t>(base, 1), 1, 4);
    const std::size_t nd = n_disk * boost, nh = n_halo * boost;

    std::mt19937_64 rng(seed ^ 0x9E3779B97F4A7C15ull);
    std::vector<double> x, y, m;
    x.reserve(nd + nh + 1);
    y.reserve(nd + nh + 1);
    m.reserve(nd + nh + 1);
    auto add = [&](double r, double mass) {
      const double th = kTwoPi * uniform01(rng);
      x.push_back(r * std::cos(th));
      y.push_back(r * std::sin(th));
      m.push_back(mass);
    };
    for (std::size_t k = 0; k < nd; ++k) add(sample_disk_radius(rng), m_disk / static_cast<double>(nd));
    for (std::size_t k = 0; k < nh; ++k) add(sample_halo_radius(rng), m_halo / static_cast<double>(nh));
    const std::size_t n_test = x.size();
    if (m_bh > 0.0) {
      x.push_back(0.0);
      y.push_back(0.0);
      m.push_back(m_bh);
    }

    BarnesHutTree tree;
    tree.build(x.data(), y.data(), m.data(), x.size());
    std::vector<std::pair<double, double>> samples(n_test);  // (r, M_eff)
#pragma omp parallel for if (n_test > 4096) schedule(dynamic, 256)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n_test); ++i) {
      double ax, ay;
      tree.accel(static_cast<std::int32_t>(i), x[i], y[i], std::min(theta, 0.4), eps2, ax, ay);
      const double r = std::hypot(x[i], y[i]);
      const double a_in = r > 0.0 ? -(ax * x[i] + ay * y[i]) / r : 0.0;
      const double d2 = r * r + eps2;
      samples[i] = {r, a_in * d2 * std::sqrt(d2) / (kG * std::max(r, 1e-30))};
    }
    std::sort(samples.begin(), samples.end());

    // Equal-count bins; each contributes (mean r, mean M_eff).
    const std::size_t per_bin = std::max<std::size_t>(200, n_test / 80);
    eff_r.clear();
    eff_m.clear();
    eff_r.push_back(0.0);
    eff_m.push_back(std::max(m_bh, 0.0));
    for (std::size_t s = 0; s < n_test; s += per_bin) {
      const std::size_t e = std::min(n_test, s + per_bin);
      if (e - s < per_bin / 2 && !eff_r.empty() && s > 0) break;  // drop a tiny tail bin
      double sr = 0, sm = 0;
      for (std::size_t k = s; k < e; ++k) {
        sr += samples[k].first;
        sm += samples[k].second;
      }
      const double cnt = static_cast<double>(e - s);
      if (sr / cnt > eff_r.back()) {
        eff_r.push_back(sr / cnt);
        eff_m.push_back(std::max(sm / cnt, 0.0));
      }
    }
  }
};

// ---------------------------------------------------------------------------
// Construction / initial conditions
// ---------------------------------------------------------------------------

NBodySimulation::NBodySimulation(const ScenarioParams& params, ProgressFn progress)
    : params_(params) {
  if (params_.mass_model != MassModel::Distributed)
    throw std::invalid_argument("NBodySimulation: scenario is not in the distributed-mass mode");
  if (const std::string err = validate_scenario(params_); !err.empty())
    throw std::invalid_argument(err);

  const DistributedParams& d = params_.distributed;
  const bool two = params_.galaxy_b.has_value();
  const GalaxyParams* galaxies[2] = {&params_.galaxy_a, two ? &*params_.galaxy_b : nullptr};
  const int n_galaxies = two ? 2 : 1;

  // Per-galaxy models.
  GalaxyModel models[2];
  for (int k = 0; k < n_galaxies; ++k) {
    const GalaxyParams& g = *galaxies[k];
    GalaxyModel& m = models[k];
    m.m_bh = d.bh_mass_fraction * g.total_mass_msun;
    m.m_halo = d.halo_mass_fraction * g.total_mass_msun;
    m.m_disk = g.total_mass_msun - m.m_bh - m.m_halo;
    m.n_disk = static_cast<std::size_t>(g.num_stars);
    m.n_halo = m.m_halo > 0.0
                   ? std::max<std::size_t>(
                         1, static_cast<std::size_t>(std::llround(
                                static_cast<double>(m.n_disk) * d.halo_particles_per_disk_particle)))
                   : 0;
    m.counterclockwise = g.counterclockwise;
    m.radius = g.radius_ly;
    m.r_in = d.inner_hole_fraction * g.radius_ly;
    m.profile = d.profile;
    if (d.profile == DiskProfile::Exponential) {
      m.h = d.scale_length_fraction * g.radius_ly;
      m.x_max = g.radius_ly / m.h;
      m.p_max = 1.0 - (1.0 + m.x_max) * std::exp(-m.x_max);
      const double x_in = m.r_in / m.h;
      m.p_in = 1.0 - (1.0 + x_in) * std::exp(-x_in);
    }
    if (m.m_halo > 0.0) {
      m.halo_core = d.halo_core_fraction * g.radius_ly;
      m.halo_radius = d.halo_radius_fraction * g.radius_ly;
      m.halo_k = std::sqrt(m.halo_radius * m.halo_radius + m.halo_core * m.halo_core) - m.halo_core;
    }
  }

  // Softening: a fraction of the mean particle spacing of the sparsest
  // component, sqrt(area / N), so macro-particles do not scatter off each other.
  double spacing = 0.0;
  for (int k = 0; k < n_galaxies; ++k) {
    const GalaxyModel& m = models[k];
    spacing = std::max(spacing, m.radius * std::sqrt(std::numbers::pi / static_cast<double>(m.n_disk)));
    if (m.n_halo > 0)
      spacing = std::max(spacing, m.halo_radius *
                                      std::sqrt(std::numbers::pi / static_cast<double>(m.n_halo)));
  }
  eps_ = d.softening_spacing_factor * spacing;
  eps2_ = eps_ * eps_;

  for (int k = 0; k < n_galaxies; ++k)
    models[k].measure_profile(eps2_, d.opening_angle, params_.seed + static_cast<std::uint64_t>(k));

  // Galaxy centres and bulk velocities: centre-of-mass frame of the two galaxies,
  // weighted by their TOTAL masses.
  Vec2 pos[2], vel[2];
  if (two) {
    const auto& e = params_.encounter;
    const double dx = std::sqrt(e.separation_ly * e.separation_ly -
                                e.impact_parameter_ly * e.impact_parameter_ly);
    const Vec2 rel_pos{dx, e.impact_parameter_ly};
    const Vec2 rel_vel{-e.relative_speed_km_s * units::kKmPerSecToLyPerYr, 0.0};
    const double ma = params_.galaxy_a.total_mass_msun, mb = params_.galaxy_b->total_mass_msun;
    pos[0] = rel_pos * (-mb / (ma + mb));
    pos[1] = rel_pos * (ma / (ma + mb));
    vel[0] = rel_vel * (-mb / (ma + mb));
    vel[1] = rel_vel * (ma / (ma + mb));
  }

  // Disks, then halos, then black holes (see the layout note in the header).
  std::mt19937_64 rng(params_.seed);
  for (int k = 0; k < n_galaxies; ++k)
    append_disk(models[k], static_cast<std::uint8_t>(k), pos[k], vel[k], rng);
  n_disk_ = x_.size();
  for (int k = 0; k < n_galaxies; ++k)
    append_halo(models[k], static_cast<std::uint8_t>(k), pos[k], vel[k], rng);
  n_halo_ = x_.size() - n_disk_;

  for (int k = 0; k < n_galaxies; ++k) {
    if (models[k].m_bh <= 0.0) continue;
    x_.push_back(pos[k].x);
    y_.push_back(pos[k].y);
    vx_.push_back(vel[k].x);
    vy_.push_back(vel[k].y);
    m_.push_back(models[k].m_bh);
    gal_.push_back(static_cast<std::uint8_t>(k));
  }
  n_bh_ = x_.size() - n_disk_ - n_halo_;

  if (d.relax_myr > 0.0) {
    // Evolve each galaxy alone, then overwrite its sampled state with the settled
    // one (re-centred on the galaxy's own centre of mass and bulk velocity).
    // Coarser opening angle and timestep: this is only to settle the structure.
    constexpr double kChunkYears = 10.0e6;
    std::size_t disk_off = 0, halo_off = n_disk_, bh_off = n_disk_ + n_halo_;
    for (int k = 0; k < n_galaxies; ++k) {
      ScenarioParams ip;
      ip.mass_model = MassModel::Distributed;
      ip.galaxy_a = *galaxies[k];
      ip.distributed = d;
      ip.distributed.relax_myr = 0.0;
      ip.distributed.opening_angle = std::max(d.opening_angle, 0.8);
      ip.step_accuracy = std::max(params_.step_accuracy, 0.1);
      ip.seed = params_.seed + static_cast<std::uint64_t>(k);
      NBodySimulation iso(ip);

      const double total = d.relax_myr * 1.0e6;
      for (double done = 0.0; done < total; done += kChunkYears) {
        iso.advance(std::min(kChunkYears, total - done), INT_MAX);
        if (progress) progress((k + (done + kChunkYears) / total) / n_galaxies);
      }

      double mtot = 0, cx = 0, cy = 0, cvx = 0, cvy = 0;
      for (std::size_t j = 0; j < iso.x_.size(); ++j) {
        mtot += iso.m_[j];
        cx += iso.m_[j] * iso.x_[j];
        cy += iso.m_[j] * iso.y_[j];
        cvx += iso.m_[j] * iso.vx_[j];
        cvy += iso.m_[j] * iso.vy_[j];
      }
      cx /= mtot; cy /= mtot; cvx /= mtot; cvy /= mtot;
      auto put = [&](std::size_t dst, std::size_t src) {
        x_[dst] = pos[k].x + (iso.x_[src] - cx);
        y_[dst] = pos[k].y + (iso.y_[src] - cy);
        vx_[dst] = vel[k].x + (iso.vx_[src] - cvx);
        vy_[dst] = vel[k].y + (iso.vy_[src] - cvy);
      };
      for (std::size_t j = 0; j < models[k].n_disk; ++j) put(disk_off + j, j);
      for (std::size_t j = 0; j < models[k].n_halo; ++j)
        put(halo_off + j, models[k].n_disk + j);
      if (models[k].m_bh > 0.0) put(bh_off++, models[k].n_disk + models[k].n_halo);
      disk_off += models[k].n_disk;
      halo_off += models[k].n_halo;
    }
  }

  ax_.assign(x_.size(), 0.0);
  ay_.assign(x_.size(), 0.0);

  const double amax2 = compute_forces();
  dt_ = compute_initial_dt(amax2);
  aref2_ = amax2;
  e0_ = total_energy();

  // Camera extents: visible bodies (disk + black holes) and everything.
  const Vec2 com = center_of_mass();
  auto farthest2 = [&](std::size_t i) {
    const double dx = x_[i] - com.x, dy = y_[i] - com.y;
    return dx * dx + dy * dy;
  };
  double vis2 = 0.0, all2 = 0.0;
  std::vector<double> vis_r2, all_r2;
  const bool robust = d.relax_myr > 0.0;  // relaxation flings a few stragglers far out
  for (std::size_t i = 0; i < x_.size(); ++i) {
    const double r2 = farthest2(i);
    const bool is_halo = i >= n_disk_ && i < n_disk_ + n_halo_;
    all2 = std::max(all2, r2);
    if (!is_halo) vis2 = std::max(vis2, r2);
    if (robust) {
      all_r2.push_back(r2);
      if (!is_halo) vis_r2.push_back(r2);
    }
  }
  if (robust) {
    // 99.5th percentile (radius is monotonic in r2), with a little margin.
    auto pct = [](std::vector<double>& v) {
      const std::size_t k = std::min(v.size() - 1, static_cast<std::size_t>(0.995 * static_cast<double>(v.size())));
      std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
      return v[k] * 1.1025;
    };
    vis2 = pct(vis_r2);
    all2 = pct(all_r2);
  }
  extent_ = std::sqrt(vis2);
  dm_extent_ = std::sqrt(all2);
}

// The visible disk. Circular speeds come from the analytic enclosed mass of ALL
// components (black hole + disk + halo) in the SOFTENED force law that is
// actually integrated. Dispersion adds a Gaussian velocity, as a fraction of the
// local circular speed, to each component.
void NBodySimulation::append_disk(const GalaxyModel& g, std::uint8_t gid, Vec2 centre,
                                  Vec2 bulk_vel, std::mt19937_64& rng) {
  const DistributedParams& d = params_.distributed;
  const std::size_t n = g.n_disk;
  const double m_particle = g.m_disk / static_cast<double>(n);

  const std::size_t base = x_.size();
  for (auto* v : {&x_, &y_, &vx_, &vy_, &m_}) v->resize(base + n);
  gal_.resize(base + n, gid);

  const double spin = g.counterclockwise ? 1.0 : -1.0;
  double sx = 0, sy = 0, svx = 0, svy = 0;

  for (std::size_t k = 0; k < n; ++k) {
    const double r = g.sample_disk_radius(rng);
    const double th = kTwoPi * uniform01(rng);
    const double st = std::sin(th), ct = std::cos(th);

    const double d2 = r * r + eps2_;
    const double v_circ = std::sqrt(kG * g.effective_mass(r) * r * r / (d2 * std::sqrt(d2)));

    double vx = -spin * v_circ * st;
    double vy = spin * v_circ * ct;
    if (d.dispersion_fraction > 0.0) {
      const double sigma = d.dispersion_fraction * v_circ;
      vx += sigma * gaussian(rng);
      vy += sigma * gaussian(rng);
    }

    const std::size_t i = base + k;
    x_[i] = r * ct;
    y_[i] = r * st;
    vx_[i] = vx;
    vy_[i] = vy;
    m_[i] = m_particle;
    sx += x_[i];
    sy += y_[i];
    svx += vx;
    svy += vy;
  }

  // Remove the sampling noise in the group's mean position and velocity so the
  // black hole sits at the disk's centre of mass, then place the galaxy.
  const double inv = 1.0 / static_cast<double>(n);
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t i = base + k;
    x_[i] += centre.x - sx * inv;
    y_[i] += centre.y - sy * inv;
    vx_[i] += bulk_vel.x - svx * inv;
    vy_[i] += bulk_vel.y - svy * inv;
  }
}

// The dark-matter halo: surface density ~ 1/sqrt(R^2 + core^2) out to the halo
// radius, no net rotation, isotropic in-plane dispersion from the radial Jeans
// equation of a non-rotating 2D system,
//     d(Sigma sigma^2)/dR = -Sigma * a_R(R)
// integrated inward from the truncation radius (where the pressure is zero).
// a_R is the radial acceleration from the total enclosed mass.
void NBodySimulation::append_halo(const GalaxyModel& g, std::uint8_t gid, Vec2 centre,
                                  Vec2 bulk_vel, std::mt19937_64& rng) {
  const std::size_t n = g.n_halo;
  if (n == 0) return;
  const double m_particle = g.m_halo / static_cast<double>(n);
  const double a = g.halo_core, rh = g.halo_radius;

  // Tabulate sigma^2(R) on a uniform grid by trapezoid, from the outside in.
  constexpr int kGrid = 1024;
  const double dr = rh / kGrid;
  std::vector<double> sigma2(kGrid + 1, 0.0);
  auto integrand = [&](double r) {
    const double d2 = r * r + eps2_;
    const double a_r = kG * g.effective_mass(r) * r / (d2 * std::sqrt(d2));
    return a_r / std::sqrt(r * r + a * a);  // Sigma(R) * a_R(R), Sigma up to a constant
  };
  double integral = 0.0;
  for (int k = kGrid - 1; k >= 0; --k) {
    const double r0 = k * dr, r1 = (k + 1) * dr;
    integral += 0.5 * (integrand(r0) + integrand(r1)) * dr;
    sigma2[k] = integral * std::sqrt(r0 * r0 + a * a);  // divide by Sigma(r0)
  }

  const std::size_t base = x_.size();
  for (auto* v : {&x_, &y_, &vx_, &vy_, &m_}) v->resize(base + n);
  gal_.resize(base + n, gid);

  double sx = 0, sy = 0, svx = 0, svy = 0;
  for (std::size_t k = 0; k < n; ++k) {
    // Invert the cumulative mass fraction (outside the central hole).
    const double r = g.sample_halo_radius(rng);
    const double th = kTwoPi * uniform01(rng);

    const double t = std::min(r / dr, static_cast<double>(kGrid) - 1e-9);
    const auto idx = static_cast<std::size_t>(t);
    const double frac = t - static_cast<double>(idx);
    const double sig = std::sqrt(std::max(0.0, sigma2[idx] * (1.0 - frac) + sigma2[idx + 1] * frac));

    const std::size_t i = base + k;
    x_[i] = r * std::cos(th);
    y_[i] = r * std::sin(th);
    vx_[i] = sig * gaussian(rng);
    vy_[i] = sig * gaussian(rng);
    m_[i] = m_particle;
    sx += x_[i];
    sy += y_[i];
    svx += vx_[i];
    svy += vy_[i];
  }

  const double inv = 1.0 / static_cast<double>(n);
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t i = base + k;
    x_[i] += centre.x - sx * inv;
    y_[i] += centre.y - sy * inv;
    vx_[i] += bulk_vel.x - svx * inv;
    vy_[i] += bulk_vel.y - svy * inv;
  }
}

double NBodySimulation::compute_initial_dt(double max_acc2) const {
  const double eta = params_.step_accuracy;
  // Acceleration criterion (standard N-body form): dt = eta * sqrt(eps / |a|max).
  double dt = eta * std::sqrt(eps_ / std::sqrt(std::max(max_acc2, 1e-300)));

  // Two galaxies: also resolve the pericentre of the galaxy pair, treating each
  // galaxy as a point mass. The floor keeps a near head-on setup affordable;
  // step() shrinks dt on the fly if forces turn out stronger than assumed.
  if (params_.galaxy_b) {
    const auto& e = params_.encounter;
    const double ma = params_.galaxy_a.total_mass_msun, mb = params_.galaxy_b->total_mass_msun;
    const double mu = kG * (ma + mb);
    const double r = e.separation_ly;
    const double v = e.relative_speed_km_s * units::kKmPerSecToLyPerYr;
    const double h = std::abs(e.impact_parameter_ly) * v;
    const double en = 0.5 * v * v - mu / r;
    const double ecc = std::sqrt(std::max(0.0, 1.0 + 2.0 * en * h * h / (mu * mu)));
    const double r_peri = (h * h / mu) / (1.0 + ecc);
    const double r_floor = std::max(
        eps_, 0.01 * std::min(params_.galaxy_a.radius_ly, params_.galaxy_b->radius_ly));
    dt = std::min(dt, eta * dynamical_time(std::max(r_peri, r_floor), mu));
  }
  return dt;
}

// ---------------------------------------------------------------------------
// Forces and integration
// ---------------------------------------------------------------------------

double NBodySimulation::compute_forces() {
  tree_.build(x_.data(), y_.data(), m_.data(), x_.size());

  const double theta = params_.distributed.opening_angle;
  const double eps2 = eps2_;
  const double* x = x_.data();
  const double* y = y_.data();
  double* ax = ax_.data();
  double* ay = ay_.data();
  const auto n = static_cast<std::ptrdiff_t>(x_.size());
  const BarnesHutTree& tree = tree_;
  const std::int32_t* order = tree.leaf_order().data();  // spatially coherent target order

  double amax2 = 0.0;
#pragma omp parallel for reduction(max : amax2) if (n > 2048) schedule(dynamic, 256)
  for (std::ptrdiff_t k = 0; k < n; ++k) {
    const std::int32_t i = order[k];
    double fx, fy;
    tree.accel(i, x[i], y[i], theta, eps2, fx, fy);
    ax[i] = fx;
    ay[i] = fy;
    amax2 = std::max(amax2, fx * fx + fy * fy);
  }
  return amax2;
}

void NBodySimulation::step(double dt) {
  const double half = 0.5 * dt;
  const auto n = static_cast<std::ptrdiff_t>(x_.size());
  double* x = x_.data();
  double* y = y_.data();
  double* vx = vx_.data();
  double* vy = vy_.data();
  const double* ax = ax_.data();
  const double* ay = ay_.data();

  // Kick (cached acceleration) + drift.
#pragma omp parallel for if (n > 8192) schedule(static)
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    vx[i] += half * ax[i];
    vy[i] += half * ay[i];
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
  }

  // New tree and accelerations at the new positions, then the second kick.
  const double amax2 = compute_forces();
#pragma omp parallel for if (n > 8192) schedule(static)
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    vx[i] += half * ax[i];
    vy[i] += half * ay[i];
  }

  time_ += dt;
  ++step_count_;

  // A fixed step sized from the starting state can become too coarse when two
  // black holes meet or a dense clump forms. Because dt ~ 1/sqrt(|a|), halve it
  // whenever the peak acceleration has grown 4x past what it was sized for. dt
  // only ever shrinks, so the integrator stays close to a fixed-step leapfrog.
  if (amax2 > 16.0 * aref2_) {
    dt_ *= 0.5;
    aref2_ = amax2;
  }
}

double NBodySimulation::run_queued(int max_steps) {
  auto n = static_cast<long long>(carry_ / dt_);
  if (n > max_steps) n = max_steps;
  double advanced = 0.0;
  // dt_ can shrink mid-run, so re-read it each step and account for the time
  // actually simulated.
  for (long long i = 0; i < n; ++i) {
    const double dt = dt_;
    step(dt);
    advanced += dt;
    carry_ -= dt;
  }
  return advanced;
}

void NBodySimulation::drop_backlog() { carry_ = std::min(carry_, dt_); }

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

const std::vector<double>& NBodySimulation::potentials() const {
  if (phi_step_ != step_count_ || phi_.size() != x_.size()) {
    phi_.resize(x_.size());
    const double theta = params_.distributed.opening_angle;
    const double eps2 = eps2_;
    const double* x = x_.data();
    const double* y = y_.data();
    double* phi = phi_.data();
    const auto n = static_cast<std::ptrdiff_t>(x_.size());
    const BarnesHutTree& tree = tree_;
    const std::int32_t* order = tree.leaf_order().data();
#pragma omp parallel for if (n > 2048) schedule(dynamic, 256)
    for (std::ptrdiff_t k = 0; k < n; ++k) {
      const std::int32_t i = order[k];
      phi[i] = tree.potential(i, x[i], y[i], theta, eps2);
    }
    phi_step_ = step_count_;
  }
  return phi_;
}

double NBodySimulation::total_energy() const {
  const std::vector<double>& phi = potentials();
  double ke = 0.0, pe = 0.0;
  for (std::size_t i = 0; i < x_.size(); ++i) {
    ke += 0.5 * m_[i] * (vx_[i] * vx_[i] + vy_[i] * vy_[i]);
    pe += m_[i] * phi[i];
  }
  return ke + 0.5 * pe;  // each pair appears twice in sum_i m_i phi_i
}

double NBodySimulation::total_energy_direct() const {
  const std::size_t n = x_.size();
  double ke = 0.0, pe = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    ke += 0.5 * m_[i] * (vx_[i] * vx_[i] + vy_[i] * vy_[i]);
    for (std::size_t j = i + 1; j < n; ++j) {
      const double dx = x_[i] - x_[j], dy = y_[i] - y_[j];
      pe -= kG * m_[i] * m_[j] / std::sqrt(dx * dx + dy * dy + eps2_);
    }
  }
  return ke + pe;
}

BoundFractions NBodySimulation::bound_fractions() const {
  BoundFractions out;
  out.galaxy = std::numeric_limits<double>::quiet_NaN();
  const std::size_t n = n_disk_;
  if (n == 0) return out;

  // Disk particles only (the "stars"; neither the halo nor the black holes).
  // Energy per unit mass in the system's centre-of-mass frame, against the
  // tree potential of everything else.
  const std::vector<double>& phi = potentials();
  const Vec2 p = total_momentum();
  double mtot = 0.0;
  for (double m : m_) mtot += m;
  const double vcx = p.x / mtot, vcy = p.y / mtot;

  std::size_t bound = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double rvx = vx_[i] - vcx, rvy = vy_[i] - vcy;
    if (0.5 * (rvx * rvx + rvy * rvy) + phi[i] < 0.0) ++bound;
  }
  out.system = static_cast<double>(bound) / static_cast<double>(n);
  return out;
}

Vec2 NBodySimulation::center_of_mass() const {
  double sx = 0, sy = 0, m = 0;
  for (std::size_t i = 0; i < x_.size(); ++i) {
    sx += m_[i] * x_[i];
    sy += m_[i] * y_[i];
    m += m_[i];
  }
  return {sx / m, sy / m};
}

Vec2 NBodySimulation::total_momentum() const {
  double px = 0, py = 0;
  for (std::size_t i = 0; i < x_.size(); ++i) {
    px += m_[i] * vx_[i];
    py += m_[i] * vy_[i];
  }
  return {px, py};
}

double NBodySimulation::angular_momentum() const {
  double l = 0.0;
  for (std::size_t i = 0; i < x_.size(); ++i) l += m_[i] * (x_[i] * vy_[i] - y_[i] * vx_[i]);
  return l;
}

double NBodySimulation::central_separation_ly() const {
  if (!params_.galaxy_b) return 0.0;
  double sx[2] = {0, 0}, sy[2] = {0, 0}, m[2] = {0, 0};
  for (std::size_t i = 0; i < x_.size(); ++i) {
    const int g = gal_[i];
    sx[g] += m_[i] * x_[i];
    sy[g] += m_[i] * y_[i];
    m[g] += m_[i];
  }
  const double dx = sx[1] / m[1] - sx[0] / m[0];
  const double dy = sy[1] / m[1] - sy[0] / m[0];
  return std::sqrt(dx * dx + dy * dy);
}

void NBodySimulation::fill_render_points(Vec2 view_center, double lag, float* xy) const {
  // Output order: disk particles, black holes, halo (see RenderCounts).
  std::size_t out = 0;
  auto put = [&](std::size_t i) {
    xy[2 * out] = static_cast<float>(x_[i] + vx_[i] * lag - view_center.x);
    xy[2 * out + 1] = static_cast<float>(y_[i] + vy_[i] * lag - view_center.y);
    ++out;
  };
  for (std::size_t i = 0; i < n_disk_; ++i) put(i);
  for (std::size_t i = n_disk_ + n_halo_; i < x_.size(); ++i) put(i);  // black holes
  for (std::size_t i = n_disk_; i < n_disk_ + n_halo_; ++i) put(i);    // halo
}

}  // namespace galaxy
