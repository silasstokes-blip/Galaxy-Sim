#include "sim/simulation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace galaxy {

namespace {

using units::kG;
using units::kStarMassMsun;

constexpr double kTwoPi = 2.0 * std::numbers::pi;

// Portable uniform [0,1) from raw mt19937_64 bits (std::uniform_real_distribution
// is not guaranteed to give identical results across standard libraries).
double uniform01(std::mt19937_64& rng) {
  return static_cast<double>(rng() >> 11) * (1.0 / 9007199254740992.0);  // 2^-53
}


// Acceleration at (x, y) from up to a handful of Plummer-softened point masses.
inline void accel_from_sources(const double* sx, const double* sy, const double* sgm,
                               int n, double eps2, double x, double y, double& ax,
                               double& ay) {
  double axs = 0.0, ays = 0.0;
  for (int c = 0; c < n; ++c) {
    const double dx = sx[c] - x;
    const double dy = sy[c] - y;
    const double r2 = dx * dx + dy * dy + eps2;
    const double inv = sgm[c] / (r2 * std::sqrt(r2));
    axs += dx * inv;
    ays += dy * inv;
  }
  ax = axs;
  ay = ays;
}

}  // namespace

// ---------------------------------------------------------------------------
// Construction / initial conditions
// ---------------------------------------------------------------------------

Simulation::Simulation(const ScenarioParams& params) : params_(params) {
  if (params_.mass_model != MassModel::Centralized)
    throw std::invalid_argument("Simulation: scenario is not in the centralized-mass mode");
  if (const std::string err = validate_scenario(params_); !err.empty())
    throw std::invalid_argument(err);
  eps2_ = params_.softening_ly * params_.softening_ly;

  build_centrals();

  std::mt19937_64 rng(params_.seed);
  populate_galaxy(params_.galaxy_a, centrals_[0], rng);
  if (params_.galaxy_b) populate_galaxy(*params_.galaxy_b, centrals_[1], rng);

  compute_central_acc();

  // Prime the cached star accelerations for the first half-kick.
  double sx[2], sy[2], sgm[2];
  const int nc = static_cast<int>(centrals_.size());
  for (int c = 0; c < nc; ++c) {
    sx[c] = centrals_[c].pos.x;
    sy[c] = centrals_[c].pos.y;
    sgm[c] = kG * centrals_[c].mass;
  }
  for (std::size_t i = 0; i < stars_.size(); ++i)
    accel_from_sources(sx, sy, sgm, nc, eps2_, stars_.x[i], stars_.y[i], stars_.ax[i],
                       stars_.ay[i]);

  dt_ = compute_dt();
  e0_ = total_energy();
  extent_ = compute_extent();
}

void Simulation::build_centrals() {
  CentralBody a;
  a.mass = central_mass_msun(params_.galaxy_a);
  centrals_.push_back(a);

  if (!params_.galaxy_b) return;

  CentralBody b;
  b.mass = central_mass_msun(*params_.galaxy_b);

  // Relative orbit: B relative to A starts at (sqrt(d^2 - b^2), b) moving
  // along -x. Then shift into the centre-of-mass frame.
  const auto& e = params_.encounter;
  const double dx = std::sqrt(e.separation_ly * e.separation_ly -
                              e.impact_parameter_ly * e.impact_parameter_ly);
  const Vec2 rel_pos{dx, e.impact_parameter_ly};
  const Vec2 rel_vel{-e.relative_speed_km_s * units::kKmPerSecToLyPerYr, 0.0};

  const double m_tot = a.mass + b.mass;
  a.pos = rel_pos * (-b.mass / m_tot);
  b.pos = rel_pos * (a.mass / m_tot);
  a.vel = rel_vel * (-b.mass / m_tot);
  b.vel = rel_vel * (a.mass / m_tot);

  centrals_[0] = a;
  centrals_.push_back(b);
}

void Simulation::populate_galaxy(const GalaxyParams& g, const CentralBody& c,
                                 std::mt19937_64& rng) {
  const auto n = static_cast<std::size_t>(g.num_stars);
  const std::size_t base = stars_.size();
  for (auto* v : {&stars_.x, &stars_.y, &stars_.vx, &stars_.vy, &stars_.ax, &stars_.ay})
    v->resize(base + n);

  const double r_in = g.inner_fraction * g.radius_ly;
  const double r_in2 = r_in * r_in;
  const double r_out2 = g.radius_ly * g.radius_ly;
  const double spin = g.counterclockwise ? 1.0 : -1.0;

  for (std::size_t i = 0; i < n; ++i) {
    // sqrt of a uniform variate in r^2 gives a uniform surface density.
    const double r = std::sqrt(r_in2 + uniform01(rng) * (r_out2 - r_in2));
    const double th = kTwoPi * uniform01(rng);
    const double s = std::sin(th), co = std::cos(th);

    // Circular speed in the *softened* potential, so the orbit is truly circular
    // for the force law we integrate.
    const double d2 = r * r + eps2_;
    const double v_circ = std::sqrt(kG * c.mass * r * r / (d2 * std::sqrt(d2)));

    const std::size_t k = base + i;
    stars_.x[k] = c.pos.x + r * co;
    stars_.y[k] = c.pos.y + r * s;
    stars_.vx[k] = c.vel.x - spin * v_circ * s;
    stars_.vy[k] = c.vel.y + spin * v_circ * co;
  }
}

// Shortest dynamical time in the setup, times step_accuracy. Considers the
// innermost star of each galaxy and, for two galaxies, the pericentre of the
// central pair's Kepler orbit (computed analytically from the initial state).
double Simulation::compute_dt() const {
  auto t_dyn = [](double r, double gm) { return std::sqrt(r * r * r / gm); };
  const double eps = params_.softening_ly;
  double t = std::numeric_limits<double>::infinity();

  auto consider_galaxy = [&](const GalaxyParams& g, const CentralBody& c) {
    if (g.num_stars == 0) return;
    const double r_in = std::max(g.inner_fraction * g.radius_ly, eps);
    t = std::min(t, t_dyn(r_in, kG * c.mass));
  };
  consider_galaxy(params_.galaxy_a, centrals_[0]);

  if (params_.galaxy_b) {
    consider_galaxy(*params_.galaxy_b, centrals_[1]);

    const Vec2 rp = centrals_[1].pos - centrals_[0].pos;
    const Vec2 rv = centrals_[1].vel - centrals_[0].vel;
    const double mu = kG * (centrals_[0].mass + centrals_[1].mass);
    const double r = rp.length();
    const double h = std::abs(cross(rp, rv));
    const double en = 0.5 * rv.length2() - mu / r;
    const double ecc = std::sqrt(std::max(0.0, 1.0 + 2.0 * en * h * h / (mu * mu)));
    const double r_peri = (h * h / mu) / (1.0 + ecc);
    t = std::min(t, t_dyn(std::max(r_peri, eps), mu));
  }

  // Zero stars, single galaxy: nothing constrains dt; pick something sane.
  if (!std::isfinite(t)) t = t_dyn(10.0 * eps, kG * centrals_[0].mass);
  return params_.step_accuracy * t;
}

double Simulation::compute_extent() const {
  const Vec2 com = center_of_mass();
  double m2 = 0.0;
  for (const auto& c : centrals_) m2 = std::max(m2, (c.pos - com).length2());
  for (std::size_t i = 0; i < stars_.size(); ++i) {
    const double dx = stars_.x[i] - com.x, dy = stars_.y[i] - com.y;
    m2 = std::max(m2, dx * dx + dy * dy);
  }
  return std::sqrt(m2);
}

// ---------------------------------------------------------------------------
// Integration (kick-drift-kick leapfrog)
// ---------------------------------------------------------------------------

void Simulation::compute_central_acc() {
  for (auto& c : centrals_) c.acc = {};
  if (centrals_.size() < 2) return;
  auto& a = centrals_[0];
  auto& b = centrals_[1];
  const Vec2 d = b.pos - a.pos;
  const double r2 = d.length2() + eps2_;
  const double inv = kG / (r2 * std::sqrt(r2));
  a.acc = d * (inv * b.mass);
  b.acc = d * (-inv * a.mass);
}

void Simulation::step_centrals(double dt) {
  const double half = 0.5 * dt;
  for (auto& c : centrals_) {
    c.vel += c.acc * half;
    c.pos += c.vel * dt;
  }
  compute_central_acc();
  for (auto& c : centrals_) c.vel += c.acc * half;
}

void Simulation::step(double dt) {
  // Stars never act back on the central bodies, so advance the centrals first.
  // Then every star needs just ONE pass: half-kick with the cached
  // acceleration, drift, recompute acceleration against the new central
  // positions, half-kick, cache. Equivalent to textbook KDK leapfrog.
  step_centrals(dt);

  double sx[2], sy[2], sgm[2];
  const int nc = static_cast<int>(centrals_.size());
  for (int c = 0; c < nc; ++c) {
    sx[c] = centrals_[c].pos.x;
    sy[c] = centrals_[c].pos.y;
    sgm[c] = kG * centrals_[c].mass;
  }

  const double half = 0.5 * dt;
  const double eps2 = eps2_;
  double* x = stars_.x.data();
  double* y = stars_.y.data();
  double* vx = stars_.vx.data();
  double* vy = stars_.vy.data();
  double* ax = stars_.ax.data();
  double* ay = stars_.ay.data();
  const auto n = static_cast<std::ptrdiff_t>(stars_.size());

#pragma omp parallel for if (n > 8192) schedule(static)
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    double vxi = vx[i] + half * ax[i];
    double vyi = vy[i] + half * ay[i];
    const double xi = x[i] + vxi * dt;
    const double yi = y[i] + vyi * dt;
    double axn, ayn;
    accel_from_sources(sx, sy, sgm, nc, eps2, xi, yi, axn, ayn);
    vxi += half * axn;
    vyi += half * ayn;
    x[i] = xi;
    y[i] = yi;
    vx[i] = vxi;
    vy[i] = vyi;
    ax[i] = axn;
    ay[i] = ayn;
  }

  time_ += dt;
}

double Simulation::run_queued(int max_steps) {
  auto n = static_cast<long long>(carry_ / dt_);
  if (n > max_steps) n = max_steps;
  for (long long i = 0; i < n; ++i) step(dt_);
  carry_ -= static_cast<double>(n) * dt_;
  return static_cast<double>(n) * dt_;
}

// Keep the sub-step remainder, but never let a backlog build up.
void Simulation::drop_backlog() { carry_ = std::min(carry_, dt_); }

void Simulation::fill_render_points(Vec2 view_center, double lag, float* xy) const {
  const std::size_t n = stars_.size();
  for (std::size_t i = 0; i < n; ++i) {
    xy[2 * i] = static_cast<float>(stars_.x[i] + stars_.vx[i] * lag - view_center.x);
    xy[2 * i + 1] = static_cast<float>(stars_.y[i] + stars_.vy[i] * lag - view_center.y);
  }
  for (std::size_t c = 0; c < centrals_.size(); ++c) {
    const CentralBody& b = centrals_[c];
    xy[2 * (n + c)] = static_cast<float>(b.pos.x + b.vel.x * lag - view_center.x);
    xy[2 * (n + c) + 1] = static_cast<float>(b.pos.y + b.vel.y * lag - view_center.y);
  }
}

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

double Simulation::total_energy() const {
  double e = 0.0;

  // Central bodies: kinetic + mutual (softened) potential.
  for (const auto& c : centrals_) e += 0.5 * c.mass * c.vel.length2();
  if (centrals_.size() == 2) {
    const double d2 = (centrals_[1].pos - centrals_[0].pos).length2() + eps2_;
    e -= kG * centrals_[0].mass * centrals_[1].mass / std::sqrt(d2);
  }

  // Stars: kinetic + potential in the field of the central bodies.
  double es = 0.0;
  for (std::size_t i = 0; i < stars_.size(); ++i) {
    double pot = 0.0;
    for (const auto& c : centrals_) {
      const double dx = stars_.x[i] - c.pos.x, dy = stars_.y[i] - c.pos.y;
      pot -= kG * c.mass / std::sqrt(dx * dx + dy * dy + eps2_);
    }
    es += kStarMassMsun * (0.5 * (stars_.vx[i] * stars_.vx[i] + stars_.vy[i] * stars_.vy[i]) + pot);
  }
  return e + es;
}

// Instantaneous energy test. In an encounter the potential is time-dependent,
// so "unbound" here means "unbound right now", not a guarantee of escape.
BoundFractions Simulation::bound_fractions() const {
  const std::size_t total = stars_.size();
  if (total == 0) return {};

  const int nc = static_cast<int>(centrals_.size());
  double cx[2], cy[2], cvx[2], cvy[2], cgm[2];
  for (int c = 0; c < nc; ++c) {
    cx[c] = centrals_[c].pos.x;
    cy[c] = centrals_[c].pos.y;
    cvx[c] = centrals_[c].vel.x;
    cvy[c] = centrals_[c].vel.y;
    cgm[c] = kG * centrals_[c].mass;
  }
  const Vec2 vcom = center_of_mass_velocity();
  const double eps2 = eps2_;
  const double* x = stars_.x.data();
  const double* y = stars_.y.data();
  const double* vx = stars_.vx.data();
  const double* vy = stars_.vy.data();
  const auto n = static_cast<std::ptrdiff_t>(total);

  long long n_system = 0, n_galaxy = 0;
#pragma omp parallel for reduction(+ : n_system, n_galaxy) if (n > 8192) schedule(static)
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    double phi = 0.0;  // combined potential at the star
    bool bound_to_a_galaxy = false;
    for (int c = 0; c < nc; ++c) {
      const double dx = x[i] - cx[c], dy = y[i] - cy[c];
      const double pot = cgm[c] / std::sqrt(dx * dx + dy * dy + eps2);
      phi -= pot;
      const double rvx = vx[i] - cvx[c], rvy = vy[i] - cvy[c];
      if (0.5 * (rvx * rvx + rvy * rvy) < pot) bound_to_a_galaxy = true;
    }
    const double svx = vx[i] - vcom.x, svy = vy[i] - vcom.y;
    if (0.5 * (svx * svx + svy * svy) + phi < 0.0) ++n_system;
    if (bound_to_a_galaxy) ++n_galaxy;
  }
  const auto dn = static_cast<double>(total);
  return {static_cast<double>(n_system) / dn, static_cast<double>(n_galaxy) / dn};
}

double Simulation::central_separation_ly() const {
  if (centrals_.size() < 2) return 0.0;
  return (centrals_[1].pos - centrals_[0].pos).length();
}

double Simulation::angular_momentum() const {
  double l = 0.0;
  for (const auto& c : centrals_) l += c.mass * cross(c.pos, c.vel);
  for (std::size_t i = 0; i < stars_.size(); ++i)
    l += kStarMassMsun * (stars_.x[i] * stars_.vy[i] - stars_.y[i] * stars_.vx[i]);
  return l;
}

Vec2 Simulation::center_of_mass() const {
  Vec2 s;
  double m = 0.0;
  for (const auto& c : centrals_) {
    s += c.pos * c.mass;
    m += c.mass;
  }
  return s * (1.0 / m);
}

Vec2 Simulation::center_of_mass_velocity() const {
  Vec2 s;
  double m = 0.0;
  for (const auto& c : centrals_) {
    s += c.vel * c.mass;
    m += c.mass;
  }
  return s * (1.0 / m);
}

}  // namespace galaxy
