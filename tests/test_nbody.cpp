// Tests for the Barnes-Hut tree and the distributed-mass (self-gravitating) mode.

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

#include "sim/nbody_simulation.hpp"
#include "sim/quadtree.hpp"

using namespace galaxy;

static int g_failures = 0;

#define CHECK(cond)                                                   \
  do {                                                                \
    if (!(cond)) {                                                    \
      std::printf("    FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      ++g_failures;                                                   \
    }                                                                 \
  } while (0)

#define CHECK_THROWS(expr)                                                         \
  do {                                                                             \
    bool threw_ = false;                                                           \
    try { (void)(expr); } catch (const std::invalid_argument&) { threw_ = true; }  \
    if (!threw_) {                                                                 \
      std::printf("    FAIL %s:%d: expected throw: %s\n", __FILE__, __LINE__, #expr); \
      ++g_failures;                                                                \
    }                                                                              \
  } while (0)

static ScenarioParams distributed(int n, std::uint64_t seed = 1) {
  ScenarioParams p;
  p.mass_model = MassModel::Distributed;
  p.galaxy_a.num_stars = n;
  p.seed = seed;
  return p;
}

static ScenarioParams distributed_two(int n) {
  ScenarioParams p = distributed(n);
  p.galaxy_b = GalaxyParams{};
  p.galaxy_b->num_stars = n;
  return p;
}

// --- tree ------------------------------------------------------------------

struct Cloud {
  std::vector<double> x, y, m;
};

static Cloud random_cloud(int n, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> pos(-1.0e4, 1.0e4), mass(1.0e6, 1.0e8);
  Cloud c;
  for (int i = 0; i < n; ++i) {
    // A dense core plus a halo, so the tree has both deep and shallow regions.
    const double scale = (i % 3 == 0) ? 0.05 : 1.0;
    c.x.push_back(pos(rng) * scale);
    c.y.push_back(pos(rng) * scale);
    c.m.push_back(mass(rng));
  }
  return c;
}

// Relative RMS force error of the tree against direct summation.
static double tree_force_error(const Cloud& c, double theta, double eps2) {
  const int n = static_cast<int>(c.x.size());
  BarnesHutTree tree;
  tree.build(c.x.data(), c.y.data(), c.m.data(), c.x.size());
  double err2 = 0.0, ref2 = 0.0;
  for (int i = 0; i < n; ++i) {
    double ax = 0, ay = 0;
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      const double dx = c.x[j] - c.x[i], dy = c.y[j] - c.y[i];
      const double d2 = dx * dx + dy * dy + eps2;
      const double inv = units::kG * c.m[j] / (d2 * std::sqrt(d2));
      ax += dx * inv;
      ay += dy * inv;
    }
    double tx, ty;
    tree.accel(i, c.x[i], c.y[i], theta, eps2, tx, ty);
    err2 += (tx - ax) * (tx - ax) + (ty - ay) * (ty - ay);
    ref2 += ax * ax + ay * ay;
  }
  return std::sqrt(err2 / ref2);
}

static void test_tree_matches_direct_sum() {
  const Cloud c = random_cloud(1500, 42);
  const double eps2 = 25.0 * 25.0;
  const double e_tiny = tree_force_error(c, 1e-6, eps2);  // opens every cell: exact
  const double e025 = tree_force_error(c, 0.25, eps2);
  const double e05 = tree_force_error(c, 0.5, eps2);
  const double e10 = tree_force_error(c, 1.0, eps2);
  std::printf("    rms relative force error: theta~0: %.2e, 0.25: %.2e, 0.5: %.2e, 1.0: %.2e\n",
              e_tiny, e025, e05, e10);
  CHECK(e_tiny < 1e-12);
  CHECK(e025 < e05 && e05 < e10);  // accuracy improves as theta shrinks
  CHECK(e05 < 0.02);               // under 2% rms at the default opening angle
}

static void test_tree_potential_matches_direct_sum() {
  const Cloud c = random_cloud(800, 7);
  const double eps2 = 400.0;
  BarnesHutTree tree;
  tree.build(c.x.data(), c.y.data(), c.m.data(), c.x.size());
  double worst = 0.0;
  for (int i = 0; i < 800; i += 37) {
    double direct = 0.0;
    for (int j = 0; j < 800; ++j) {
      if (j == i) continue;
      const double dx = c.x[j] - c.x[i], dy = c.y[j] - c.y[i];
      direct -= units::kG * c.m[j] / std::sqrt(dx * dx + dy * dy + eps2);
    }
    const double exact = tree.potential(i, c.x[i], c.y[i], 1e-6, eps2);
    const double approx = tree.potential(i, c.x[i], c.y[i], 0.5, eps2);
    worst = std::max(worst, std::abs(approx / direct - 1.0));
    CHECK(std::abs(exact / direct - 1.0) < 1e-12);
  }
  std::printf("    worst potential error at theta 0.5: %.2e\n", worst);
  CHECK(worst < 0.01);
}

static void test_tree_degenerate_inputs() {
  // Coincident bodies must not recurse forever; they act as one merged mass.
  std::vector<double> x = {5.0, 5.0, 5.0, 100.0}, y = {5.0, 5.0, 5.0, 100.0}, m = {1, 2, 3, 4};
  BarnesHutTree tree;
  tree.build(x.data(), y.data(), m.data(), x.size());
  double ax, ay;
  tree.accel(3, 100.0, 100.0, 0.0001, 1.0, ax, ay);  // must open down to the merged leaf
  CHECK(std::isfinite(ax) && std::isfinite(ay));
  CHECK(ax < 0.0 && ay < 0.0);  // pulled toward the cluster at (5, 5)

  // Empty and single-body trees.
  BarnesHutTree empty;
  empty.build(nullptr, nullptr, nullptr, 0);
  empty.accel(-1, 0.0, 0.0, 0.5, 1.0, ax, ay);
  CHECK(ax == 0.0 && ay == 0.0);
  double one_x = 1.0, one_y = 2.0, one_m = 3.0;
  BarnesHutTree one;
  one.build(&one_x, &one_y, &one_m, 1);
  one.accel(0, one_x, one_y, 0.5, 1.0, ax, ay);  // no force from itself
  CHECK(ax == 0.0 && ay == 0.0);
}

static void test_tree_leaf_order_is_a_permutation() {
  const Cloud c = random_cloud(2000, 3);
  BarnesHutTree tree;
  tree.build(c.x.data(), c.y.data(), c.m.data(), c.x.size());
  std::vector<std::int32_t> order = tree.leaf_order();
  std::sort(order.begin(), order.end());
  bool perm = order.size() == c.x.size();
  for (std::size_t i = 0; perm && i < order.size(); ++i) perm = order[i] == static_cast<std::int32_t>(i);
  CHECK(perm);
}

// --- initial conditions ----------------------------------------------------

static void test_mass_bookkeeping() {
  ScenarioParams p = distributed(2000);
  p.distributed.bh_mass_fraction = 0.02;
  p.distributed.halo_mass_fraction = 0.80;
  p.distributed.halo_particles_per_disk_particle = 2.0;
  NBodySimulation sim(p);
  CHECK(sim.particle_count() == 2000);
  CHECK(sim.halo_count() == 4000);
  CHECK(sim.black_hole_count() == 1);
  CHECK(sim.body_count() == 6001);

  double m_disk = 0.0, m_halo = 0.0;
  for (std::size_t i = 0; i < sim.particle_count(); ++i) m_disk += sim.mass()[i];
  for (std::size_t i = sim.particle_count(); i < sim.particle_count() + sim.halo_count(); ++i)
    m_halo += sim.mass()[i];
  const double m_bh = sim.mass().back();
  std::printf("    BH %.3e, disk %.3e (%.3e each), halo %.3e (%.3e each) Msun\n", m_bh, m_disk,
              sim.mass()[0], m_halo, sim.mass()[sim.particle_count()]);
  CHECK(std::abs(m_bh / 1.0e11 - 0.02) < 1e-12);
  CHECK(std::abs(m_halo / 1.0e11 - 0.80) < 1e-12);
  CHECK(std::abs(m_disk / 1.0e11 - 0.18) < 1e-12);

  const RenderCounts rc = sim.render_counts();
  CHECK(rc.regular == 2000 && rc.highlighted == 1 && rc.halo == 4000);

  // No halo: the disk and black hole only.
  p.distributed.halo_mass_fraction = 0.0;
  NBodySimulation no_halo(p);
  CHECK(no_halo.halo_count() == 0);
  CHECK(no_halo.body_count() == 2001);

  // No halo and no black hole: particles only.
  p.distributed.bh_mass_fraction = 0.0;
  NBodySimulation bare(p);
  CHECK(bare.black_hole_count() == 0);
  CHECK(bare.body_count() == 2000);
}

static void test_profiles() {
  // Uniform: <r^2> = R^2 / 2.
  {
    NBodySimulation sim(distributed(40000));
    double s = 0.0;
    for (std::size_t i = 0; i < sim.particle_count(); ++i)
      s += sim.x()[i] * sim.x()[i] + sim.y()[i] * sim.y()[i];
    const double got = s / sim.particle_count(), expect = 0.5 * 5.0e4 * 5.0e4;
    std::printf("    uniform <r^2> = %.4e (expected %.4e)\n", got, expect);
    CHECK(std::abs(got / expect - 1.0) < 0.03);
  }
  // Exponential (h = 0.2 R, truncated at R = 5 h): the fraction inside r < h is
  // (1 - 2/e) / (1 - 6 e^-5).
  {
    ScenarioParams p = distributed(40000);
    p.distributed.profile = DiskProfile::Exponential;
    NBodySimulation sim(p);
    const double h = 0.2 * 5.0e4;
    std::size_t inside = 0;
    for (std::size_t i = 0; i < sim.particle_count(); ++i)
      if (std::hypot(sim.x()[i], sim.y()[i]) < h) ++inside;
    const double got = static_cast<double>(inside) / sim.particle_count();
    const double expect = (1.0 - 2.0 / std::exp(1.0)) / (1.0 - 6.0 * std::exp(-5.0));
    std::printf("    exponential fraction within 1 scale length = %.4f (expected %.4f)\n", got,
                expect);
    CHECK(std::abs(got - expect) < 0.01);
  }
}

static void test_halo_structure() {
  NBodySimulation sim(distributed(3000));  // defaults: halo 90%, radius 2R, core 0.25R, 3 per disk
  const std::size_t d0 = sim.particle_count(), nh = sim.halo_count();
  CHECK(nh == 9000);

  const double R = 5.0e4, rh = 2.0 * R, a = 0.25 * R;
  const double K = std::sqrt(rh * rh + a * a) - a;
  double rmax = 0.0, px = 0.0, py = 0.0, lz = 0.0, scale = 0.0;
  std::size_t inside = 0;
  for (std::size_t i = d0; i < d0 + nh; ++i) {
    const double r = std::hypot(sim.x()[i], sim.y()[i]);
    rmax = std::max(rmax, r);
    if (r < R) ++inside;
    px += sim.vx()[i];
    py += sim.vy()[i];
    lz += sim.x()[i] * sim.vy()[i] - sim.y()[i] * sim.vx()[i];
    scale += r * std::hypot(sim.vx()[i], sim.vy()[i]);
  }
  const double got = static_cast<double>(inside) / static_cast<double>(nh);
  const double expect = (std::sqrt(R * R + a * a) - a) / K;
  std::printf("    halo: mass fraction inside R = %.4f (expected %.4f), max radius %.0f, |L|/scale = %.3f\n",
              got, expect, rmax, std::abs(lz) / scale);
  CHECK(rmax <= rh * 1.02);
  CHECK(std::abs(got - expect) < 0.02);
  CHECK(std::hypot(px, py) < 1e-9);   // no net motion
  CHECK(std::abs(lz) < 0.03 * scale);  // no net rotation (random-walk level only)
}

static void test_initial_virial_ratio() {
  // A galaxy set up in equilibrium has 2K / |W| close to 1 (W: exact direct sum).
  for (bool exponential : {false, true}) {
    ScenarioParams p = distributed(1000);
    if (exponential) p.distributed.profile = DiskProfile::Exponential;
    NBodySimulation sim(p);
    double ke = 0.0;
    for (std::size_t i = 0; i < sim.body_count(); ++i)
      ke += 0.5 * sim.mass()[i] * (sim.vx()[i] * sim.vx()[i] + sim.vy()[i] * sim.vy()[i]);
    const double w = sim.total_energy_direct() - ke;
    const double ratio = 2.0 * ke / std::abs(w);
    std::printf("    %s disk: 2K/|W| = %.3f\n", exponential ? "exponential" : "uniform", ratio);
    CHECK(ratio > 0.7 && ratio < 1.3);
  }
}

static void test_cold_disk_is_circular_and_bound() {
  ScenarioParams p = distributed(3000);
  p.distributed.dispersion_fraction = 0.0;
  NBodySimulation sim(p);
  // With no dispersion every disk velocity is tangential, up to the small uniform
  // shift that removes the sampling noise in the mean velocity (~v/sqrt(N)).
  // That shift only matters where the rotation speed is itself small, so check
  // the outer half of the disk.
  double worst_radial = 0.0;
  for (std::size_t i = 0; i < sim.particle_count(); ++i) {
    const double r = std::hypot(sim.x()[i], sim.y()[i]);
    if (r < 0.5 * 5.0e4) continue;
    const double vr = (sim.x()[i] * sim.vx()[i] + sim.y()[i] * sim.vy()[i]) / r;
    const double vt = std::hypot(sim.vx()[i], sim.vy()[i]);
    worst_radial = std::max(worst_radial, std::abs(vr) / vt);
  }
  std::printf("    worst |v_radial| / |v| in the outer half: %.3e\n", worst_radial);
  CHECK(worst_radial < 0.05);
  CHECK(sim.bound_fractions().system > 0.99);
  CHECK(std::isnan(sim.bound_fractions().galaxy));
}

static void test_two_galaxy_frame() {
  NBodySimulation sim(distributed_two(1500));
  CHECK(sim.black_hole_count() == 2);
  CHECK(sim.center_of_mass().length() < 1e-6);
  CHECK(sim.total_momentum().length() / 2.0e11 < 1e-12);  // per unit mass, at rest
  const auto& e = sim.params().encounter;
  CHECK(std::abs(sim.central_separation_ly() / e.separation_ly - 1.0) < 1e-9);
  CHECK(sim.initial_extent() > e.separation_ly * 0.4);
  // The halo extends past the disks, so the dark-matter view is wider.
  CHECK(sim.dark_matter_extent() > sim.initial_extent());

  // Render order: disk particles, black holes, then halo.
  const RenderCounts rc = sim.render_counts();
  std::vector<float> xy(2 * (rc.regular + rc.highlighted + rc.halo));
  sim.fill_render_points(Vec2{}, 0.0, xy.data());
  const std::size_t bh0 = rc.regular;  // first black hole in the output
  CHECK(std::abs(xy[2 * bh0] - static_cast<float>(sim.x()[sim.body_count() - 2])) < 1.0f);
  const std::size_t halo0 = rc.regular + rc.highlighted;
  CHECK(std::abs(xy[2 * halo0] - static_cast<float>(sim.x()[sim.particle_count()])) < 1.0f);
}

// --- dynamics --------------------------------------------------------------

static void test_isolated_galaxy_conservation() {
  ScenarioParams p = distributed(1000);
  p.distributed.profile = DiskProfile::Exponential;
  NBodySimulation sim(p);
  const double e_direct0 = sim.total_energy_direct();
  const double l0 = sim.angular_momentum();
  const double e_tree0 = sim.initial_energy();
  sim.advance(3.0e8, INT_MAX);
  const double dE = sim.total_energy_direct() / e_direct0 - 1.0;
  const double dL = sim.angular_momentum() / l0 - 1.0;
  std::printf("    dt = %.0f yr, tree E0 vs direct E0: %.2e, 300 Myr: dE (direct) = %.2e, dL = %.2e\n",
              sim.dt(), e_tree0 / e_direct0 - 1.0, dE, dL);
  CHECK(std::abs(e_tree0 / e_direct0 - 1.0) < 0.02);  // tree PE within 2% of exact
  CHECK(std::abs(dE) < 0.02);
  CHECK(std::abs(dL) < 0.05);
}

// How lumpy the visible disk is: the second moment of particle counts in a 24x24
// grid over the cells inside the original disk radius, normalised so a uniform
// distribution (plus shot noise) gives ~1. Also the fraction of particles that
// have wandered beyond 1.5x the original radius.
struct DiskStats {
  double clumpiness = 0.0;
  double beyond = 0.0;
};

static DiskStats disk_stats(const NBodySimulation& sim) {
  constexpr double R = 5.0e4;
  constexpr int n = 24;
  const double span = 3.0 * R, cell = span / n;
  std::vector<double> counts(n * n, 0.0);
  std::size_t beyond = 0;
  for (std::size_t i = 0; i < sim.particle_count(); ++i) {
    const double x = sim.x()[i], y = sim.y()[i];
    if (std::hypot(x, y) > 1.5 * R) {
      ++beyond;
      continue;
    }
    const int cx = static_cast<int>((x + 1.5 * R) / cell), cy = static_cast<int>((y + 1.5 * R) / cell);
    if (cx >= 0 && cx < n && cy >= 0 && cy < n) counts[cy * n + cx] += 1.0;
  }
  double s1 = 0.0, s2 = 0.0;
  int cells = 0;
  for (int cy = 0; cy < n; ++cy)
    for (int cx = 0; cx < n; ++cx) {
      const double mx = -1.5 * R + (cx + 0.5) * cell, my = -1.5 * R + (cy + 0.5) * cell;
      if (std::hypot(mx, my) >= R) continue;  // only cells inside the disk radius
      const double c = counts[cy * n + cx];
      s1 += c;
      s2 += c * c;
      ++cells;
    }
  DiskStats out;
  out.clumpiness = cells * s2 / (s1 * s1);
  out.beyond = static_cast<double>(beyond) / static_cast<double>(sim.particle_count());
  return out;
}

static void test_halo_keeps_the_disk_from_fragmenting() {
  // The reason the halo exists: a disk holding ~all of a galaxy's mass has
  // Toomre Q ~ 0.1 and tears itself into clumps. Same disk, with and without
  // the halo, over ~5 dynamical times.
  auto run = [](double halo_fraction) {
    ScenarioParams p = distributed(1200);
    p.distributed.halo_mass_fraction = halo_fraction;
    NBodySimulation sim(p);
    const DiskStats before = disk_stats(sim);
    sim.advance(5.0e8, INT_MAX);
    return std::make_pair(before, disk_stats(sim));
  };
  const auto with_halo = run(0.90);
  const auto no_halo = run(0.0);
  std::printf("    with halo: clumpiness %.2f -> %.2f, %.1f%% beyond 1.5R\n", with_halo.first.clumpiness,
              with_halo.second.clumpiness, 100.0 * with_halo.second.beyond);
  std::printf("    no halo:   clumpiness %.2f -> %.2f, %.1f%% beyond 1.5R\n", no_halo.first.clumpiness,
              no_halo.second.clumpiness, 100.0 * no_halo.second.beyond);
  CHECK(with_halo.second.clumpiness < 2.4);  // still a smooth disk (no halo: ~11)
  CHECK(with_halo.second.beyond < 0.10);     // and still mostly where it started
  CHECK(no_halo.second.beyond > 0.25 || no_halo.second.clumpiness > 2.0 * with_halo.second.clumpiness);
}

// Median disk radius swing over `watch_myr` after an optional pre-relaxation.
static double disk_radius_swing(double relax_myr, double watch_myr, double* r0_out = nullptr) {
  ScenarioParams p = distributed(1000, 7);
  p.distributed.relax_myr = relax_myr;
  NBodySimulation sim(p);
  auto median_r = [&] {
    std::vector<double> r;
    for (std::size_t i = 0; i < sim.particle_count(); ++i) r.push_back(std::hypot(sim.x()[i], sim.y()[i]));
    std::nth_element(r.begin(), r.begin() + static_cast<std::ptrdiff_t>(r.size() / 2), r.end());
    return r[r.size() / 2];
  };
  const double r0 = median_r();
  if (r0_out) *r0_out = r0;
  double worst = 0.0;
  for (int k = 0; k < 8; ++k) {
    sim.advance(watch_myr / 8.0 * 1.0e6, INT_MAX);
    worst = std::max(worst, std::abs(median_r() / r0 - 1.0));
  }
  return worst;
}

static void test_pre_relaxation_removes_the_breathing() {
  // A freshly sampled galaxy is not an equilibrium: the whole thing collapses
  // and rebounds (median radius -40%). After pre-relaxing it, it stays put.
  const double raw = disk_radius_swing(0.0, 400.0);
  const double relaxed = disk_radius_swing(kDefaultRelaxMyr, 400.0);
  std::printf("    median-radius swing over 400 Myr: raw %.1f%%, pre-relaxed %.1f%%\n", 100.0 * raw,
              100.0 * relaxed);
  CHECK(raw > 0.15);
  CHECK(relaxed < 0.5 * raw);
  CHECK(relaxed < 0.2);

  // Relaxation keeps the layout, masses and centre of the galaxy.
  ScenarioParams p = distributed_two(500);
  p.distributed.relax_myr = 100.0;
  NBodySimulation sim(p);
  CHECK(sim.particle_count() == 1000 && sim.halo_count() == 3000 && sim.black_hole_count() == 2);
  CHECK(std::abs(sim.central_separation_ly() / p.encounter.separation_ly - 1.0) < 0.05);
}

static void test_central_hole() {
  // The hole empties the centre of disk AND halo, keeps counts and masses, and
  // allows a much larger timestep (the fastest motion sets dt).
  ScenarioParams base = distributed(2000, 5);
  NBodySimulation full(base);
  ScenarioParams hp = base;
  hp.distributed.inner_hole_fraction = 0.2;
  NBodySimulation holed(hp);

  const double r_hole = 0.2 * hp.galaxy_a.radius_ly;
  double rmin = 1e30;
  for (std::size_t i = 0; i + holed.black_hole_count() < holed.body_count(); ++i)
    rmin = std::min(rmin, std::hypot(holed.x()[i], holed.y()[i]));
  std::printf("    hole radius %.0f ly, innermost disk/halo particle %.0f ly, dt %.0f -> %.0f yr\n",
              r_hole, rmin, full.dt(), holed.dt());
  CHECK(rmin > 0.85 * r_hole);  // re-centring on the sampled mean shifts the hole by ~sampling noise
  CHECK(holed.particle_count() == full.particle_count());
  CHECK(holed.halo_count() == full.halo_count());
  CHECK(holed.dt() > 2.0 * full.dt());

  // Exponential disk too.
  hp.distributed.profile = DiskProfile::Exponential;
  NBodySimulation expo(hp);
  double rmin_e = 1e30;
  for (std::size_t i = 0; i < expo.particle_count(); ++i)
    rmin_e = std::min(rmin_e, std::hypot(expo.x()[i], expo.y()[i]));
  CHECK(rmin_e > 0.85 * r_hole);
  hp.distributed.inner_hole_fraction = 0.95;
  CHECK_THROWS(NBodySimulation(hp));
}

static void test_two_galaxy_run_is_sane() {
  ScenarioParams p = distributed_two(500);
  p.encounter.impact_parameter_ly = 60000.0;
  NBodySimulation sim(p);
  const double sep0 = sim.central_separation_ly();
  sim.advance(2.0e8, INT_MAX);
  const double sep1 = sim.central_separation_ly();
  const auto b = sim.bound_fractions();
  std::printf("    separation %.0f -> %.0f ly after 200 Myr, bound %.3f, E/E0 - 1 = %.2e\n", sep0,
              sep1, b.system, sim.energy_ratio() - 1.0);
  CHECK(sep1 < sep0);  // they fall toward each other
  CHECK(b.system >= 0.0 && b.system <= 1.0);
  CHECK(sim.center_of_mass().length() < 100.0);  // COM stays put (tree momentum error is small)
}

static void test_validation() {
  ScenarioParams p = distributed(100);
  p.distributed.bh_mass_fraction = 0.9;
  CHECK_THROWS(NBodySimulation{p});
  p = distributed(100);
  p.distributed.opening_angle = 0.0;
  CHECK_THROWS(NBodySimulation{p});
  p = distributed(100);
  p.distributed.dispersion_fraction = -0.1;
  CHECK_THROWS(NBodySimulation{p});
  p = distributed(100);
  p.distributed.halo_mass_fraction = 0.99;  // above the 0.95 cap
  CHECK_THROWS(NBodySimulation{p});
  p = distributed(100);
  p.distributed.bh_mass_fraction = 0.2;
  p.distributed.halo_mass_fraction = 0.9;  // BH + halo leave < 2% for the disk
  CHECK_THROWS(NBodySimulation{p});
  p = distributed(100);
  p.distributed.halo_core_fraction = 0.0;
  CHECK_THROWS(NBodySimulation{p});
  p = distributed(0);
  CHECK_THROWS(NBodySimulation{p});
  p = ScenarioParams{};  // centralized scenario handed to the wrong backend
  CHECK_THROWS(NBodySimulation{p});
}

#define RUN(test)                                                              \
  do {                                                                         \
    const int before = g_failures;                                             \
    std::printf("[ RUN  ] %s\n", #test);                                       \
    test();                                                                    \
    std::printf("[ %s ] %s\n", g_failures == before ? " OK " : "FAIL", #test); \
  } while (0)

int main() {
  RUN(test_tree_matches_direct_sum);
  RUN(test_tree_potential_matches_direct_sum);
  RUN(test_tree_degenerate_inputs);
  RUN(test_tree_leaf_order_is_a_permutation);
  RUN(test_mass_bookkeeping);
  RUN(test_profiles);
  RUN(test_halo_structure);
  RUN(test_initial_virial_ratio);
  RUN(test_cold_disk_is_circular_and_bound);
  RUN(test_two_galaxy_frame);
  RUN(test_isolated_galaxy_conservation);
  RUN(test_halo_keeps_the_disk_from_fragmenting);
  RUN(test_pre_relaxation_removes_the_breathing);
  RUN(test_central_hole);
  RUN(test_two_galaxy_run_is_sane);
  RUN(test_validation);
  std::printf("\n%s (%d failed check%s)\n", g_failures ? "FAILED" : "ALL PASSED", g_failures,
              g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
