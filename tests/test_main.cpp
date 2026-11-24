// Minimal dependency-free test runner. Each test prints its key numbers so the
// output doubles as a record of how well the integrator is behaving.

#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>

#include "sim/simulation.hpp"

using namespace galaxy;

static int g_failures = 0;

#define CHECK(cond)                                                        \
  do {                                                                     \
    if (!(cond)) {                                                         \
      std::printf("    FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);      \
      ++g_failures;                                                        \
    }                                                                      \
  } while (0)

#define CHECK_THROWS(expr)                                                 \
  do {                                                                     \
    bool threw_ = false;                                                   \
    try { (void)(expr); } catch (const std::invalid_argument&) { threw_ = true; } \
    if (!threw_) {                                                         \
      std::printf("    FAIL %s:%d: expected throw: %s\n", __FILE__, __LINE__, #expr); \
      ++g_failures;                                                        \
    }                                                                      \
  } while (0)

static ScenarioParams single_galaxy(int n, std::uint64_t seed = 1) {
  ScenarioParams p;
  p.galaxy_a.num_stars = n;
  p.seed = seed;
  return p;
}

static ScenarioParams two_galaxies(int n) {
  ScenarioParams p = single_galaxy(n);
  p.galaxy_b = GalaxyParams{};
  p.galaxy_b->num_stars = n;
  return p;
}

// --- tests -----------------------------------------------------------------

static void test_units() {
  // Earth: 1 AU around 1 Msun should take ~1 Julian year.
  const double au_ly = 1.495978707e11 / units::kLightYear_m;
  const double period = 2.0 * std::numbers::pi * std::sqrt(au_ly * au_ly * au_ly / units::kG);
  std::printf("    Earth period from G: %.6f yr, G = %.5e ly^3/(Msun yr^2)\n", period, units::kG);
  CHECK(std::abs(period - 1.0) < 1e-4);
}

static void test_distribution() {
  ScenarioParams p = single_galaxy(40000);
  Simulation sim(p);
  const double r_in = p.galaxy_a.inner_fraction * p.galaxy_a.radius_ly;
  const double r_out = p.galaxy_a.radius_ly;
  const auto& s = sim.stars();
  CHECK(s.size() == 40000);

  double sum_r2 = 0.0, rmin = 1e300, rmax = 0.0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const double r = std::hypot(s.x[i], s.y[i]);
    sum_r2 += r * r;
    rmin = std::min(rmin, r);
    rmax = std::max(rmax, r);
  }
  // Uniform in area  =>  <r^2> = (r_in^2 + r_out^2) / 2.
  const double expect = 0.5 * (r_in * r_in + r_out * r_out);
  const double got = sum_r2 / static_cast<double>(s.size());
  std::printf("    <r^2> = %.4e (expected %.4e), r in [%.0f, %.0f]\n", got, expect, rmin, rmax);
  CHECK(std::abs(got / expect - 1.0) < 0.02);
  CHECK(rmin >= r_in * 0.999999 && rmax <= r_out * 1.000001);

  // Central mass = total - stars.
  CHECK(std::abs(sim.centrals()[0].mass - (1e11 - 40000.0)) < 1e-3);
  // Camera extent covers the whole disk.
  CHECK(sim.initial_extent() >= rmax - 1e-6);
}

static void test_determinism() {
  Simulation a(single_galaxy(1000, 7)), b(single_galaxy(1000, 7)), c(single_galaxy(1000, 8));
  bool same = true, differs = false;
  for (std::size_t i = 0; i < 1000; ++i) {
    same = same && a.stars().x[i] == b.stars().x[i] && a.stars().vy[i] == b.stars().vy[i];
    differs = differs || a.stars().x[i] != c.stars().x[i];
  }
  CHECK(same);
  CHECK(differs);
}

static void test_validation() {
  ScenarioParams p = single_galaxy(100);
  p.galaxy_a.total_mass_msun = 50.0;  // less than the stars' own mass
  CHECK_THROWS(Simulation{p});

  p = single_galaxy(100);
  p.galaxy_a.radius_ly = 0.0;
  CHECK_THROWS(Simulation{p});

  p = two_galaxies(100);
  p.encounter.impact_parameter_ly = p.encounter.separation_ly;  // not allowed
  CHECK_THROWS(Simulation{p});

  p = single_galaxy(100);
  p.softening_ly = 0.0;
  CHECK_THROWS(Simulation{p});
}

static void test_single_galaxy_conservation() {
  Simulation sim(single_galaxy(2000));
  const double L0 = sim.angular_momentum();
  const double years = 3.0e9;
  const double adv = sim.advance(years, INT_MAX);
  const double dE = sim.energy_ratio() - 1.0;
  const double dL = sim.angular_momentum() / L0 - 1.0;
  std::printf("    dt = %.1f yr, advanced %.3e yr, E/E0 - 1 = %.3e, L/L0 - 1 = %.3e\n",
              sim.dt(), adv, dE, dL);
  CHECK(std::abs(dE) < 1e-5);
  CHECK(std::abs(dL) < 1e-9);
  CHECK(std::abs(adv - years) <= sim.dt());
}

static void test_two_galaxy_conservation() {
  ScenarioParams p = two_galaxies(500);
  Simulation sim(p);
  // The COM frame must start at rest at the origin.
  CHECK(sim.center_of_mass().length() < 1e-6);
  CHECK(sim.center_of_mass_velocity().length() < 1e-15);

  const double years = 8.0e9;
  sim.advance(years, INT_MAX);
  const double dE = sim.energy_ratio() - 1.0;
  const double com = sim.center_of_mass().length();
  const double comv = sim.center_of_mass_velocity().length();
  const double sep = (sim.centrals()[1].pos - sim.centrals()[0].pos).length();
  std::printf("    dt = %.1f yr, E/E0 - 1 = %.3e, |COM| = %.3e ly, |V_com| = %.3e, final sep = %.0f ly\n",
              sim.dt(), dE, com, comv, sep);
  CHECK(std::abs(dE) < 1e-4);
  CHECK(com < 1e-3);  // momentum is conserved to round-off
  CHECK(comv < 1e-12);
}

static void test_second_order_convergence() {
  // A single star on a circular orbit; after exactly one period it should be
  // back where it started. Leapfrog's phase error scales as dt^2, so halving
  // dt must cut the error by ~4.
  auto error_after_one_period = [](int n_steps) {
    ScenarioParams p = single_galaxy(1);
    Simulation sim(p);
    const double x0 = sim.stars().x[0], y0 = sim.stars().y[0];
    const double r = std::hypot(x0, y0);
    const double d2 = r * r + p.softening_ly * p.softening_ly;
    const double omega = std::sqrt(units::kG * sim.centrals()[0].mass / (d2 * std::sqrt(d2)));
    const double period = 2.0 * std::numbers::pi / omega;
    const double dt = period / n_steps;
    for (int i = 0; i < n_steps; ++i) sim.step(dt);
    return std::hypot(sim.stars().x[0] - x0, sim.stars().y[0] - y0) / r;
  };
  const double e1 = error_after_one_period(100);
  const double e2 = error_after_one_period(200);
  const double e3 = error_after_one_period(400);
  std::printf("    rel. position error after 1 orbit: n=100: %.3e, n=200: %.3e, n=400: %.3e\n", e1, e2, e3);
  std::printf("    ratios: %.2f, %.2f (2nd order => 4)\n", e1 / e2, e2 / e3);
  CHECK(e1 / e2 > 3.5 && e1 / e2 < 4.5);
  CHECK(e2 / e3 > 3.5 && e2 / e3 < 4.5);
}

static void test_bound_fractions() {
  // A single galaxy in a static field: every star is bound, always.
  {
    Simulation sim(single_galaxy(2000));
    auto b0 = sim.bound_fractions();
    sim.advance(2.0e9, INT_MAX);
    auto b1 = sim.bound_fractions();
    std::printf("    single galaxy: system %.4f -> %.4f, galaxy %.4f -> %.4f\n", b0.system,
                b1.system, b0.galaxy, b1.galaxy);
    CHECK(b0.system == 1.0 && b0.galaxy == 1.0);
    CHECK(b1.system == 1.0 && b1.galaxy == 1.0);
  }
  // Two galaxies at rest, far apart: everything is bound both ways.
  {
    ScenarioParams p = two_galaxies(500);
    p.encounter.relative_speed_km_s = 0.0;
    p.encounter.impact_parameter_ly = 0.0;
    p.encounter.separation_ly = 600000.0;
    Simulation sim(p);
    auto b = sim.bound_fractions();
    CHECK(b.system == 1.0 && b.galaxy == 1.0);
  }
  // A very fast (unbound) pass: untouched stars are still bound to their own
  // galaxy, but the system-frame test already calls many of them unbound.
  {
    ScenarioParams p = two_galaxies(500);
    p.encounter.relative_speed_km_s = 800.0;
    Simulation sim(p);
    auto b = sim.bound_fractions();
    std::printf("    800 km/s pass at t=0: system %.3f, galaxy %.3f\n", b.system, b.galaxy);
    CHECK(b.galaxy == 1.0);
    CHECK(b.system < 1.0);
  }
  // Default flyby over a long run: fractions stay valid.
  {
    Simulation sim(two_galaxies(500));
    sim.advance(8.0e9, INT_MAX);
    auto b = sim.bound_fractions();
    std::printf("    default flyby after 8 Gyr: system %.3f, galaxy %.3f\n", b.system, b.galaxy);
    CHECK(b.system >= 0.0 && b.system <= 1.0);
    CHECK(b.galaxy >= 0.0 && b.galaxy <= 1.0);
  }
}

static void test_advance_cap() {
  Simulation sim(single_galaxy(100));
  const double dt = sim.dt();
  // Ask for 1000 steps' worth but allow only 10: must report honest progress.
  const double adv = sim.advance(1000.0 * dt, 10);
  CHECK(std::abs(adv - 10.0 * dt) < 1e-9 * dt);
  // No backlog should build up: a following small request advances ~0-1 steps.
  const double adv2 = sim.advance(0.5 * dt, 10);
  CHECK(adv2 <= dt * 1.0000001);
}

// --- runner ----------------------------------------------------------------

#define RUN(test)                                   \
  do {                                              \
    const int before = g_failures;                  \
    std::printf("[ RUN  ] %s\n", #test);            \
    test();                                         \
    std::printf("[ %s ] %s\n", g_failures == before ? " OK " : "FAIL", #test); \
  } while (0)

int main() {
  RUN(test_units);
  RUN(test_distribution);
  RUN(test_determinism);
  RUN(test_validation);
  RUN(test_single_galaxy_conservation);
  RUN(test_two_galaxy_conservation);
  RUN(test_second_order_convergence);
  RUN(test_bound_fractions);
  RUN(test_advance_cap);
  std::printf("\n%s (%d failed check%s)\n", g_failures ? "FAILED" : "ALL PASSED", g_failures,
              g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
