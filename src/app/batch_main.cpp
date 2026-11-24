// Headless two-galaxy runner: no window, no GL. Runs the simulation for a fixed
// span of simulated time and writes the same CSV the GUI logger writes.
//
//   galaxy_batch --impact-frac 0.2 --gyr 5 --out logs/run.csv
//   galaxy_batch --impact 30000 --vel 100 --gyr 5 --out logs/run.csv
//
// Options:
//   --impact-frac <f>   impact parameter as a fraction of the starting distance
//   --impact <ly>       impact parameter in light-years (overrides --impact-frac)
//   --vel <km/s>        relative speed (default: the library default)
//   --gyr <Gyr>         simulated duration (default 5)
//   --out <path>        CSV output path (required)
//   --distributed       distributed-mass mode (self-gravitating, Barnes-Hut)
//   --stars <N>         particles / stars per galaxy
//   --bh-frac <f>       black hole mass fraction (distributed mode)
//   --halo-frac <f>     dark-matter halo mass fraction, 0 = no halo (distributed mode)
//   --halo-ratio <f>    halo particles per disk particle (distributed mode)
//   --profile <name>    uniform | exp (distributed mode)
//   --disp <f>          velocity dispersion fraction (distributed mode)
//   --theta <f>         Barnes-Hut opening angle (distributed mode)
//   --hole <f>          central hole radius as a fraction of the galaxy radius (both modes)
//   --relax <Myr>       pre-relax each galaxy alone first, 0 = off (default 800)
//   --one               a single galaxy (no encounter)
//   --snapshot <path>   at the end, write every regular body's "x,y" to a CSV

#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "app/csv_logger.hpp"
#include "sim/isimulation.hpp"

using namespace galaxy;

namespace {
constexpr double kLogIntervalYears = 1.0e6;  // one row per Myr, same as the GUI
constexpr double kProgressYears = 500.0e6;
}  // namespace

int main(int argc, char** argv) {
  ScenarioParams p;
  p.galaxy_b = GalaxyParams{};
  p.distributed.relax_myr = kDefaultRelaxMyr;

  std::optional<double> impact_frac, impact_ly;
  double gyr = 5.0;
  std::string out, snapshot;
  bool one_galaxy = false;

  for (int i = 1; i < argc; ++i) {
    const bool has_value = i + 1 < argc;
    if (!std::strcmp(argv[i], "--impact-frac") && has_value) impact_frac = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--impact") && has_value) impact_ly = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--vel") && has_value) p.encounter.relative_speed_km_s = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--gyr") && has_value) gyr = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--out") && has_value) out = argv[++i];
    else if (!std::strcmp(argv[i], "--one")) one_galaxy = true;
    else if (!std::strcmp(argv[i], "--snapshot") && has_value) snapshot = argv[++i];
    else if (!std::strcmp(argv[i], "--distributed")) p.mass_model = MassModel::Distributed;
    else if (!std::strcmp(argv[i], "--stars") && has_value) {
      p.galaxy_a.num_stars = p.galaxy_b->num_stars = std::atoi(argv[++i]);
    }
    else if (!std::strcmp(argv[i], "--bh-frac") && has_value) p.distributed.bh_mass_fraction = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--halo-frac") && has_value) p.distributed.halo_mass_fraction = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--halo-ratio") && has_value) p.distributed.halo_particles_per_disk_particle = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--disp") && has_value) p.distributed.dispersion_fraction = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--theta") && has_value) p.distributed.opening_angle = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--hole") && has_value) {
      p.distributed.inner_hole_fraction = std::atof(argv[++i]);
      p.galaxy_a.inner_fraction = p.galaxy_b->inner_fraction = p.distributed.inner_hole_fraction;
    }
    else if (!std::strcmp(argv[i], "--relax") && has_value) p.distributed.relax_myr = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--profile") && has_value) {
      const char* name = argv[++i];
      p.distributed.profile = !std::strcmp(name, "exp") ? DiskProfile::Exponential : DiskProfile::Uniform;
    }
    else {
      std::fprintf(stderr, "unknown or incomplete option: %s\n", argv[i]);
      return 2;
    }
  }
  if (one_galaxy) p.galaxy_b.reset();
  if (impact_ly) p.encounter.impact_parameter_ly = *impact_ly;
  else if (impact_frac) p.encounter.impact_parameter_ly = *impact_frac * p.encounter.separation_ly;
  if (out.empty() || !(gyr > 0.0)) {
    std::fprintf(stderr, "usage: galaxy_batch [--impact-frac f | --impact ly] [--vel km/s] [--gyr G] --out file.csv\n");
    return 2;
  }
  if (const std::string err = validate_scenario(p); !err.empty()) {
    std::fprintf(stderr, "invalid scenario: %s\n", err.c_str());
    return 2;
  }

  const std::unique_ptr<ISimulation> sim_owner = make_simulation(p);
  ISimulation& sim = *sim_owner;
  CsvLogger log;
  if (std::string err; !log.open(out, err)) {
    std::fprintf(stderr, "%s\n", err.c_str());
    return 1;
  }

  std::printf("%s: impact %.0f ly (%.1f%% of %.0f ly), vel %.0f km/s, %zu bodies, dt %.0f yr, %.2f Gyr -> %s\n",
              p.mass_model == MassModel::Distributed ? "distributed" : "centralized",
              p.encounter.impact_parameter_ly,
              100.0 * p.encounter.impact_parameter_ly / p.encounter.separation_ly,
              p.encounter.separation_ly, p.encounter.relative_speed_km_s,
              sim.render_counts().regular, sim.dt(), gyr, out.c_str());

  auto write_row = [&]() {
    const BoundFractions b = sim.bound_fractions();
    log.write_row(sim.time_years() / 1.0e6, sim.central_separation_ly(), 100.0 * b.system,
                  100.0 * b.galaxy);
  };

  const auto wall0 = std::chrono::steady_clock::now();
  const double end_years = gyr * 1.0e9;
  double next_progress = kProgressYears;
  write_row();  // t = 0

  while (sim.time_years() < end_years) {
    sim.queue_time(kLogIntervalYears);
    sim.run_queued(INT_MAX);  // sub-step remainder stays queued for the next row
    write_row();
    if (sim.time_years() >= next_progress) {
      const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
      std::printf("  t = %6.0f Myr  E/E0 - 1 = %+.2e  (%.0f s)\n", sim.time_years() / 1.0e6,
                  sim.energy_ratio() - 1.0, wall);
      std::fflush(stdout);
      next_progress += kProgressYears;
    }
  }

  if (!snapshot.empty()) {
    const RenderCounts rc = sim.render_counts();
    std::vector<float> xy(2 * (rc.regular + rc.highlighted + rc.halo));  // all groups are written
    sim.fill_render_points(Vec2{}, 0.0, xy.data());
    if (std::FILE* f = std::fopen(snapshot.c_str(), "w")) {
      std::fputs("x,y\n", f);
      for (std::size_t i = 0; i < rc.regular; ++i) std::fprintf(f, "%.1f,%.1f\n", xy[2 * i], xy[2 * i + 1]);
      std::fclose(f);
    } else {
      std::fprintf(stderr, "cannot write snapshot %s\n", snapshot.c_str());
    }
  }

  const BoundFractions b = sim.bound_fractions();
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
  std::printf("done: t = %.0f Myr, E/E0 - 1 = %+.3e, bound system %.2f%%, bound galaxy %.2f%%, wall %.0f s\n",
              sim.time_years() / 1.0e6, sim.energy_ratio() - 1.0, 100.0 * b.system,
              100.0 * b.galaxy, wall);
  return 0;
}
