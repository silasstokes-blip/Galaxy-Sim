#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include "sim/units.hpp"

namespace galaxy {

// One galaxy: a single massive central body plus N one-solar-mass test stars
// on circular orbits around it.
struct GalaxyParams {
  // Total galaxy mass. The central body gets (total - num_stars * 1 Msun).
  double total_mass_msun = 1.0e11;
  int num_stars = 50000;
  // Outer radius of the star disk.
  double radius_ly = 50000.0;
  // Stars are distributed uniformly in AREA between inner_fraction*radius and
  // radius. A small hole in the middle keeps the fixed timestep affordable:
  // the innermost orbit sets dt for the whole sim.
  double inner_fraction = 0.05;
  // Spin direction of the disk (counter-clockwise when viewed with +y up).
  bool counterclockwise = true;
};

// Only used when two galaxies are present. Everything is set up in the
// centre-of-mass frame of the two central bodies.
struct EncounterParams {
  // Initial distance between the two central bodies.
  double separation_ly = 150000.0;
  // Initial relative speed (B relative to A), directed along -x.
  double relative_speed_km_s = 100.0;
  // Perpendicular offset of the approach line; 0 = head-on (radial) plunge.
  double impact_parameter_ly = 30000.0;
};

// How a galaxy's mass is represented.
enum class MassModel {
  // One massive central body per galaxy; stars are 1 Msun test particles that
  // feel it but exert no force. Two-body (Kepler) orbit for the pair. O(N).
  Centralized,
  // A central black hole plus N equal-mass macro-particles that carry the rest
  // of the galaxy's mass and all attract each other (Barnes-Hut tree).
  Distributed,
};

enum class DiskProfile { Uniform, Exponential };

// Used only when mass_model == Distributed.
struct DistributedParams {
  // Fraction of each galaxy's total mass in its central black hole.
  double bh_mass_fraction = 0.01;
  DiskProfile profile = DiskProfile::Uniform;
  // Exponential disk scale length as a fraction of the galaxy radius.
  double scale_length_fraction = 0.2;
  // Random velocity per component, as a fraction of the local circular speed.
  // 0 gives a cold disk (pure circular orbits). With the default halo, 20% keeps
  // the disk's Toomre Q above ~1 across the disk; a disk with no halo would need
  // 60-90%, which is why a halo exists.
  double dispersion_fraction = 0.20;

  // Live dark-matter halo: invisible particles that carry most of the mass and
  // keep the disk from fragmenting. They are in-plane (2D) like everything
  // else, with surface density ~ 1/sqrt(R^2 + core^2) out to the halo radius
  // (a roughly flat rotation curve), and a pressure-supported, non-rotating
  // velocity distribution from the isotropic Jeans equation. 0 disables the halo.
  double halo_mass_fraction = 0.90;          // of the galaxy's TOTAL mass
  double halo_radius_fraction = 2.0;         // truncation radius / galaxy radius
  double halo_core_fraction = 0.25;          // core radius / galaxy radius
  // Halo particles per disk particle (so each halo particle is
  // (halo mass / disk mass) / this heavier than a disk particle).
  double halo_particles_per_disk_particle = 3.0;
  // Barnes-Hut opening angle theta: a cell is treated as one point mass when
  // (cell width / distance) < theta. Smaller = more accurate, slower.
  double opening_angle = 0.5;
  // Softening = this factor x the mean particle spacing of the largest disk.
  double softening_spacing_factor = 0.3;
  // Pre-relaxation: before t = 0 each galaxy is evolved ALONE for this many Myr,
  // and the settled state becomes the starting condition. Freshly sampled
  // particles are not an exact equilibrium (a 2D Gaussian halo is not a true
  // distribution function), so the whole galaxy would otherwise collapse and
  // rebound for ~700 Myr. 0 = start from the raw sampled state.
  double relax_myr = 0.0;
  // Empty hole at the galaxy's centre, as a fraction of the galaxy radius: disk
  // AND halo particles are placed outside hole*radius (same mass, same count).
  // The fastest motion, close to the black hole, sets the timestep, so a bigger
  // hole allows a larger dt. (Halo-only is not enough: the halo particles that
  // happen to start near the black hole have ~20x the disk's peak acceleration.)
  double inner_hole_fraction = 0.0;
};

// What the application and the batch runner use unless told otherwise. The raw
// struct default is 0 so unit tests can inspect the sampled initial state.
constexpr double kDefaultRelaxMyr = 800.0;

struct ScenarioParams {
  GalaxyParams galaxy_a;
  std::optional<GalaxyParams> galaxy_b;  // empty = single-galaxy scenario
  EncounterParams encounter;             // ignored for a single galaxy

  MassModel mass_model = MassModel::Centralized;
  DistributedParams distributed;  // used only in the Distributed mode

  // Plummer softening length, applied to every interaction (Centralized mode;
  // the Distributed mode derives its own from the particle spacing).
  double softening_ly = 20.0;
  // Timestep = step_accuracy * (shortest dynamical time in the setup).
  // 0.05 is roughly 125 steps per orbit at the innermost star.
  double step_accuracy = 0.05;
  std::uint64_t seed = 12345;
};

inline double central_mass_msun(const GalaxyParams& g) {
  return g.total_mass_msun - g.num_stars * units::kStarMassMsun;
}

// Empty string if the scenario is valid, otherwise a human-readable reason.
// Simulation's constructor throws std::invalid_argument with the same message.
std::string validate_scenario(const ScenarioParams& p);

}  // namespace galaxy
