#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include "sim/isimulation.hpp"
#include "sim/nbody_simulation.hpp"
#include "sim/simulation.hpp"

namespace galaxy {

namespace {

std::string galaxy_error(const GalaxyParams& g, const std::string& name, bool distributed) {
  if (g.num_stars < 0) return name + ": number of stars must be >= 0";
  if (!(g.radius_ly > 0.0)) return name + ": radius must be > 0";
  if (distributed) {
    if (g.num_stars < 1) return name + ": need at least 1 particle";
    if (!(g.total_mass_msun > 0.0)) return name + ": total mass must be > 0";
  } else {
    if (!(g.inner_fraction >= 0.0 && g.inner_fraction < 1.0))
      return name + ": inner_fraction must be in [0, 1)";
    if (!(central_mass_msun(g) > 0.0))
      return name + ": total mass must exceed the combined star mass";
  }
  return {};
}

}  // namespace

std::string validate_scenario(const ScenarioParams& p) {
  const bool distributed = p.mass_model == MassModel::Distributed;

  if (auto e = galaxy_error(p.galaxy_a, "galaxy 1", distributed); !e.empty()) return e;
  if (p.galaxy_b) {
    if (auto e = galaxy_error(*p.galaxy_b, "galaxy 2", distributed); !e.empty()) return e;
    const auto& enc = p.encounter;
    if (!(enc.separation_ly > 0.0)) return "separation must be > 0";
    if (!(std::abs(enc.impact_parameter_ly) < enc.separation_ly))
      return "|impact parameter| must be less than the separation";
    if (!(enc.relative_speed_km_s >= 0.0)) return "relative speed must be >= 0";
  }

  if (distributed) {
    const auto& d = p.distributed;
    if (!(d.bh_mass_fraction >= 0.0 && d.bh_mass_fraction <= 0.5))
      return "black hole mass fraction must be in [0, 0.5]";
    if (!(d.dispersion_fraction >= 0.0 && d.dispersion_fraction <= 1.0))
      return "velocity dispersion must be in [0, 1]";
    if (!(d.halo_mass_fraction >= 0.0 && d.halo_mass_fraction <= 0.95))
      return "dark matter halo fraction must be in [0, 0.95]";
    if (d.bh_mass_fraction + d.halo_mass_fraction > 0.98)
      return "black hole + halo must leave at least 2% of the mass for the disk";
    if (d.halo_mass_fraction > 0.0) {
      if (!(d.halo_radius_fraction >= 0.5 && d.halo_radius_fraction <= 10.0))
        return "halo radius must be 0.5 to 10 x the galaxy radius";
      if (!(d.halo_core_fraction > 0.0 && d.halo_core_fraction <= 5.0))
        return "halo core must be in (0, 5] x the galaxy radius";
      if (!(d.halo_particles_per_disk_particle >= 0.1 && d.halo_particles_per_disk_particle <= 20.0))
        return "halo particles per disk particle must be 0.1 to 20";
    }
    if (!(d.inner_hole_fraction >= 0.0 && d.inner_hole_fraction < 0.9))
      return "inner hole must be in [0, 0.9) x the galaxy radius";
    if (!(d.relax_myr >= 0.0 && d.relax_myr <= 5000.0))
      return "pre-relaxation must be 0 to 5000 Myr";
    if (!(d.opening_angle > 0.0 && d.opening_angle <= 2.0))
      return "opening angle must be in (0, 2]";
    if (!(d.softening_spacing_factor > 0.0)) return "softening factor must be > 0";
    if (d.profile == DiskProfile::Exponential &&
        !(d.scale_length_fraction > 0.0 && d.scale_length_fraction <= 2.0))
      return "scale length must be in (0, 2] x the radius";
  } else {
    if (!(p.softening_ly > 0.0)) return "softening must be > 0";
  }
  if (!(p.step_accuracy > 0.0)) return "step_accuracy must be > 0";
  return {};
}

std::unique_ptr<ISimulation> make_simulation(const ScenarioParams& p, ProgressFn progress) {
  if (p.mass_model == MassModel::Distributed)
    return std::make_unique<NBodySimulation>(p, std::move(progress));
  return std::make_unique<Simulation>(p);
}

}  // namespace galaxy
