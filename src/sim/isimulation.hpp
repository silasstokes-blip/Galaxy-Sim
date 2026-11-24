#pragma once

#include <cstddef>
#include <functional>
#include <memory>

#include "sim/params.hpp"
#include "sim/vec2.hpp"

namespace galaxy {

// Fractions of bodies (0..1) that are gravitationally bound, by two definitions.
struct BoundFractions {
  // Total energy < 0 in the centre-of-mass frame, against the combined
  // potential of everything else: "has not reached escape velocity".
  double system = 1.0;
  // Energy < 0 relative to at least ONE central body, in that body's rest
  // frame: "still bound to a galaxy". NaN when the mode cannot measure it
  // (the distributed-mass mode reports the system measure only).
  double galaxy = 1.0;
};

// Bodies are drawn in three groups, in this order: `regular` (stars / disk
// particles), `highlighted` (black holes / central bodies), then `halo` (dark
// matter, drawn only on request).
struct RenderCounts {
  std::size_t regular = 0;
  std::size_t highlighted = 0;
  std::size_t halo = 0;
};

// Everything the app, batch runner and CSV logging need from a simulation,
// independent of how gravity is computed.
class ISimulation {
 public:
  virtual ~ISimulation() = default;
  ISimulation() = default;
  ISimulation(const ISimulation&) = delete;
  ISimulation& operator=(const ISimulation&) = delete;

  // Time stepping in chunks: queue simulated time, run up to max_steps of it
  // (repeatable within a frame), then drop any backlog beyond one step.
  virtual void queue_time(double years) = 0;
  virtual double run_queued(int max_steps) = 0;  // returns years advanced
  virtual void drop_backlog() = 0;
  // Queued time not yet simulated; after drop_backlog() it is < dt(). A
  // renderer extrapolates positions by pos + vel * backlog_years().
  virtual double backlog_years() const = 0;

  virtual double dt() const = 0;
  virtual double time_years() const = 0;

  // Diagnostics (sim unit system; see units.hpp).
  virtual double total_energy() const = 0;
  virtual double initial_energy() const = 0;
  virtual BoundFractions bound_fractions() const = 0;
  // Distance between the two galaxies' centres (centres of mass in the
  // distributed mode); 0 for a single galaxy.
  virtual double central_separation_ly() const = 0;

  // Camera helpers.
  virtual Vec2 center_of_mass() const = 0;     // of the whole system
  // Max |pos - COM| at t = 0 over the VISIBLE bodies (not dark matter).
  virtual double initial_extent() const = 0;
  // Same, including dark matter; larger when a halo extends past the disks.
  virtual double dark_matter_extent() const { return initial_extent(); }

  virtual RenderCounts render_counts() const = 0;
  // Writes 2 floats per body in render_counts() order (regular, highlighted,
  // halo): position extrapolated by lag_years, relative to view_center.
  virtual void fill_render_points(Vec2 view_center, double lag_years, float* xy) const = 0;

  virtual const ScenarioParams& params() const = 0;

  double energy_ratio() const { return total_energy() / initial_energy(); }

  // Convenience: queue, run (capped) and drop the backlog in one call. Returns
  // the years actually advanced; a shortfall means the sim is running slower
  // than requested and is not made up later.
  double advance(double target_years, int max_steps) {
    queue_time(target_years);
    const double advanced = run_queued(max_steps);
    drop_backlog();
    return advanced;
  }
};

// Builds the backend selected by params.mass_model. Throws
// std::invalid_argument for invalid parameters.
// Called with a fraction in [0, 1] while a long constructor (the distributed
// mode's pre-relaxation) is running. May throw to abort construction.
using ProgressFn = std::function<void(double)>;

std::unique_ptr<ISimulation> make_simulation(const ScenarioParams& params, ProgressFn progress = {});

}  // namespace galaxy
