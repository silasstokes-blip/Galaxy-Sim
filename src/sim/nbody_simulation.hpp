#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "sim/isimulation.hpp"
#include "sim/params.hpp"
#include "sim/quadtree.hpp"
#include "sim/vec2.hpp"

namespace galaxy {

// MassModel::Distributed backend. Each galaxy is
//   - a central black hole,
//   - a visible disk of equal-mass particles ("stars"), and
//   - a live dark-matter halo of heavier particles (invisible unless asked for).
// EVERY body attracts every other body; forces come from a Barnes-Hut tree
// rebuilt each step, and the integrator is the same kick-drift-kick leapfrog as
// the centralized mode. The halo carries most of the mass, which is what keeps
// the disk from fragmenting (see DistributedParams).
//
// Arrays are laid out [disk of galaxy 1, disk of galaxy 2, halo 1, halo 2,
// BH 1, BH 2], so the visible particles are one contiguous range at the front.
class NBodySimulation final : public ISimulation {
 public:
  // Throws std::invalid_argument on bad parameters.
  explicit NBodySimulation(const ScenarioParams& params, ProgressFn progress = {});

  // The tree keeps pointers into this object's arrays, so it cannot be copied.
  NBodySimulation(const NBodySimulation&) = delete;
  NBodySimulation& operator=(const NBodySimulation&) = delete;

  // One kick-drift-kick leapfrog step of length dt (years).
  void step(double dt);

  void queue_time(double years) override { carry_ += years; }
  double run_queued(int max_steps) override;
  void drop_backlog() override;
  double backlog_years() const override { return carry_; }

  double dt() const override { return dt_; }
  double time_years() const override { return time_; }

  // Energy with the potential energy taken from the tree (so it includes the
  // tree's approximation error); total_energy_direct() is the exact O(N^2) sum.
  double total_energy() const override;
  double initial_energy() const override { return e0_; }
  double total_energy_direct() const;
  // Disk particles only (the "stars"). Only the system measure exists here;
  // `galaxy` is NaN.
  BoundFractions bound_fractions() const override;
  double central_separation_ly() const override;  // between galaxy centres of mass

  Vec2 center_of_mass() const override;
  Vec2 total_momentum() const;
  double angular_momentum() const;  // z, about the origin
  double initial_extent() const override { return extent_; }
  double dark_matter_extent() const override { return dm_extent_; }

  RenderCounts render_counts() const override { return {n_disk_, n_bh_, n_halo_}; }
  void fill_render_points(Vec2 view_center, double lag_years, float* xy) const override;

  const ScenarioParams& params() const override { return params_; }

  // Inspection (tests, tools).
  std::size_t body_count() const { return x_.size(); }
  std::size_t particle_count() const { return n_disk_; }  // visible disk particles
  std::size_t halo_count() const { return n_halo_; }
  std::size_t black_hole_count() const { return n_bh_; }
  double softening_ly() const { return eps_; }
  const std::vector<double>& x() const { return x_; }
  const std::vector<double>& y() const { return y_; }
  const std::vector<double>& vx() const { return vx_; }
  const std::vector<double>& vy() const { return vy_; }
  const std::vector<double>& mass() const { return m_; }
  const std::vector<double>& ax() const { return ax_; }
  const std::vector<double>& ay() const { return ay_; }

 private:
  struct GalaxyModel;  // per-galaxy masses and analytic enclosed-mass profile

  void append_disk(const GalaxyModel& g, std::uint8_t gid, Vec2 centre, Vec2 bulk_vel,
                   std::mt19937_64& rng);
  void append_halo(const GalaxyModel& g, std::uint8_t gid, Vec2 centre, Vec2 bulk_vel,
                   std::mt19937_64& rng);
  // Rebuild the tree at the current positions and refresh every acceleration.
  // Returns the largest |a|^2.
  double compute_forces();
  const std::vector<double>& potentials() const;
  double compute_initial_dt(double max_acc2) const;

  ScenarioParams params_;
  double eps_ = 0.0, eps2_ = 0.0;
  std::size_t n_disk_ = 0, n_halo_ = 0, n_bh_ = 0;
  std::vector<double> x_, y_, vx_, vy_, ax_, ay_, m_;
  std::vector<std::uint8_t> gal_;  // which galaxy each body started in (0 or 1)
  BarnesHutTree tree_;

  double dt_ = 0.0;
  double aref2_ = 0.0;  // |a|^2 that dt_ was sized for; see step()
  double time_ = 0.0;
  double carry_ = 0.0;
  double e0_ = 0.0;
  double extent_ = 0.0;     // visible bodies (disk + black holes)
  double dm_extent_ = 0.0;  // including the halo

  // Potentials are O(N log N) to get, shared by total_energy() and
  // bound_fractions() and valid for one step count.
  mutable std::vector<double> phi_;
  mutable std::uint64_t phi_step_ = ~std::uint64_t{0};
  std::uint64_t step_count_ = 0;
};

}  // namespace galaxy
