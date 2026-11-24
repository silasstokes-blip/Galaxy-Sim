#pragma once

#include <cstddef>
#include <random>
#include <vector>

#include "sim/isimulation.hpp"
#include "sim/params.hpp"
#include "sim/vec2.hpp"

namespace galaxy {

// A massive body that attracts everything. There are one or two of these.
struct CentralBody {
  Vec2 pos;
  Vec2 vel;
  Vec2 acc;
  double mass = 0.0;  // Msun
};

// Stars are massless-for-dynamics test particles: they feel the central
// bodies but exert no force on anything. Structure-of-arrays for speed.
struct StarArrays {
  std::vector<double> x, y, vx, vy;
  std::vector<double> ax, ay;  // cached acceleration at the current position
  std::size_t size() const { return x.size(); }
};

// MassModel::Centralized backend: one or two central bodies on a Kepler orbit,
// stars as test particles. O(N) per step.
class Simulation final : public ISimulation {
 public:
  // Throws std::invalid_argument on bad parameters.
  explicit Simulation(const ScenarioParams& params);

  // One kick-drift-kick leapfrog step of length dt (years).
  void step(double dt);

  void queue_time(double years) override { carry_ += years; }
  double run_queued(int max_steps) override;
  void drop_backlog() override;
  double backlog_years() const override { return carry_; }

  double dt() const override { return dt_; }
  double time_years() const override { return time_; }

  double total_energy() const override;
  double initial_energy() const override { return e0_; }
  BoundFractions bound_fractions() const override;  // O(N)
  double central_separation_ly() const override;    // 0 for a single galaxy
  double angular_momentum() const;                   // z, about origin

  Vec2 center_of_mass() const override;  // of the central bodies
  Vec2 center_of_mass_velocity() const;
  double initial_extent() const override { return extent_; }

  RenderCounts render_counts() const override { return {stars_.size(), centrals_.size(), 0}; }
  void fill_render_points(Vec2 view_center, double lag_years, float* xy) const override;

  const std::vector<CentralBody>& centrals() const { return centrals_; }
  const StarArrays& stars() const { return stars_; }
  const ScenarioParams& params() const override { return params_; }

 private:
  struct SourcePoint {
    double x, y, gm;
  };

  void build_centrals();
  void populate_galaxy(const GalaxyParams& g, const CentralBody& c,
                       std::mt19937_64& rng);
  void compute_central_acc();
  void step_centrals(double dt);
  double compute_dt() const;
  double compute_extent() const;

  ScenarioParams params_;
  double eps2_ = 0.0;
  std::vector<CentralBody> centrals_;
  StarArrays stars_;
  double dt_ = 0.0;
  double time_ = 0.0;
  double carry_ = 0.0;  // sub-step remainder of requested time
  double e0_ = 0.0;
  double extent_ = 0.0;
};

}  // namespace galaxy
