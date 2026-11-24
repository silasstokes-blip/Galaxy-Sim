#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace galaxy {

// Barnes-Hut quadtree over 2D point masses.
//
// Gravity is the 3D Newtonian 1/r^2 force with Plummer softening evaluated for
// bodies on a plane (same convention as the rest of the sim). Distant cells
// are replaced by a single point mass at their centre of mass (monopole only),
// so the force error scales roughly as theta^2.
//
// The tree stores indices, not positions: the arrays passed to build() must
// stay alive and unmoved until the next build(). Rebuild every step.
class BarnesHutTree {
 public:
  void build(const double* x, const double* y, const double* m, std::size_t n);

  // Acceleration at (x, y). `self` is the index of the body at that point (its
  // own mass is excluded), or -1 for an arbitrary point.
  void accel(std::int32_t self, double x, double y, double theta, double eps2, double& ax,
             double& ay) const;

  // Gravitational potential at (x, y), excluding `self`. (Negative.)
  double potential(std::int32_t self, double x, double y, double theta, double eps2) const;

  std::size_t node_count() const { return nodes_.size(); }

  // Body indices in depth-first (spatial) order. Looping over targets in this
  // order makes neighbouring bodies walk the same cells, which keeps the tree in
  // cache. Falls back to 0..n-1 if some bodies were merged into shared leaves.
  const std::vector<std::int32_t>& leaf_order() const { return order_; }

 private:
  struct Node {
    double px = 0, py = 0, half = 0;  // cell centre and half-width
    double mass = 0;
    double cx = 0, cy = 0;            // sum of m*pos while building, then the COM
    std::int32_t child = -1;          // first of 4 contiguous children, or -1 (leaf)
    std::int32_t body = -1;           // leaf: single body index; -2 = merged bodies; -1 = empty
  };

  void insert(std::int32_t i);
  template <class F>
  void walk(std::int32_t self, double x, double y, double theta2, F&& visit) const;

  std::vector<Node> nodes_;
  std::vector<std::int32_t> order_;
  const double* x_ = nullptr;
  const double* y_ = nullptr;
  const double* m_ = nullptr;
  std::size_t n_ = 0;
};

}  // namespace galaxy
