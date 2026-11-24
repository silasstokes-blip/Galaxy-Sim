#include "sim/quadtree.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "sim/units.hpp"

namespace galaxy {

namespace {
// Bodies that share a position (to within 2^-48 of the root width) are merged
// into one leaf instead of subdividing forever.
constexpr int kMaxDepth = 48;

inline int quadrant(double px, double py, double x, double y) {
  return (x > px ? 1 : 0) | (y > py ? 2 : 0);
}
}  // namespace

void BarnesHutTree::build(const double* x, const double* y, const double* m, std::size_t n) {
  x_ = x;
  y_ = y;
  m_ = m;
  n_ = n;
  nodes_.clear();
  if (n == 0) return;

  double minx = x[0], maxx = x[0], miny = y[0], maxy = y[0];
  for (std::size_t i = 1; i < n; ++i) {
    minx = std::min(minx, x[i]);
    maxx = std::max(maxx, x[i]);
    miny = std::min(miny, y[i]);
    maxy = std::max(maxy, y[i]);
  }
  double half = 0.5 * std::max(maxx - minx, maxy - miny);
  half = half > 0.0 ? half * (1.0 + 1e-9) : 1.0;  // slightly padded so edge bodies are inside

  nodes_.reserve(3 * n + 16);
  Node root;
  root.px = 0.5 * (minx + maxx);
  root.py = 0.5 * (miny + maxy);
  root.half = half;
  nodes_.push_back(root);

  for (std::size_t i = 0; i < n; ++i) insert(static_cast<std::int32_t>(i));

  // cx, cy held sums of m*pos; turn them into centres of mass.
  for (Node& nd : nodes_) {
    if (nd.mass > 0.0) {
      nd.cx /= nd.mass;
      nd.cy /= nd.mass;
    }
  }

  // Depth-first leaf order (see leaf_order()).
  order_.clear();
  order_.reserve(n);
  std::int32_t stack[256];
  int sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const Node& nd = nodes_[stack[--sp]];
    if (nd.child >= 0) {
      for (int q = 3; q >= 0; --q) stack[sp++] = nd.child + q;  // visit quadrant 0 first
    } else if (nd.body >= 0) {
      order_.push_back(nd.body);
    }
  }
  if (order_.size() != n) {  // merged leaves hide some indices: use the plain order
    order_.resize(n);
    for (std::size_t i = 0; i < n; ++i) order_[i] = static_cast<std::int32_t>(i);
  }
}

void BarnesHutTree::insert(std::int32_t i) {
  const double xi = x_[i], yi = y_[i], mi = m_[i];
  std::int32_t node = 0;
  bool add = true;  // false right after subdividing a node that already counted body i
  int depth = 0;

  for (;;) {
    if (add) {
      Node& n = nodes_[node];
      n.mass += mi;
      n.cx += mi * xi;
      n.cy += mi * yi;
    }
    Node& n = nodes_[node];

    if (n.child >= 0) {  // internal: descend
      node = n.child + quadrant(n.px, n.py, xi, yi);
      add = true;
      ++depth;
      continue;
    }
    if (n.body == -1) {  // empty leaf: take the body
      n.body = i;
      return;
    }
    if (n.body == -2 || depth >= kMaxDepth) {  // merged leaf
      n.body = -2;
      return;
    }

    // Occupied leaf: subdivide, move the resident body into its child, then
    // loop again to place body i.
    const std::int32_t j = n.body;
    const double h = 0.5 * n.half, px = n.px, py = n.py;
    const auto first = static_cast<std::int32_t>(nodes_.size());
    for (int q = 0; q < 4; ++q) {
      Node c;
      c.px = px + ((q & 1) ? h : -h);
      c.py = py + ((q & 2) ? h : -h);
      c.half = h;
      nodes_.push_back(c);  // may reallocate: no references held across this
    }
    nodes_[node].child = first;
    nodes_[node].body = -1;

    Node& cj = nodes_[first + quadrant(px, py, x_[j], y_[j])];
    cj.mass = m_[j];
    cj.cx = m_[j] * x_[j];
    cj.cy = m_[j] * y_[j];
    cj.body = j;
    add = false;
  }
}

template <class F>
void BarnesHutTree::walk(std::int32_t self, double x, double y, double theta2, F&& visit) const {
  if (nodes_.empty()) return;
  std::int32_t stack[256];  // depth <= 48, 3 pending siblings per level
  int sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const Node& n = nodes_[stack[--sp]];
    if (n.mass <= 0.0) continue;
    const double dx = n.cx - x, dy = n.cy - y;
    const double r2 = dx * dx + dy * dy;

    bool as_point;
    if (n.child < 0) {
      if (n.body == self) continue;  // never act on yourself
      as_point = true;
    } else {
      const double s = 2.0 * n.half;
      // Accept the cell only if it is far enough AND the target is not inside it.
      as_point = s * s < theta2 * r2 &&
                 !(std::fabs(x - n.px) <= n.half && std::fabs(y - n.py) <= n.half);
    }
    if (as_point) {
      visit(n.mass, dx, dy, r2);
    } else {
      for (int q = 0; q < 4; ++q) stack[sp++] = n.child + q;
    }
  }
}

void BarnesHutTree::accel(std::int32_t self, double x, double y, double theta, double eps2,
                          double& ax, double& ay) const {
  double axs = 0.0, ays = 0.0;
  walk(self, x, y, theta * theta, [&](double m, double dx, double dy, double r2) {
    const double d2 = r2 + eps2;
    const double inv = m / (d2 * std::sqrt(d2));
    axs += dx * inv;
    ays += dy * inv;
  });
  ax = units::kG * axs;
  ay = units::kG * ays;
}

double BarnesHutTree::potential(std::int32_t self, double x, double y, double theta,
                                double eps2) const {
  double p = 0.0;
  walk(self, x, y, theta * theta,
       [&](double m, double, double, double r2) { p += m / std::sqrt(r2 + eps2); });
  return -units::kG * p;
}

}  // namespace galaxy
