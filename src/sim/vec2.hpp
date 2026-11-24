#pragma once

#include <cmath>

namespace galaxy {

struct Vec2 {
  double x = 0.0;
  double y = 0.0;

  constexpr Vec2 operator+(Vec2 o) const { return {x + o.x, y + o.y}; }
  constexpr Vec2 operator-(Vec2 o) const { return {x - o.x, y - o.y}; }
  constexpr Vec2 operator-() const { return {-x, -y}; }
  constexpr Vec2 operator*(double s) const { return {x * s, y * s}; }
  constexpr Vec2& operator+=(Vec2 o) { x += o.x; y += o.y; return *this; }

  constexpr double length2() const { return x * x + y * y; }
  double length() const { return std::sqrt(length2()); }
};

constexpr Vec2 operator*(double s, Vec2 v) { return v * s; }
constexpr double cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }

}  // namespace galaxy
