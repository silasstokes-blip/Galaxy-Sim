// Spec for src/render/camera.cpp. These fail until the camera is implemented.

#include <cmath>
#include <cstdio>

#include "render/camera.hpp"

using namespace galaxy;

static int g_failures = 0;

#define CHECK_NEAR(a, b, tol)                                                     \
  do {                                                                            \
    const double a_ = (a), b_ = (b);                                              \
    if (!(std::abs(a_ - b_) <= (tol))) {                                          \
      std::printf("    FAIL %s:%d: %s = %.9g, expected %.9g\n", __FILE__, __LINE__, \
                  #a, a_, b_);                                                    \
      ++g_failures;                                                               \
    }                                                                             \
  } while (0)

static void test_fit_keeps_center_and_scale() {
  const Vec2 c{3.0, -4.0};
  // Square viewport: disc radius 1e5 * margin 1.0 spans half of 1000 px.
  Camera sq = fit_camera(c, 1.0e5, 1000, 1000, 1.0);
  CHECK_NEAR(sq.center.x, 3.0, 1e-12);
  CHECK_NEAR(sq.center.y, -4.0, 1e-12);
  CHECK_NEAR(sq.ly_per_pixel, 1.0e5 / 500.0, 1e-9);

  // Margin enlarges the visible disc.
  Camera m = fit_camera(c, 1.0e5, 1000, 1000, 1.05);
  CHECK_NEAR(m.ly_per_pixel, 1.0e5 * 1.05 / 500.0, 1e-9);

  // The shorter dimension limits, for wide and tall windows alike.
  Camera wide = fit_camera(c, 1.0e5, 2000, 1000, 1.05);
  Camera tall = fit_camera(c, 1.0e5, 1000, 2000, 1.05);
  CHECK_NEAR(wide.ly_per_pixel, m.ly_per_pixel, 1e-9);
  CHECK_NEAR(tall.ly_per_pixel, m.ly_per_pixel, 1e-9);
}

static void test_world_to_ndc() {
  const Vec2 c{10.0, 20.0};
  Camera cam = fit_camera(c, 1.0e5, 1000, 1000, 1.0);

  Vec2 mid = world_to_ndc(cam, c, 1000, 1000);
  CHECK_NEAR(mid.x, 0.0, 1e-12);
  CHECK_NEAR(mid.y, 0.0, 1e-12);

  // Edge of the disc lands on the edge of the viewport; +y is up.
  Vec2 right = world_to_ndc(cam, {c.x + 1.0e5, c.y}, 1000, 1000);
  Vec2 up = world_to_ndc(cam, {c.x, c.y + 1.0e5}, 1000, 1000);
  CHECK_NEAR(right.x, 1.0, 1e-9);
  CHECK_NEAR(right.y, 0.0, 1e-9);
  CHECK_NEAR(up.x, 0.0, 1e-9);
  CHECK_NEAR(up.y, 1.0, 1e-9);

  // In a 2000 x 1000 window the same x offset covers only half the width.
  Camera wide = fit_camera(c, 1.0e5, 2000, 1000, 1.0);
  Vec2 r2 = world_to_ndc(wide, {c.x + 1.0e5, c.y}, 2000, 1000);
  Vec2 u2 = world_to_ndc(wide, {c.x, c.y + 1.0e5}, 2000, 1000);
  CHECK_NEAR(r2.x, 0.5, 1e-9);
  CHECK_NEAR(u2.y, 1.0, 1e-9);
}

static void test_ndc_scale_matches() {
  const Vec2 c{-7.0, 5.0};
  Camera cam = fit_camera(c, 5.0e4, 1600, 900, 1.05);
  const Vec2 s = ndc_scale(cam, 1600, 900);
  const Vec2 offsets[] = {{1234.5, -678.9}, {-4.0e4, 3.0e4}, {0.0, 0.0}};
  for (const Vec2& d : offsets) {
    Vec2 n = world_to_ndc(cam, c + d, 1600, 900);
    CHECK_NEAR(n.x, d.x * s.x, 1e-9);
    CHECK_NEAR(n.y, d.y * s.y, 1e-9);
  }
}

#define RUN(test)                                                             \
  do {                                                                        \
    const int before = g_failures;                                            \
    std::printf("[ RUN  ] %s\n", #test);                                      \
    test();                                                                   \
    std::printf("[ %s ] %s\n", g_failures == before ? " OK " : "FAIL", #test); \
  } while (0)

int main() {
  RUN(test_fit_keeps_center_and_scale);
  RUN(test_world_to_ndc);
  RUN(test_ndc_scale_matches);
  std::printf("\n%s (%d failed check%s)\n", g_failures ? "FAILED" : "ALL PASSED",
              g_failures, g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
