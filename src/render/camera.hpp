#pragma once

#include "sim/vec2.hpp"

namespace galaxy {

// A fixed top-down 2D camera. World units are light-years, +y is up on screen.
// NDC (normalized device coordinates) run from -1 to +1 on both axes, with
// (0, 0) at the middle of the viewport.
struct Camera {
  Vec2 center;                // world point shown at the middle of the viewport
  double ly_per_pixel = 1.0;  // world distance covered by one screen pixel
};

// YOUR PIECE (src/render/camera.cpp). Specs are pinned down by
// tests/test_camera.cpp.

// Build a camera centered on `center` that shows a disc of radius
// extent_ly * margin entirely inside a width_px x height_px viewport. The
// limiting dimension is the SHORTER one, so the disc touches (margin == 1) the
// top/bottom edges of a wide window, or the left/right edges of a tall one.
Camera fit_camera(Vec2 center, double extent_ly, int width_px, int height_px,
                  double margin = 1.05);

// Where a world point lands in NDC for this camera and viewport size.
Vec2 world_to_ndc(const Camera& cam, Vec2 world, int width_px, int height_px);

// Per-axis factor such that  ndc = (world - cam.center) * ndc_scale(...).
// The renderer hands this to the vertex shader, so the GPU does the transform.
Vec2 ndc_scale(const Camera& cam, int width_px, int height_px);

}  // namespace galaxy
