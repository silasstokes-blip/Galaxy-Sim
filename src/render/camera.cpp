#include "render/camera.hpp"

#include <algorithm>

namespace galaxy {

Camera fit_camera(Vec2 center, double extent_ly, int width_px, int height_px,
                  double margin) {
  // The disc must fit in the shorter dimension, so it spans half of that many
  // pixels from the center to the edge.
  const double half_px = 0.5 * static_cast<double>(std::min(width_px, height_px));
  Camera cam;
  cam.center = center;
  if (half_px > 0.0 && extent_ly > 0.0) cam.ly_per_pixel = extent_ly * margin / half_px;
  return cam;  // otherwise keep the default ly_per_pixel of 1.0 (degenerate input)
}

Vec2 ndc_scale(const Camera& cam, int width_px, int height_px) {
  // One pixel is 2/size of NDC on that axis, and one ly is 1/ly_per_pixel pixels.
  return {2.0 / (static_cast<double>(width_px) * cam.ly_per_pixel),
          2.0 / (static_cast<double>(height_px) * cam.ly_per_pixel)};
}

Vec2 world_to_ndc(const Camera& cam, Vec2 world, int width_px, int height_px) {
  const Vec2 s = ndc_scale(cam, width_px, height_px);
  const Vec2 d = world - cam.center;
  return {d.x * s.x, d.y * s.y};  // +y world is +y NDC, which is up on screen
}

}  // namespace galaxy
