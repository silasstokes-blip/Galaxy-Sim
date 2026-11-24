#pragma once

#include <cstddef>
#include <string>

#include "render/gl_functions.hpp"

namespace galaxy {

// Draws 1-pixel points with additive blending. Vertex positions are 2D floats
// RELATIVE TO THE CAMERA CENTER (the CPU subtracts the center in double
// precision before converting to float); the shader multiplies by the camera
// scale to get NDC.
class PointRenderer {
 public:
  PointRenderer() = default;
  PointRenderer(const PointRenderer&) = delete;
  PointRenderer& operator=(const PointRenderer&) = delete;
  ~PointRenderer() { shutdown(); }

  // Needs gl::load() to have succeeded. On failure returns false and fills `error`.
  bool init(std::string& error);
  void shutdown();

  // Replace the vertex data: `count` points, 2 floats each.
  void upload(const float* xy, std::size_t count);

  // Draw points [first, first + count) with the given NDC scale. Each point
  // adds (r, g, b) to the framebuffer, so overlapping points brighten.
  void draw(std::size_t first, std::size_t count, float scale_x, float scale_y, float r, float g,
            float b) const;

 private:
  GLuint vao_ = 0, vbo_ = 0, program_ = 0;
  GLint u_scale_ = -1, u_color_ = -1;
  std::size_t capacity_ = 0;  // points the VBO can currently hold
};

}  // namespace galaxy
