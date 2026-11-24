#include "render/point_renderer.hpp"

#include <vector>

namespace galaxy {

namespace {

const char* kVertexSrc = R"(#version 330 core
layout(location = 0) in vec2 a_pos;   // world offset from camera center, ly
uniform vec2 u_scale;                 // NDC per ly, per axis
void main() {
  gl_Position = vec4(a_pos * u_scale, 0.0, 1.0);
}
)";

const char* kFragmentSrc = R"(#version 330 core
uniform vec3 u_color;
out vec4 frag_color;
void main() {
  frag_color = vec4(u_color, 1.0);
}
)";

GLuint compile(GLenum type, const char* src, std::string& error) {
  const GLuint s = gl::CreateShader(type);
  gl::ShaderSource(s, 1, &src, nullptr);
  gl::CompileShader(s);
  GLint ok = 0;
  gl::GetShaderiv(s, GL_COMPILE_STATUS, &ok);
  if (!ok) {
    char log[1024];
    gl::GetShaderInfoLog(s, sizeof log, nullptr, log);
    error = std::string(type == GL_VERTEX_SHADER ? "vertex" : "fragment") +
            " shader: " + log;
    gl::DeleteShader(s);
    return 0;
  }
  return s;
}

}  // namespace

bool PointRenderer::init(std::string& error) {
  const GLuint vs = compile(GL_VERTEX_SHADER, kVertexSrc, error);
  if (!vs) return false;
  const GLuint fs = compile(GL_FRAGMENT_SHADER, kFragmentSrc, error);
  if (!fs) {
    gl::DeleteShader(vs);
    return false;
  }

  program_ = gl::CreateProgram();
  gl::AttachShader(program_, vs);
  gl::AttachShader(program_, fs);
  gl::LinkProgram(program_);
  gl::DeleteShader(vs);
  gl::DeleteShader(fs);

  GLint ok = 0;
  gl::GetProgramiv(program_, GL_LINK_STATUS, &ok);
  if (!ok) {
    char log[1024];
    gl::GetProgramInfoLog(program_, sizeof log, nullptr, log);
    error = std::string("program link: ") + log;
    gl::DeleteProgram(program_);
    program_ = 0;
    return false;
  }
  u_scale_ = gl::GetUniformLocation(program_, "u_scale");
  u_color_ = gl::GetUniformLocation(program_, "u_color");

  gl::GenVertexArrays(1, &vao_);
  gl::GenBuffers(1, &vbo_);
  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);
  gl::BindVertexArray(0);
  return true;
}

void PointRenderer::shutdown() {
  if (vbo_) gl::DeleteBuffers(1, &vbo_);
  if (vao_) gl::DeleteVertexArrays(1, &vao_);
  if (program_) gl::DeleteProgram(program_);
  vbo_ = vao_ = program_ = 0;
  capacity_ = 0;
}

void PointRenderer::upload(const float* xy, std::size_t count) {
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  const auto bytes = static_cast<GLsizeiptr>(count * 2 * sizeof(float));
  if (count > capacity_) {
    gl::BufferData(GL_ARRAY_BUFFER, bytes, xy, GL_STREAM_DRAW);
    capacity_ = count;
  } else {
    gl::BufferSubData(GL_ARRAY_BUFFER, 0, bytes, xy);
  }
}

void PointRenderer::draw(std::size_t first, std::size_t count, float scale_x, float scale_y,
                         float r, float g, float b) const {
  if (count == 0) return;
  gl::UseProgram(program_);
  gl::Uniform2f(u_scale_, scale_x, scale_y);
  gl::Uniform3f(u_color_, r, g, b);
  gl::BindVertexArray(vao_);
  glDrawArrays(GL_POINTS, static_cast<GLint>(first), static_cast<GLsizei>(count));
  gl::BindVertexArray(0);
}

}  // namespace galaxy
