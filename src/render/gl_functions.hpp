#pragma once

// Minimal hand-rolled OpenGL function loader. Windows' opengl32 only exports
// GL 1.1, so everything newer is fetched through SDL_GL_GetProcAddress. The
// function-pointer typedefs come from SDL_opengl.h.
//
// Add a function by adding one line to GL_FUNCTION_LIST.

#include <SDL_opengl.h>

#define GL_FUNCTION_LIST(X)                                          \
  X(PFNGLGENVERTEXARRAYSPROC, GenVertexArrays, "glGenVertexArrays")  \
  X(PFNGLBINDVERTEXARRAYPROC, BindVertexArray, "glBindVertexArray")  \
  X(PFNGLDELETEVERTEXARRAYSPROC, DeleteVertexArrays, "glDeleteVertexArrays") \
  X(PFNGLGENBUFFERSPROC, GenBuffers, "glGenBuffers")                 \
  X(PFNGLBINDBUFFERPROC, BindBuffer, "glBindBuffer")                 \
  X(PFNGLBUFFERDATAPROC, BufferData, "glBufferData")                 \
  X(PFNGLBUFFERSUBDATAPROC, BufferSubData, "glBufferSubData")        \
  X(PFNGLDELETEBUFFERSPROC, DeleteBuffers, "glDeleteBuffers")        \
  X(PFNGLENABLEVERTEXATTRIBARRAYPROC, EnableVertexAttribArray, "glEnableVertexAttribArray") \
  X(PFNGLVERTEXATTRIBPOINTERPROC, VertexAttribPointer, "glVertexAttribPointer") \
  X(PFNGLCREATESHADERPROC, CreateShader, "glCreateShader")           \
  X(PFNGLSHADERSOURCEPROC, ShaderSource, "glShaderSource")           \
  X(PFNGLCOMPILESHADERPROC, CompileShader, "glCompileShader")        \
  X(PFNGLGETSHADERIVPROC, GetShaderiv, "glGetShaderiv")              \
  X(PFNGLGETSHADERINFOLOGPROC, GetShaderInfoLog, "glGetShaderInfoLog") \
  X(PFNGLDELETESHADERPROC, DeleteShader, "glDeleteShader")           \
  X(PFNGLCREATEPROGRAMPROC, CreateProgram, "glCreateProgram")        \
  X(PFNGLATTACHSHADERPROC, AttachShader, "glAttachShader")           \
  X(PFNGLLINKPROGRAMPROC, LinkProgram, "glLinkProgram")              \
  X(PFNGLGETPROGRAMIVPROC, GetProgramiv, "glGetProgramiv")           \
  X(PFNGLGETPROGRAMINFOLOGPROC, GetProgramInfoLog, "glGetProgramInfoLog") \
  X(PFNGLDELETEPROGRAMPROC, DeleteProgram, "glDeleteProgram")        \
  X(PFNGLUSEPROGRAMPROC, UseProgram, "glUseProgram")                 \
  X(PFNGLGETUNIFORMLOCATIONPROC, GetUniformLocation, "glGetUniformLocation") \
  X(PFNGLUNIFORM1FPROC, Uniform1f, "glUniform1f")                    \
  X(PFNGLUNIFORM2FPROC, Uniform2f, "glUniform2f")                    \
  X(PFNGLUNIFORM3FPROC, Uniform3f, "glUniform3f")

namespace gl {

#define X(type, name, str) extern type name;
GL_FUNCTION_LIST(X)
#undef X

// Needs a current GL context. Returns false (and fills `missing`) if any
// function could not be resolved.
bool load(const char** missing);

}  // namespace gl
