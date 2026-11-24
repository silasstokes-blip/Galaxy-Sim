#include "render/gl_functions.hpp"

#include <SDL.h>

namespace gl {

#define X(type, name, str) type name = nullptr;
GL_FUNCTION_LIST(X)
#undef X

bool load(const char** missing) {
#define X(type, name, str)                                              \
  name = reinterpret_cast<type>(SDL_GL_GetProcAddress(str));            \
  if (!name) {                                                          \
    if (missing) *missing = str;                                        \
    return false;                                                       \
  }
  GL_FUNCTION_LIST(X)
#undef X
  return true;
}

}  // namespace gl
