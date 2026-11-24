cmake_minimum_required(VERSION 3.20)
project(GalaxySim LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(GALAXY_USE_OPENMP "Parallelise the star update loop with OpenMP" ON)
option(GALAXY_BUILD_APP "Build the SDL2 + ImGui application (downloads dependencies)" ON)

function(galaxy_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

# ---- simulation library (no graphics dependencies) -------------------------
# NOTE: deliberately no -ffast-math. We want reproducible floating point so
# energy-conservation tests mean something.
add_library(galaxy_sim STATIC
  src/sim/simulation.cpp
  src/sim/scenario.cpp
  src/sim/quadtree.cpp
  src/sim/nbody_simulation.cpp)
target_include_directories(galaxy_sim PUBLIC src)
galaxy_warnings(galaxy_sim)

if(GALAXY_USE_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(galaxy_sim PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

# ---- camera math (pure, no GL/SDL) -----------------------------------------
add_library(galaxy_camera STATIC src/render/camera.cpp)
target_link_libraries(galaxy_camera PUBLIC galaxy_sim)
galaxy_warnings(galaxy_camera)

# ---- CSV logging (plain file I/O, no GL/SDL) --------------------------------
add_library(galaxy_io STATIC src/app/csv_logger.cpp)
target_include_directories(galaxy_io PUBLIC src)
galaxy_warnings(galaxy_io)

# ---- headless batch runner (no GL/SDL) --------------------------------------
add_executable(galaxy_batch src/app/batch_main.cpp)
target_link_libraries(galaxy_batch PRIVATE galaxy_sim galaxy_io)
galaxy_warnings(galaxy_batch)

# ---- tests -----------------------------------------------------------------
enable_testing()
add_executable(sim_tests tests/test_main.cpp)
target_link_libraries(sim_tests PRIVATE galaxy_sim)
add_test(NAME sim_tests COMMAND sim_tests)

add_executable(nbody_tests tests/test_nbody.cpp)
target_link_libraries(nbody_tests PRIVATE galaxy_sim)
add_test(NAME nbody_tests COMMAND nbody_tests)

add_executable(camera_tests tests/test_camera.cpp)
target_link_libraries(camera_tests PRIVATE galaxy_camera)
add_test(NAME camera_tests COMMAND camera_tests)

add_executable(logger_tests tests/test_logger.cpp)
target_link_libraries(logger_tests PRIVATE galaxy_io)
add_test(NAME logger_tests COMMAND logger_tests)

# ---- application -----------------------------------------------------------
if(GALAXY_BUILD_APP)
  include(FetchContent)

  set(SDL_SHARED OFF CACHE BOOL "" FORCE)
  set(SDL_STATIC ON CACHE BOOL "" FORCE)
  set(SDL_TEST OFF CACHE BOOL "" FORCE)
  set(SDL_TESTS OFF CACHE BOOL "" FORCE)
  set(SDL2_DISABLE_INSTALL ON CACHE BOOL "" FORCE)

  FetchContent_Declare(SDL2
    GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
    GIT_TAG release-2.30.9
    GIT_SHALLOW TRUE)
  FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.5
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(SDL2 imgui)

  add_library(imgui_lib STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_sdl2.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)
  target_include_directories(imgui_lib PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
  target_link_libraries(imgui_lib PUBLIC SDL2::SDL2-static)

  add_library(galaxy_render STATIC
    src/render/gl_functions.cpp
    src/render/point_renderer.cpp)
  target_include_directories(galaxy_render PUBLIC src)
  target_link_libraries(galaxy_render PUBLIC SDL2::SDL2-static)
  if(WIN32)
    target_link_libraries(galaxy_render PUBLIC opengl32)
  else()
    find_package(OpenGL REQUIRED)
    target_link_libraries(galaxy_render PUBLIC OpenGL::GL)
  endif()
  galaxy_warnings(galaxy_render)

  add_executable(galaxy_app src/app/main.cpp src/app/setup_panel.cpp)
  target_link_libraries(galaxy_app PRIVATE
    SDL2::SDL2main galaxy_render galaxy_camera galaxy_io galaxy_sim imgui_lib)
  galaxy_warnings(galaxy_app)
endif()
