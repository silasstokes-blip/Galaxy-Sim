#include <SDL.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "app/csv_logger.hpp"
#include "app/setup_panel.hpp"
#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_sdl2.h"
#include "render/camera.hpp"
#include "render/gl_functions.hpp"
#include "render/point_renderer.hpp"
#include "sim/isimulation.hpp"

using namespace galaxy;

namespace {

constexpr float kMinSpeedMyr = 1.0f;     // slider range, Myr of sim time per second
constexpr float kMaxSpeedMyr = 100.0f;
constexpr double kStepBudgetSeconds = 0.012;  // max wall time spent stepping per frame
constexpr float kStarIntensity = 0.45f;
constexpr float kCentralIntensity = 1.0f;
// Dark matter: blue, dim per particle so the (many) halo particles read as a haze.
constexpr float kHaloR = 0.06f, kHaloG = 0.22f, kHaloB = 0.55f;
constexpr double kLogIntervalYears = 1.0e6;  // one CSV row per Myr of sim time

std::string make_log_path() {
  const std::time_t t = std::time(nullptr);
  char stamp[32];
  std::strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", std::localtime(&t));
  return std::string("logs/galaxy_") + stamp + ".csv";
}

double now_seconds() {
  return static_cast<double>(SDL_GetPerformanceCounter()) /
         static_cast<double>(SDL_GetPerformanceFrequency());
}

// A running simulation plus everything the renderer needs to draw it.
struct Session {
  std::unique_ptr<ISimulation> sim;  // null while the setup panel is showing
  ScenarioParams params;
  Vec2 view_center;
  double view_extent = 1.0;     // fits the visible bodies
  double dm_extent = 1.0;       // fits the dark-matter halo too
  std::vector<float> verts;  // regular bodies, highlighted ones, then halo; camera-relative
  RenderCounts counts;

  bool two_galaxies() const { return params.galaxy_b.has_value(); }
  bool distributed() const { return params.mass_model == MassModel::Distributed; }

  void start(const ScenarioParams& p, const ProgressFn& progress = {}) {
    auto s = make_simulation(p, progress);  // throws on invalid parameters; leaves the old session intact
    params = p;
    sim = std::move(s);
    view_center = sim->center_of_mass();  // fixed for the whole run
    view_extent = sim->initial_extent();
    dm_extent = sim->dark_matter_extent();
    counts = sim->render_counts();
    verts.assign(2 * (counts.regular + counts.highlighted + counts.halo), 0.0f);
  }
};

}  // namespace

int main(int argc, char** argv) {
  // --- setup state, optionally prefilled from the command line ---
  // Any command-line argument skips the setup panel and starts immediately.
  //   --two  --distributed  --stars <N>  --vel <km/s>  --speed <Myr/s>  --log
  //   --halo-frac <f>  --dm (draw the dark matter in blue)  --relax <Myr> (0 = off)
  SetupState setup;
  float speed_myr = kMinSpeedMyr;
  bool start_show_dm = false;
  bool stars_given = false;
  const bool auto_start = argc > 1;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--two") == 0) setup.two_galaxies = true;
    if (std::strcmp(argv[i], "--dm") == 0) start_show_dm = true;
    if (std::strcmp(argv[i], "--halo-frac") == 0 && i + 1 < argc) {
      const double f = std::atof(argv[++i]);
      if (f >= 0.0) setup.distributed.halo_mass_fraction = f;
    }
    if (std::strcmp(argv[i], "--relax") == 0 && i + 1 < argc) {
      const double r = std::atof(argv[++i]);
      if (r >= 0.0) setup.distributed.relax_myr = r;
    }
    if (std::strcmp(argv[i], "--log") == 0) setup.log_csv = true;
    if (std::strcmp(argv[i], "--distributed") == 0) setup.mass_model = MassModel::Distributed;
    if (std::strcmp(argv[i], "--stars") == 0 && i + 1 < argc) {
      const int n = std::atoi(argv[++i]);
      if (n > 0) {
        setup.a.num_stars = setup.b.num_stars = n;
        stars_given = true;
      }
    }
    if (std::strcmp(argv[i], "--vel") == 0 && i + 1 < argc) {
      const double v = std::atof(argv[++i]);
      if (v >= 0.0) setup.encounter.relative_speed_km_s = v;
    }
    if (std::strcmp(argv[i], "--speed") == 0 && i + 1 < argc) {
      const double v = std::atof(argv[++i]);
      if (v > 0.0) speed_myr = std::clamp(static_cast<float>(v), kMinSpeedMyr, kMaxSpeedMyr);
    }
  }

  // Distributed mode defaults to far fewer particles than the test-star mode.
  if (setup.mass_model == MassModel::Distributed && !stars_given)
    setup.a.num_stars = setup.b.num_stars = 5000;

  // --- window + GL context ---
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    std::fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
    return 1;
  }
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

  SDL_Window* window = SDL_CreateWindow(
      "Galaxy Sim", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720,
      SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
  if (!window) {
    std::fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
    return 1;
  }
  SDL_GLContext ctx = SDL_GL_CreateContext(window);
  if (!ctx) {
    std::fprintf(stderr, "SDL_GL_CreateContext failed: %s\n", SDL_GetError());
    return 1;
  }
  SDL_GL_SetSwapInterval(1);  // vsync

  const char* missing = nullptr;
  if (!gl::load(&missing)) {
    std::fprintf(stderr, "Missing OpenGL function: %s\n", missing);
    return 1;
  }
  PointRenderer renderer;
  {
    std::string err;
    if (!renderer.init(err)) {
      std::fprintf(stderr, "Renderer init failed: %s\n", err.c_str());
      return 1;
    }
  }

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui::GetIO().IniFilename = nullptr;
  ImGui::StyleColorsDark();
  ImGui_ImplSDL2_InitForOpenGL(window, ctx);
  ImGui_ImplOpenGL3_Init("#version 330 core");

  // --- run state ---
  Session session;
  bool paused = false;
  bool running = true;
  bool reset_requested = false;
  std::optional<ScenarioParams> pending_start;  // set inside a UI frame, run before the next one
  bool show_dark_matter = start_show_dm;  // draw the halo in blue (zooms out to fit it)

  double speed_acc_years = 0.0, speed_acc_wall = 0.0, actual_speed = 0.0;
  double energy_ratio = 1.0, energy_timer = 1e9;
  BoundFractions bound;  // refreshed together with E/E0
  double fps = 0.0;
  double last = now_seconds();

  // --- CSV logging (one row per kLogIntervalYears of simulated time) ---
  bool log_enabled = false;
  CsvLogger logger;
  double next_log_time = 0.0;
  std::string log_error;

  auto write_log_row = [&]() {
    const ISimulation& s = *session.sim;
    const BoundFractions bf = s.bound_fractions();
    std::optional<double> sep;
    if (session.two_galaxies()) sep = s.central_separation_ly();
    logger.write_row(s.time_years() / 1.0e6, sep, 100.0 * bf.system, 100.0 * bf.galaxy);
    next_log_time = (std::floor(s.time_years() / kLogIntervalYears) + 1.0) * kLogIntervalYears;
  };
  auto start_logging = [&]() {
    log_error.clear();
    if (!logger.open(make_log_path(), log_error)) {
      std::fprintf(stderr, "CSV logging disabled: %s\n", log_error.c_str());
      log_enabled = false;
      return;
    }
    write_log_row();  // row for the current time
  };
  auto stop_logging = [&]() { logger.close(); };

  // Draws a "relaxing" screen from inside the (long) constructor so the window
  // stays responsive. Closing the window aborts the build.
  auto show_progress = [&](double fraction) {
    SDL_Event ev;
    while (SDL_PollEvent(&ev)) {
      ImGui_ImplSDL2_ProcessEvent(&ev);
      if (ev.type == SDL_QUIT) {
        running = false;
        throw std::runtime_error("cancelled");
      }
    }
    int pw = 0, ph = 0;
    SDL_GL_GetDrawableSize(window, &pw, &ph);
    glViewport(0, 0, pw, ph);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x * 0.5f, ImGui::GetIO().DisplaySize.y * 0.5f),
                            ImGuiCond_Always, ImVec2(0.5f, 0.5f));
    ImGui::Begin("Relaxing", nullptr,
                 ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
                     ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoMove);
    ImGui::Text("Relaxing each galaxy so it starts in equilibrium...");
    ImGui::ProgressBar(static_cast<float>(fraction), ImVec2(360, 0));
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    SDL_GL_SwapWindow(window);
  };

  auto begin_session = [&](const ScenarioParams& p) {
    stop_logging();
    try {
      session.start(p, show_progress);
    } catch (const std::exception& e) {  // the panel validates, so this is a backstop
      std::fprintf(stderr, "Could not start: %s\n", e.what());
      session.sim.reset();
      return;
    }
    if (log_enabled) start_logging();
    paused = false;
    speed_acc_years = speed_acc_wall = actual_speed = 0.0;
    energy_ratio = 1.0;
    bound = BoundFractions{};
    energy_timer = 1e9;  // refresh immediately
    last = now_seconds();
  };

  if (auto_start) {
    log_enabled = setup.log_csv;
    begin_session(setup.to_params());
  }

  while (running) {
    ImGuiIO& io = ImGui::GetIO();

    SDL_Event e;
    while (SDL_PollEvent(&e)) {
      ImGui_ImplSDL2_ProcessEvent(&e);
      if (e.type == SDL_QUIT) running = false;
      if (e.type == SDL_KEYDOWN && !e.key.repeat) {
        if (e.key.keysym.sym == SDLK_ESCAPE) running = false;
        if (!io.WantCaptureKeyboard && session.sim) {
          if (e.key.keysym.sym == SDLK_SPACE) paused = !paused;
          if (e.key.keysym.sym == SDLK_r) reset_requested = true;
        }
      }
    }

    if (pending_start) {
      const ScenarioParams p = *pending_start;
      pending_start.reset();
      begin_session(p);
    } else if (reset_requested && session.sim) {
      begin_session(session.params);
    }
    reset_requested = false;

    const double t_now = now_seconds();
    const double frame_dt = std::min(t_now - last, 0.1);  // clamp after stalls
    last = t_now;
    fps = 0.9 * fps + 0.1 * (frame_dt > 0 ? 1.0 / frame_dt : 0.0);

    int w = 0, h = 0;
    SDL_GL_GetDrawableSize(window, &w, &h);
    if (w > 0 && h > 0) {
      glViewport(0, 0, w, h);
      glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT);
    }

    if (session.sim) {
      ISimulation& sim = *session.sim;
      const double speed = static_cast<double>(speed_myr) * 1.0e6;  // sim yr per wall s

      // --- advance within a wall-clock budget ---
      // One step per call: self-gravity steps can take tens of ms, so the budget
      // is checked after every step. Logged inside the loop, not per frame, so
      // no rows are skipped when the sim runs many Myr per frame.
      double advanced = 0.0;
      if (!paused) {
        sim.queue_time(speed * frame_dt);
        const double t0 = now_seconds();
        while (sim.backlog_years() >= sim.dt()) {
          advanced += sim.run_queued(1);
          if (logger.is_open() && sim.time_years() >= next_log_time) write_log_row();
          if (now_seconds() - t0 > kStepBudgetSeconds) break;
        }
        sim.drop_backlog();
      }
      speed_acc_years += advanced;
      speed_acc_wall += paused ? 0.0 : frame_dt;
      if (speed_acc_wall >= 0.5) {
        actual_speed = speed_acc_years / speed_acc_wall;
        speed_acc_years = speed_acc_wall = 0.0;
      }
      if (paused) actual_speed = 0.0;

      // E/E0 and the bound fractions are O(N) or worse; refresh a few times per second.
      energy_timer += frame_dt;
      if (energy_timer >= 0.25) {
        energy_ratio = sim.total_energy() / sim.initial_energy();
        bound = sim.bound_fractions();
        energy_timer = 0.0;
      }

      // --- draw the bodies ---
      if (w > 0 && h > 0) {
        sim.fill_render_points(session.view_center, paused ? 0.0 : sim.backlog_years(),
                               session.verts.data());
        const bool dm = show_dark_matter && session.counts.halo > 0;
        const Camera cam = fit_camera(session.view_center,
                                      dm ? session.dm_extent : session.view_extent, w, h);
        const Vec2 scale = ndc_scale(cam, w, h);
        const auto sx = static_cast<float>(scale.x), sy = static_cast<float>(scale.y);
        const std::size_t n_vis = session.counts.regular + session.counts.highlighted;

        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);  // additive
        glPointSize(1.0f);
        renderer.upload(session.verts.data(), n_vis + (dm ? session.counts.halo : 0));
        // Halo first, so the stars and black holes brighten on top of the haze.
        if (dm) renderer.draw(n_vis, session.counts.halo, sx, sy, kHaloR, kHaloG, kHaloB);
        renderer.draw(0, session.counts.regular, sx, sy, kStarIntensity, kStarIntensity,
                      kStarIntensity);
        renderer.draw(session.counts.regular, session.counts.highlighted, sx, sy,
                      kCentralIntensity, kCentralIntensity, kCentralIntensity);
      }
    }

    // --- UI ---
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();

    if (!session.sim) {
      if (draw_setup_panel(setup)) {
        log_enabled = setup.log_csv;
        pending_start = setup.to_params();
      }
    } else {
      ISimulation& sim = *session.sim;
      bool back_to_setup = false;

      ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Always);
      ImGui::SetNextWindowBgAlpha(0.55f);
      ImGui::Begin("HUD", nullptr,
                   ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
                       ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
                       ImGuiWindowFlags_NoMove);
      ImGui::Text("%s", session.distributed() ? "distributed mass + Barnes-Hut"
                                               : "centralized mass + Kepler");
      ImGui::Text("t = %.1f Myr%s", sim.time_years() / 1.0e6, paused ? "  [PAUSED]" : "");
      ImGui::SetNextItemWidth(220);
      ImGui::SliderFloat("Speed (Myr/s)", &speed_myr, kMinSpeedMyr, kMaxSpeedMyr, "%.1f",
                         ImGuiSliderFlags_Logarithmic);
      ImGui::Text("actual: %.1f Myr/s", actual_speed / 1.0e6);
      if (std::abs(sim.initial_energy()) > 1e-30) {
        ImGui::Text("E / E0 = %.10f  (1 %+.2e)", energy_ratio, energy_ratio - 1.0);
      } else {
        ImGui::Text("E / E0 = n/a (E0 ~ 0)");
      }
      ImGui::Text("bound to system: %.3f%%", 100.0 * bound.system);
      if (!std::isnan(bound.galaxy)) ImGui::Text("bound to a galaxy: %.3f%%", 100.0 * bound.galaxy);
      if (session.distributed()) {
        ImGui::Text("bodies: %zu disk + %zu halo + %zu black holes", session.counts.regular,
                    session.counts.halo, session.counts.highlighted);
        if (session.counts.halo > 0) ImGui::Checkbox("Show dark matter (blue)", &show_dark_matter);
      } else {
        ImGui::Text("bodies: %zu stars + %zu central", session.counts.regular,
                    session.counts.highlighted);
      }
      ImGui::Text("dt = %.0f yr   %.0f fps", sim.dt(), fps);

      if (ImGui::Button(paused ? "Resume" : "Pause")) paused = !paused;
      ImGui::SameLine();
      if (ImGui::Button("Reset")) reset_requested = true;
      ImGui::SameLine();
      if (ImGui::Button("New setup")) back_to_setup = true;

      const bool was_logging = log_enabled;
      ImGui::Checkbox("Log CSV (1 row / Myr)", &log_enabled);
      if (log_enabled != was_logging) {
        if (log_enabled) start_logging(); else stop_logging();
      }
      if (logger.is_open()) ImGui::TextDisabled("-> %s", logger.path().c_str());
      if (!log_error.empty()) ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "%s", log_error.c_str());

      ImGui::TextDisabled("Space: pause   R: reset   Esc: quit");
      ImGui::End();

      if (back_to_setup) {
        stop_logging();
        session.sim.reset();
      }
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    SDL_GL_SwapWindow(window);
  }

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();
  renderer.shutdown();
  SDL_GL_DeleteContext(ctx);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
