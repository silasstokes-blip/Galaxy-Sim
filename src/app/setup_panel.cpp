#include "app/setup_panel.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

#include "imgui.h"
#include "sim/units.hpp"

namespace galaxy {

namespace {

constexpr float kItemWidth = 240.0f;

void galaxy_controls(GalaxyParams& g, const DistributedParams& d, bool distributed) {
  static const double kMassMin = 1.0e9, kMassMax = 1.0e12;
  static const double kRadiusMin = 5.0e3, kRadiusMax = 2.0e5;

  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::SliderScalar("Total mass (Msun)", ImGuiDataType_Double, &g.total_mass_msun,
                      &kMassMin, &kMassMax, "%.3e", ImGuiSliderFlags_Logarithmic);
  ImGui::SetNextItemWidth(kItemWidth);
  if (distributed) {
    ImGui::SliderInt("Disk particles", &g.num_stars, 1000, 100000, "%d",
                     ImGuiSliderFlags_Logarithmic);
  } else {
    ImGui::SliderInt("Stars", &g.num_stars, 1000, 500000, "%d", ImGuiSliderFlags_Logarithmic);
  }
  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::SliderScalar("Radius (ly)", ImGuiDataType_Double, &g.radius_ly, &kRadiusMin,
                      &kRadiusMax, "%.0f");
  ImGui::Checkbox("Counter-clockwise spin", &g.counterclockwise);

  if (distributed) {
    const double m_bh = d.bh_mass_fraction * g.total_mass_msun;
    const double m_halo = d.halo_mass_fraction * g.total_mass_msun;
    const double m_disk = g.total_mass_msun - m_bh - m_halo;
    const double n_halo = d.halo_mass_fraction > 0.0 ? g.num_stars * d.halo_particles_per_disk_particle : 0.0;
    ImGui::TextDisabled("Black hole %.2e   disk %.2e (%.2e each)", m_bh, m_disk,
                        g.num_stars > 0 ? m_disk / g.num_stars : 0.0);
    if (n_halo > 0.0)
      ImGui::TextDisabled("Halo %.2e Msun in %.0f particles (%.2e each)", m_halo, n_halo,
                          m_halo / n_halo);
  } else {
    double hole_pct = 100.0 * g.inner_fraction;
    const double hole_min = 0.0, hole_max = 60.0;
    ImGui::SetNextItemWidth(kItemWidth);
    if (ImGui::SliderScalar("Central hole (% of radius)", ImGuiDataType_Double, &hole_pct, &hole_min,
                            &hole_max, "%.1f"))
      g.inner_fraction = hole_pct / 100.0;
    ImGui::TextDisabled("Central body: %.4e Msun (total - stars)", central_mass_msun(g));
  }
}

void distributed_controls(DistributedParams& d, const GalaxyParams& a, const GalaxyParams* b) {
  // Black hole fraction, edited as a percent (0 allowed: particles only).
  double bh_pct = 100.0 * d.bh_mass_fraction;
  const double pct_min = 0.0, pct_max = 10.0;
  ImGui::SetNextItemWidth(kItemWidth);
  if (ImGui::SliderScalar("Black hole mass (% of galaxy)", ImGuiDataType_Double, &bh_pct, &pct_min,
                          &pct_max, "%.3f", ImGuiSliderFlags_Logarithmic))
    d.bh_mass_fraction = bh_pct / 100.0;

  int profile = d.profile == DiskProfile::Exponential ? 1 : 0;
  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::Combo("Radial profile", &profile, "Uniform\0Exponential\0");
  d.profile = profile == 1 ? DiskProfile::Exponential : DiskProfile::Uniform;

  if (d.profile == DiskProfile::Exponential) {
    static const double kLenMin = 0.05, kLenMax = 1.0;
    ImGui::SetNextItemWidth(kItemWidth);
    ImGui::SliderScalar("Scale length (x radius)", ImGuiDataType_Double, &d.scale_length_fraction,
                        &kLenMin, &kLenMax, "%.2f");
  }

  // Dark-matter halo (live particles that stabilize the disk).
  double halo_pct = 100.0 * d.halo_mass_fraction;
  const double halo_min = 0.0, halo_max = 95.0;
  ImGui::SetNextItemWidth(kItemWidth);
  if (ImGui::SliderScalar("Dark matter halo (% of galaxy)", ImGuiDataType_Double, &halo_pct,
                          &halo_min, &halo_max, "%.0f"))
    d.halo_mass_fraction = halo_pct / 100.0;
  if (d.halo_mass_fraction > 0.0) {
    static const double kRadMin = 0.5, kRadMax = 5.0, kCoreMin = 0.05, kCoreMax = 1.0;
    static const double kRatioMin = 0.5, kRatioMax = 8.0;
    ImGui::SetNextItemWidth(kItemWidth);
    ImGui::SliderScalar("Halo radius (x galaxy radius)", ImGuiDataType_Double,
                        &d.halo_radius_fraction, &kRadMin, &kRadMax, "%.2f");
    ImGui::SetNextItemWidth(kItemWidth);
    ImGui::SliderScalar("Halo core (x galaxy radius)", ImGuiDataType_Double,
                        &d.halo_core_fraction, &kCoreMin, &kCoreMax, "%.2f");
    ImGui::SetNextItemWidth(kItemWidth);
    ImGui::SliderScalar("Halo particles per disk particle", ImGuiDataType_Double,
                        &d.halo_particles_per_disk_particle, &kRatioMin, &kRatioMax, "%.1f",
                        ImGuiSliderFlags_Logarithmic);
  } else {
    ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f),
                       "No halo: the disk is gravitationally unstable and will fragment.");
  }

  // Dispersion, edited as a percent of the local circular speed; 0 = cold disk.
  double hole_pct = 100.0 * d.inner_hole_fraction;
  const double hole_min = 0.0, hole_max = 60.0;
  ImGui::SetNextItemWidth(kItemWidth);
  if (ImGui::SliderScalar("Central hole (% of radius)", ImGuiDataType_Double, &hole_pct, &hole_min,
                          &hole_max, "%.0f"))
    d.inner_hole_fraction = hole_pct / 100.0;
  ImGui::TextDisabled("Empties the centre of disk and halo; the fastest orbits set dt, so a bigger hole allows a larger timestep.");

  double disp_pct = 100.0 * d.dispersion_fraction;
  const double disp_min = 0.0, disp_max = 50.0;
  ImGui::SetNextItemWidth(kItemWidth);
  if (ImGui::SliderScalar("Velocity dispersion (% of circular)", ImGuiDataType_Double, &disp_pct,
                          &disp_min, &disp_max, "%.1f"))
    d.dispersion_fraction = disp_pct / 100.0;

  static const double kThetaMin = 0.1, kThetaMax = 1.0;
  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::SliderScalar("Opening angle (theta)", ImGuiDataType_Double, &d.opening_angle, &kThetaMin,
                      &kThetaMax, "%.2f", ImGuiSliderFlags_Logarithmic);

  static const double kRelaxMin = 0.0, kRelaxMax = 2000.0;
  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::SliderScalar("Pre-relax before start (Myr)", ImGuiDataType_Double, &d.relax_myr,
                      &kRelaxMin, &kRelaxMax, d.relax_myr > 0.0 ? "%.0f" : "off");
  ImGui::TextDisabled(d.relax_myr > 0.0
                          ? "Each galaxy settles alone first (takes a while); ~700 Myr removes the initial breathing."
                          : "Off: the raw initial state collapses and rebounds for ~700 Myr.");

  // Softening is derived from the particle spacing; show what it works out to.
  auto spacing = [&d](const GalaxyParams& g) {
    double s = g.num_stars > 0 ? g.radius_ly * std::sqrt(std::numbers::pi / g.num_stars) : 0.0;
    if (d.halo_mass_fraction > 0.0 && g.num_stars > 0) {
      const double n_halo = g.num_stars * d.halo_particles_per_disk_particle;
      s = std::max(s, d.halo_radius_fraction * g.radius_ly * std::sqrt(std::numbers::pi / n_halo));
    }
    return s;
  };
  double s = spacing(a);
  if (b) s = std::max(s, spacing(*b));
  ImGui::TextDisabled("Softening: %.0f ly (mean particle spacing %.0f ly)",
                      d.softening_spacing_factor * s, s);
}

void encounter_controls(EncounterParams& e) {
  static const double kSepMin = 2.0e4, kSepMax = 1.0e6;
  static const double kVelMin = 0.0, kVelMax = 1000.0;
  static const double kImpMin = -2.0e5, kImpMax = 2.0e5;

  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::SliderScalar("Starting distance (ly)", ImGuiDataType_Double, &e.separation_ly,
                      &kSepMin, &kSepMax, "%.0f", ImGuiSliderFlags_Logarithmic);
  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::SliderScalar("Relative speed (km/s)", ImGuiDataType_Double, &e.relative_speed_km_s,
                      &kVelMin, &kVelMax, "%.0f");
  ImGui::SetNextItemWidth(kItemWidth);
  ImGui::SliderScalar("Impact parameter (ly)", ImGuiDataType_Double, &e.impact_parameter_ly,
                      &kImpMin, &kImpMax, "%.0f");
}

}  // namespace

bool draw_setup_panel(SetupState& s) {
  const ImVec2 center = ImGui::GetIO().DisplaySize;
  ImGui::SetNextWindowPos(ImVec2(center.x * 0.5f, center.y * 0.5f), ImGuiCond_Always,
                          ImVec2(0.5f, 0.5f));
  ImGui::Begin("Galaxy Sim: setup", nullptr,
               ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);

  // --- mass model ---
  int model = s.mass_model == MassModel::Distributed ? 1 : 0;
  ImGui::RadioButton("Centralized mass + Kepler", &model, 0);
  ImGui::SameLine();
  ImGui::RadioButton("Distributed mass + Barnes-Hut", &model, 1);
  const MassModel previous_model = s.mass_model;
  s.mass_model = model == 1 ? MassModel::Distributed : MassModel::Centralized;
  const bool distributed = s.mass_model == MassModel::Distributed;
  if (s.mass_model != previous_model) {
    // Each mode has its own sensible particle count: test stars are cheap, but
    // self-gravitating bodies (disk + halo) cost O(N log N) per step.
    s.a.num_stars = s.b.num_stars = distributed ? 5000 : 50000;
  }
  ImGui::TextDisabled(distributed
                          ? "A black hole, a disk of stars and a live dark-matter halo, all attracting each other."
                          : "One heavy central body; 1 Msun stars are test particles (O(N)).");

  // --- galaxy count ---
  int count = s.two_galaxies ? 2 : 1;
  ImGui::RadioButton("One galaxy", &count, 1);
  ImGui::SameLine();
  ImGui::RadioButton("Two galaxies", &count, 2);
  s.two_galaxies = (count == 2);

  ImGui::SeparatorText(s.two_galaxies ? "Galaxy 1" : "Galaxy");
  ImGui::PushID("galaxy1");
  galaxy_controls(s.a, s.distributed, distributed);
  ImGui::PopID();

  if (s.two_galaxies) {
    ImGui::SeparatorText("Galaxy 2");
    ImGui::PushID("galaxy2");
    galaxy_controls(s.b, s.distributed, distributed);
    ImGui::PopID();
  }

  if (distributed) {
    ImGui::SeparatorText("Mass distribution");
    distributed_controls(s.distributed, s.a, s.two_galaxies ? &s.b : nullptr);
    const double per_disk_particle =
        1.0 + (s.distributed.halo_mass_fraction > 0.0 ? s.distributed.halo_particles_per_disk_particle : 0.0);
    const double bodies = per_disk_particle * (s.a.num_stars + (s.two_galaxies ? s.b.num_stars : 0));
    ImGui::TextDisabled("Total simulated bodies: %.0f", bodies);
    if (bodies > 50000)
      ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f),
                         "Heavy: self-gravity costs O(N log N) per step; expect slower than real time.");
  }

  if (s.two_galaxies) {
    ImGui::SeparatorText("Encounter");
    encounter_controls(s.encounter);

    // Is the starting orbit bound? Compare v^2 with the two-body escape value.
    const double m_tot = distributed ? s.a.total_mass_msun + s.b.total_mass_msun
                                     : central_mass_msun(s.a) + central_mass_msun(s.b);
    if (m_tot > 0.0 && s.encounter.separation_ly > 0.0) {
      const double v_esc_kms = std::sqrt(2.0 * units::kG * m_tot / s.encounter.separation_ly) /
                               units::kKmPerSecToLyPerYr;
      const bool bound = s.encounter.relative_speed_km_s < v_esc_kms;
      ImGui::TextDisabled("Escape speed at this distance: %.0f km/s (%s)", v_esc_kms,
                          bound ? "bound orbit" : "one-time flyby");
    }
    if (s.encounter.separation_ly < s.a.radius_ly + s.b.radius_ly)
      ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f),
                         "Warning: the disks overlap at the start.");
  }

  const std::string error = validate_scenario(s.to_params());
  if (!error.empty()) ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "%s", error.c_str());

  ImGui::Spacing();
  ImGui::Checkbox("Log separation and bound % to CSV (every 1 Myr)", &s.log_csv);
  ImGui::BeginDisabled(!error.empty());
  const bool start = ImGui::Button("Start", ImVec2(140, 0));
  ImGui::EndDisabled();

  ImGui::End();
  return start && error.empty();
}

}  // namespace galaxy
