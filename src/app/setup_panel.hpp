#pragma once

#include "sim/params.hpp"

namespace galaxy {

// Everything the setup panel edits.
struct SetupState {
  MassModel mass_model = MassModel::Centralized;
  bool two_galaxies = false;
  GalaxyParams a;
  GalaxyParams b;
  EncounterParams encounter;
  DistributedParams distributed = [] {
    DistributedParams d;
    d.relax_myr = kDefaultRelaxMyr;
    return d;
  }();
  bool log_csv = false;  // write separation and bound % to a CSV every 1 Myr

  ScenarioParams to_params() const {
    ScenarioParams p;
    p.mass_model = mass_model;
    p.galaxy_a = a;
    if (two_galaxies) p.galaxy_b = b;
    p.encounter = encounter;
    p.distributed = distributed;
    return p;
  }
};

// Draws the centered setup window (call between ImGui::NewFrame and Render).
// Returns true on the frame the user presses Start with valid settings.
bool draw_setup_panel(SetupState& state);

}  // namespace galaxy
