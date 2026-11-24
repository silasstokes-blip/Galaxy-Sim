#pragma once

// Unit system used everywhere in the simulation core:
//   length = light-year, time = Julian year, mass = solar mass.
// A nice side effect: 1 ly/yr == c exactly, so velocities are in units of c.

namespace galaxy::units {

// Defining constants (SI).
inline constexpr double kGMSun_SI      = 1.32712440018e20;   // m^3 / s^2   (IAU nominal-ish)
inline constexpr double kLightYear_m   = 9.4607304725808e15; // m
inline constexpr double kJulianYear_s  = 3.15576e7;           // s
inline constexpr double kLightSpeed_kms = 299792.458;         // km / s

// Gravitational constant in ly^3 / (Msun * yr^2).  (~1.56e-13)
inline constexpr double kG = kGMSun_SI * kJulianYear_s * kJulianYear_s /
                             (kLightYear_m * kLightYear_m * kLightYear_m);

// Convert km/s -> ly/yr.
inline constexpr double kKmPerSecToLyPerYr = 1.0 / kLightSpeed_kms;

// Every star is exactly one solar mass.
inline constexpr double kStarMassMsun = 1.0;

}  // namespace galaxy::units
