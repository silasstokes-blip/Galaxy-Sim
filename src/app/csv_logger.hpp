#pragma once

#include <fstream>
#include <optional>
#include <string>

namespace galaxy {

// Appends rows of simulation diagnostics to a CSV file:
//   time_myr,separation_ly,bound_system_pct,bound_galaxy_pct
// The separation field is left empty for single-galaxy runs, and the galaxy
// field is left empty when passed NaN (distributed-mass mode).
class CsvLogger {
 public:
  // Creates parent directories as needed and writes the header. On failure
  // returns false and fills `error`.
  bool open(const std::string& path, std::string& error);
  void close();
  bool is_open() const { return out_.is_open(); }
  const std::string& path() const { return path_; }

  // Flushes every row so the file is intact even if the app is killed.
  void write_row(double time_myr, std::optional<double> separation_ly,
                 double bound_system_pct, double bound_galaxy_pct);

 private:
  std::ofstream out_;
  std::string path_;
};

}  // namespace galaxy
