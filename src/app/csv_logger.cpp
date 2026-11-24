#include "app/csv_logger.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <system_error>

namespace galaxy {

bool CsvLogger::open(const std::string& path, std::string& error) {
  close();

  const std::filesystem::path p(path);
  if (p.has_parent_path()) {
    std::error_code ec;
    std::filesystem::create_directories(p.parent_path(), ec);
    if (ec) {
      error = "cannot create directory " + p.parent_path().string() + ": " + ec.message();
      return false;
    }
  }

  out_.open(path, std::ios::out | std::ios::trunc);
  if (!out_) {
    error = "cannot open " + path + " for writing";
    return false;
  }
  path_ = path;
  out_ << "time_myr,separation_ly,bound_system_pct,bound_galaxy_pct\n";
  out_.flush();
  return true;
}

void CsvLogger::close() {
  if (out_.is_open()) out_.close();
  path_.clear();
}

void CsvLogger::write_row(double time_myr, std::optional<double> separation_ly,
                          double bound_system_pct, double bound_galaxy_pct) {
  if (!out_.is_open()) return;
  // Fields that do not apply (single galaxy: no separation; distributed mode:
  // no per-galaxy bound fraction, passed as NaN) are left empty.
  char sep[48] = "", gal[48] = "";
  if (separation_ly) std::snprintf(sep, sizeof sep, "%.3f", *separation_ly);
  if (!std::isnan(bound_galaxy_pct)) std::snprintf(gal, sizeof gal, "%.6f", bound_galaxy_pct);
  char buf[160];
  std::snprintf(buf, sizeof buf, "%.4f,%s,%.6f,%s\n", time_myr, sep, bound_system_pct, gal);
  out_ << buf;
  out_.flush();
}

}  // namespace galaxy
