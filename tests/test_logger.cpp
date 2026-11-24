#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "app/csv_logger.hpp"

using namespace galaxy;

static int g_failures = 0;

#define CHECK(cond)                                                   \
  do {                                                                \
    if (!(cond)) {                                                    \
      std::printf("    FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      ++g_failures;                                                   \
    }                                                                 \
  } while (0)

static std::vector<std::string> read_lines(const std::string& path) {
  std::ifstream in(path);
  std::vector<std::string> lines;
  for (std::string l; std::getline(in, l);) lines.push_back(l);
  return lines;
}

int main() {
  const std::string path = "logger_test_out/nested/test.csv";  // exercises directory creation
  std::filesystem::remove_all("logger_test_out");

  CsvLogger log;
  std::string err;
  CHECK(!log.is_open());
  CHECK(log.open(path, err));
  CHECK(log.is_open());
  log.write_row(0.0, 150000.0, 100.0, 100.0);
  log.write_row(1.0, std::nullopt, 99.5, 98.25);
  log.write_row(2.0, 1234.5, 80.0, std::numeric_limits<double>::quiet_NaN());
  // Rows must be on disk before close() (flushed per row).
  CHECK(read_lines(path).size() == 4);
  log.close();
  CHECK(!log.is_open());

  const auto lines = read_lines(path);
  CHECK(lines.size() == 4);
  if (lines.size() == 4) {
    CHECK(lines[0] == "time_myr,separation_ly,bound_system_pct,bound_galaxy_pct");
    CHECK(lines[1] == "0.0000,150000.000,100.000000,100.000000");
    CHECK(lines[2] == "1.0000,,99.500000,98.250000");  // separation left empty
    CHECK(lines[3] == "2.0000,1234.500,80.000000,");   // NaN galaxy fraction left empty
  }

  // A path that cannot be created reports an error instead of crashing.
  CsvLogger bad;
  std::string err2;
  CHECK(!bad.open("logger_test_out/nested/test.csv/impossible/x.csv", err2));
  CHECK(!err2.empty());

  std::filesystem::remove_all("logger_test_out");
  std::printf("%s (%d failed check%s)\n", g_failures ? "FAILED" : "ALL PASSED", g_failures,
              g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
