#include <fmt/core.h>
#include <fmt/ostream.h>

#include "nlohmann/json.hpp"

#include "TObject.h"

void print_hline() {
  std::cout << "-----------------------------------------------------------------------------------"
               "----------------------------------------------------------";
}
void vcs_summary() {
  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("database/status.json");
    json_input_file >> j;
  }

  // std::cout << j.dump(2);
  auto print_header = []() {
    // print header
    std::cout << "\n";
    fmt::print(" {:>9} ", "Setting");
    fmt::print(" {:^7} ", "target");
    fmt::print(" {:>6} ", "P_hms");
    fmt::print(" {:>7} ", "th_hms");
    fmt::print(" {:>6} ", "P_shms");
    fmt::print(" {:>7} ", "th_shms");
    fmt::print(" {:>5} ", "n_runs");
    fmt::print(" {:>7} ", "n_pi0");
    fmt::print(" {:^20} ", "pi0 yield / mC");
    fmt::print(" {:>9} ", "est. n_VCS");
    fmt::print(" {:>7} ", "cut1");
    fmt::print(" {:>7} ", "cut2");
    fmt::print(" {:>7} ", "cut3");
    fmt::print(" {:>7} ", "Q [mC]");
    std::cout << "\n";
  };

  std::string old_name = "";

  print_header();
  print_hline();
  for (auto setting : j) {
    const std::string name      = setting["name"];
    const std::string target    = setting["target"];
    const double      p_hms     = setting["hms_momentum"].get<double>();
    const double      th_hms    = setting["hms_angle"].get<double>();
    const double      p_shms    = setting["shms_momentum"].get<double>();
    const double      th_shms   = setting["shms_angle"].get<double>();
    const int         n_runs    = setting["good_run_list"].size();
    const double      n_pi0     = setting["pi0_count"].get<double>();
    const double      n_vcs     = setting["projections"]["vcs_count"].get<double>();
    const int         n_cut1    = setting["projections"]["vcs_count_cut1"].get<double>();
    const int         n_cut2    = setting["projections"]["vcs_count_cut2"].get<double>();
    const int         n_cut3    = setting["projections"]["vcs_count_cut3"].get<double>();
    const double      charge    = setting["total_charge"].get<double>();
    const double      yield     = (charge > 0) ? n_pi0 / charge : 0;
    const double      yield_err = (charge > 0) ? sqrt(n_pi0) / charge : 0;

    if (name != old_name) {
      std::cout << "\n";
    }
    old_name = name;

    fmt::print(" {:>9} ", name);
    fmt::print(" {:^7} ", target);
    fmt::print(" {:>6.3f} ", p_hms);
    fmt::print(" {:>7.2f} ", th_hms);
    fmt::print(" {:>6.3f} ", p_shms);
    fmt::print(" {:>7.2f} ", th_shms);
    fmt::print(" {:>5} ", n_runs);
    fmt::print(" {:>1.2e} ", n_pi0);
    fmt::print(" {:>1.2e} +- {:<1.2e} ", yield, yield_err);
    if (n_vcs > 0) {
      fmt::print("  {:>1.2e} ", n_vcs);
      fmt::print(" {:>7} ", n_cut1);
      fmt::print(" {:>7} ", n_cut2);
      fmt::print(" {:>7} ", n_cut3);
    } else {
      fmt::print(" {:>9} ", "");
      fmt::print(" {:>7} ", "");
      fmt::print(" {:>7} ", "");
      fmt::print(" {:>7} ", "");
    }
    fmt::print(" {:>7.0f}", charge);
    std::cout << "\n";
  }
  std::cout << std::endl;
}
