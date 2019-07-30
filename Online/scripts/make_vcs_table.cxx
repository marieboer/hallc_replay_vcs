#include <fmt/core.h>
#include <fmt/ostream.h>

#include "nlohmann/json.hpp"

#include "TObject.h"

// CHANGE THIS ONCE THE EXPERIMENT STARTS!
constexpr const int FIRST_VCS_RUN = 8819;

std::string setting_name(const nlohmann::json& spectro_info, const std::string& target,
                         const nlohmann::json& settings) {
  // loop over all settings and check if we have a setting for this run
  for (const auto& setting : settings) {
    if (fabs(setting["hms_momentum"].get<double>() - spectro_info["hms_momentum"].get<double>()) <
            0.015 &&
        fabs(setting["hms_angle"].get<double>() - spectro_info["hms_angle"].get<double>()) < 0.1 &&
        fabs(setting["shms_momentum"].get<double>() - spectro_info["shms_momentum"].get<double>()) <
            0.015 &&
        fabs(setting["shms_angle"].get<double>() - spectro_info["shms_angle"].get<double>()) <
            0.1 &&
        target == setting["target"]) {
      return setting["name"];
    }
  }
  // nothing found
  return "other";
}

void make_vcs_table() {
  using nlohmann::json;
  json vcsdb;
  {
    std::ifstream json_input_file("database/VCSdb.json");
    json_input_file >> vcsdb;
  }
  json rundb;
  {
    std::ifstream json_input_file("database/rundb_coin.json");
    json_input_file >> rundb;
  }
  json countdb_hms;
  {
    std::ifstream json_input_file("database/countdb_hms.json");
    json_input_file >> countdb_hms;
  }
  json countdb_shms;
  {
    std::ifstream json_input_file("database/countdb_shms.json");
    json_input_file >> countdb_shms;
  }
  json commentdb;
  {
    std::ifstream json_input_file("database/commentdb.json");
    json_input_file >> commentdb;
  }
  json settings;
  {
    std::ifstream json_input_file("database/settings.json");
    json_input_file >> settings;
  }
  std::ofstream good_runfile("database/good_runs.txt");

  auto print_header = []() {
    std::cout << "\n";
    fmt::print(" {:<5} ", "Run");
    fmt::print(" {:^8} ", "Setting");
    fmt::print(" {:^5} ", "Target");
    fmt::print(" {:>7} ", "P_hms ");
    fmt::print(" {:<7} ", "th_hms");
    fmt::print(" {:>7} ", "P_shms ");
    fmt::print(" {:<7} ", "th_shms");
    fmt::print(" {:^8} ", "start");
    fmt::print(" {:^17} ", "end time");
    fmt::print(" {:^18} ", "HMS p yield");
    fmt::print(" {:^18} ", "SHMS e yield");
    fmt::print(" {:^10} ", "pi0 peak");
    // fmt::print(" {:^8} ", "peak2");
    // fmt::print(" {:^8} ", "peak3");
    fmt::print(" {:>5} ", "Q(mC)");
    fmt::print(" {:<} ", "comment");
    std::cout << "\n";
  };

  std::string old_target = "";
  double      p_hms      = 0.0;
  double      th_hms     = 0.0;
  double      p_shms     = 0.0;
  double      th_shms    = 0.0;

  for (json::iterator it = rundb.begin(); it != rundb.end(); ++it) {
    auto runjs = it.value();

    if (std::stoi(it.key()) < FIRST_VCS_RUN) {
      continue;
    }

    std::string target_lab = runjs["target"]["target_label"].get<std::string>();
    if (target_lab != old_target ||
        ((p_hms != runjs["spectrometers"]["hms_momentum"].get<double>()) ||
         (th_hms != runjs["spectrometers"]["hms_angle"].get<double>()) ||
         (p_shms != runjs["spectrometers"]["shms_momentum"].get<double>()) ||
         (th_shms != runjs["spectrometers"]["shms_angle"].get<double>()))) {
      print_header();
    }

    p_hms     = runjs["spectrometers"]["hms_momentum"].get<double>();
    th_hms    = runjs["spectrometers"]["hms_angle"].get<double>();
    p_shms    = runjs["spectrometers"]["shms_momentum"].get<double>();
    th_shms   = runjs["spectrometers"]["shms_angle"].get<double>();
    auto name = setting_name(runjs["spectrometers"], target_lab, settings);

    old_target = target_lab;

    good_runfile << std::stoi(it.key()) << "\n";

    fmt::print(" {:<5} ", std::stoi(it.key()));
    fmt::print(" {:^8} ", name);
    fmt::print(" {:^5} ", target_lab);
    fmt::print(" {:>7.3f} ", runjs["spectrometers"]["hms_momentum"].get<double>());
    fmt::print(" {:<7.2f} ", runjs["spectrometers"]["hms_angle"].get<double>());
    fmt::print(" {:>7.3f} ", runjs["spectrometers"]["shms_momentum"].get<double>());
    fmt::print(" {:<7.2f} ", runjs["spectrometers"]["shms_angle"].get<double>());
    if (runjs["run_info"].find("start_time") != runjs["run_info"].end()) {
      auto start_time = runjs["run_info"]["start_time"].get<std::string>();
      fmt::print(" {:^8} ", start_time.substr(0, start_time.size() - 13));
    } else {
      fmt::print(" {:8} ", "");
    }
    if (runjs["run_info"].find("stop_time") != runjs["run_info"].end()) {
      auto stop_time = runjs["run_info"]["stop_time"].get<std::string>();
      fmt::print(" {:^17} ", stop_time.substr(0, stop_time.size() - 4));
    } else {
      fmt::print(" {:17} ", "");
    }

    double charge = -1;
    if (countdb_hms.count(it.key()) != 0) {
      auto   rl_hms    = countdb_hms[it.key()];
      double n_hms     = rl_hms["count_e"].get<double>();
      double ps_factor = rl_hms["ps_factor"].get<double>();
      charge           = rl_hms["good_total_charge"].get<double>();
      double hms_yield = n_hms * ps_factor / charge;
      double hms_unc   = sqrt(n_hms) * ps_factor / charge;
      fmt::print(" {:>1.2e}", hms_yield);
      fmt::print(" ± {:<1.2e} ", hms_unc);
    } else {
      fmt::print(" {:>1.2e}", 0.);
      fmt::print(" ± {:<1.2e} ", 0.);
    }
    if (countdb_shms.count(it.key()) != 0) {
      auto   rl_shms    = countdb_shms[it.key()];
      double n_shms     = rl_shms["count_h"].get<double>();
      double ps_factor  = rl_shms["ps_factor"].get<double>();
      charge            = rl_shms["good_total_charge"].get<double>();
      double shms_yield = n_shms * ps_factor / charge;
      double shms_unc   = sqrt(n_shms) * ps_factor / charge;
      fmt::print(" {:>1.2e}", shms_yield);
      fmt::print(" ± {:<1.2e} ", shms_unc);
    } else {
      fmt::print(" {:>1.2e}", 0.);
      fmt::print(" ± {:<1.2e} ", 0.);
    }
    if (vcsdb.count(it.key()) != 0) {
      try {
        double n_peak1 = vcsdb[it.key()]["missing_mass"]["peak1"]["integral"].get<double>();
        //        double n_peak2 =
        //        vcsdb[it.key()]["missing_mass"]["peak2"]["integral"].get<double>(); double
        //        n_peak3 = vcsdb[it.key()]["missing_mass"]["peak3"]["integral"].get<double>();
        if (charge > 0) {
          fmt::print(" {:>1.2e}  ", (n_peak1 > 1e-100) ? n_peak1 / charge : 0);
          fmt::print(" {:>5.0f} ", charge);
        } else {
          fmt::print(" {:>9} ", "");
          fmt::print(" {:>5} ", "");
        }
      } catch (std::domain_error) {
        ;  // do nothing
      } catch (nlohmann::detail::type_error) {
        // should never happen...
        std::cout << "WEIRD TYPE ERROR" << std::endl;
      }
    } else {
      fmt::print(" {:>9} ", "");
      fmt::print(" {:>5} ", "");
    }
    std::string comment;
    if (commentdb.count(it.key()) != 0) {
      if (commentdb[it.key()].find("comment") != commentdb[it.key()].end()) {
        try {
          comment = commentdb[it.key()]["comment"].get<std::string>();
        } catch (std::domain_error) {
          ;  // do nothing
        }
      }
    }
    fmt::print(" {:<} ", comment);
    std::cout << "\n";
  }
  print_header();
}
