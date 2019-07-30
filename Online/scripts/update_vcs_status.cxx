#include <cmath>
#include <iostream>
#include "nlohmann/json.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>

void update_vcs_status() {
  std::cout << "Analyzing VCS run info...\n";

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
  json jsum;
  {
    std::ifstream json_input_file("database/settings.json");
    json_input_file >> jsum;
  }

  // set all counters to zero
  for (json::iterator iset = jsum.begin(); iset != jsum.end(); ++iset) {
    auto& sumjs                            = iset.value();
    sumjs["pi0_count"]                     = 0;
    sumjs["projections"]["vcs_count"]      = 0;
    sumjs["projections"]["vcs_count_cut1"] = 0;
    sumjs["projections"]["vcs_count_cut2"] = 0;
    sumjs["projections"]["vcs_count_cut3"] = 0;
    sumjs["total_charge"]                  = 0.;
    std::vector<int> run_list              = {};
    sumjs["good_run_list"]                 = run_list;
  }

  // loop over all runs and add info to the summary where relevant
  for (json::iterator it = rundb.begin(); it != rundb.end(); ++it) {
    int  run   = std::stoi(it.key());
    auto runjs = it.value();

    const std::string target  = runjs["target"]["target_label"];
    const double      p_hms   = runjs["spectrometers"]["hms_momentum"].get<double>();
    const double      th_hms  = runjs["spectrometers"]["hms_angle"].get<double>();
    const double      p_shms  = runjs["spectrometers"]["shms_momentum"].get<double>();
    const double      th_shms = runjs["spectrometers"]["shms_angle"].get<double>();
    const int         ps6     = runjs["daq"]["ps6"].get<int>();
    if (ps6 < 0) {
      continue;
    }
    // Calculate the prescale factor
    const double ps_factor = (ps6 == 0) ? 1. : (std::pow(2, ps6 - 1) + 1);
    if (vcsdb.find(it.key()) == vcsdb.end()) {
      continue;
    }
    if (countdb_shms[it.key()].find("good_total_charge") == countdb_shms[it.key()].end()) {
      continue;
    }
    // skip "empty " runs
    const double charge = countdb_shms[it.key()]["good_total_charge"].get<double>();
    if (charge < .1) {
      continue;
    }
    // skip junk runs
    if (commentdb.count(it.key()) > 0 &&
        commentdb[it.key()].find("comment") != commentdb[it.key()].end()) {
      std::string comment = commentdb[it.key()]["comment"];
      if (comment.find("junk") != std::string::npos || comment.find("Junk") != std::string::npos ||
          comment.find("JUNK") != std::string::npos) {
        continue;
      }
    }
    // loop over all settings
    for (json::iterator iset = jsum.begin(); iset != jsum.end(); ++iset) {
      auto& sumjs = iset.value();

      if (sumjs["target"] == target && fabs(sumjs["hms_angle"].get<double>() - th_hms) < .1 &&
          fabs(sumjs["hms_momentum"].get<double>() - p_hms) < .015 &&
          fabs(sumjs["shms_angle"].get<double>() - th_shms) < .1 &&
          fabs(sumjs["shms_momentum"].get<double>() - p_shms) < .015) {
        const double n_pi0 = vcsdb[it.key()]["missing_mass"]["peak1"]["sum_withbkg"].get<double>();
        sumjs["total_charge"] = sumjs["total_charge"].get<double>() + charge;
        sumjs["pi0_count"]    = sumjs["pi0_count"].get<double>() + n_pi0 * ps_factor;
        if (sumjs.find("mc_estimates") != sumjs.end()) {
          const double vcs_over_pi0 = sumjs["mc_estimates"]["vcs_over_pi0"];
          const double eff_cut1     = sumjs["mc_estimates"]["eff_cut1"];
          const double eff_cut2     = sumjs["mc_estimates"]["eff_cut2"];
          const double eff_cut3     = sumjs["mc_estimates"]["eff_cut3"];

          sumjs["projections"]["vcs_count"] =
              sumjs["projections"]["vcs_count"].get<double>() + n_pi0 * vcs_over_pi0 * ps_factor;
          sumjs["projections"]["vcs_count_cut1"] =
              sumjs["projections"]["vcs_count_cut1"].get<double>() +
              n_pi0 * vcs_over_pi0 * eff_cut1 * ps_factor;
          sumjs["projections"]["vcs_count_cut2"] =
              sumjs["projections"]["vcs_count_cut2"].get<double>() +
              n_pi0 * vcs_over_pi0 * eff_cut2 * ps_factor;
          sumjs["projections"]["vcs_count_cut3"] =
              sumjs["projections"]["vcs_count_cut3"].get<double>() +
              n_pi0 * vcs_over_pi0 * eff_cut3 * ps_factor;
        }

        auto run_list = sumjs["good_run_list"].get<std::vector<int>>();
        run_list.push_back(run);
        sumjs["good_run_list"] = run_list;

      } else {
        ;  // do nothing
      }
    }
  }

  {
    std::ofstream json_output_file{"database/status.json"};
    json_output_file << std::setw(4) << jsum << "\n";
  }
  std::cout << "Updated status written to Online/database/status.json\n";
}
