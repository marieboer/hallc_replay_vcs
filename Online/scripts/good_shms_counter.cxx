#include <cmath>
#include <iostream>
#include "nlohmann/json.hpp"

#include "ROOT/RDataFrame.hxx"

#include "fmt/format.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

const std::string db_path  = "database";
const std::string mon_path = "monitoring";

// =================================================================================
// Some syntactic sugar to make bulk-definition of histograms easier
// =================================================================================
using RDFNode = decltype(ROOT::RDataFrame{0}.Filter(""));
using Histo1DProxy =
    decltype(ROOT::RDataFrame{0}.Histo1D(ROOT::RDF::TH1DModel{"", "", 128u, 0., 0.}, ""));

struct RDFInfo {
  RDFNode&          df;
  const std::string title;
  RDFInfo(RDFNode& df, std::string_view title) : df{df}, title{title} {}
};

// =================================================================================
// Cuts - online:
//    (1) make sure that the phi/theta/y cuts are sensible
//    (2) Make sure hadron PID cuts are as desided
//        make sure the Cherenkov cut is strict enough
// =================================================================================
std::string goodTrack = "P.gtr.dp > -10 && P.gtr.dp < 22 &&"
                        "-0.05 < P.gtr.th && P.gtr.th < 0.05 && "
                        "-0.035 < P.gtr.ph && P.gtr.ph < 0.025"
                        "&& P.gtr.y > -2.0 && P.gtr.y < 3.0";
//std::string hadronCut = "P.cal.etottracknorm < 0.80 &&"
std::string hadronCut = "P.cal.etottracknorm > 0.80 ";

// =================================================================================
// Utility function to return the first good root file for this run
// =================================================================================
std::string get_root_file(std::vector<std::string> file_names) {
  for (const auto& name : file_names) {
    if (!gSystem->AccessPathName(name.c_str())) {
      TFile file(name.c_str());
      if (file.IsZombie()) {
        std::cout << name << " is a zombie.\n";
        std::cout
            << " Did your replay finish?  Check that the it is done before running this script.\n";
        return "";
        // return;
      } else {
        std::cout << " using : " << name << "\n";
        return name;
      }
    }
  }
  return "";
}

int good_shms_counter(int RunNumber = 8585, int nevents = -1, const std::string& mode = "coin",
                      int update = 1) {
  // Sanity check the input
  if (mode != "coin" && mode != "shms") {
    std::cerr << "error: please specify a valid mode, either 'coin' or 'shms'" << std::endl;
    return -127;
  }

  // ==========================================================================================
  // Load the run info from the database
  // ==========================================================================================
  using nlohmann::json;
  const std::string run_list_fname = db_path + "/rundb_" + mode + ".json";
  std::cout << "Loading run info from: " << run_list_fname << std::endl;
  nlohmann::json db;
  {
    std::ifstream json_input_file(run_list_fname);
    try {
      json_input_file >> db;
    } catch (json::parse_error) {
      std::cerr << "error: json file, " << run_list_fname
                << ", is incomplete or has broken syntax.\n";
      return -127;
    }
  }

  auto runnum_str = std::to_string(RunNumber);
  if (db.find(runnum_str) == db.end()) {
    std::cout << "Run " << RunNumber << " not found in " << run_list_fname << "\n";
    std::cout << "Check that run number and replay exists. \n";
    std::cout << "If the problem persists please contact Sylvester (217-848-0565)\n";
    return -127;
  }

  // ==========================================================================================
  // Get the pre-scaler info from the database
  // ==========================================================================================
  // was this data taken with ps1 or ps2 or coin?
  int  ps1             = -1;    // SHMS 3/4
  int  ps2             = -1;    // SHMS ELREAL
  int  ps5             = -1;    // COIN: SHMS 3/4 & HMS ELREAL
  int  ps6             = -1;    // COIN: SHMS 3/4 & HMS 3/4
  bool singles_trigger = true;  // Was there a SHMS singles trigger for this run?
  if (db[runnum_str].find("daq") != db[runnum_str].end()) {
    ps1 = db[runnum_str]["daq"]["ps1"].get<int>();
    ps2 = db[runnum_str]["daq"]["ps2"].get<int>();
    ps5 = db[runnum_str]["daq"]["ps5"].get<int>();
    ps6 = db[runnum_str]["daq"]["ps6"].get<int>();
    std::cout << "ps1 = " << ps1 << " and ps2 = " << ps2 << "\n";
    std::cout << "ps5 = " << ps5 << " and ps6 = " << ps6 << "\n";
  } else {
    std::cout << "Error: pre-scaler unspecified in " << run_list_fname << std::endl;
    return -127;
  }

  // first try the singles trigger
  int ps = std::max(ps1, ps2);
  if (ps1 > 0 && ps2 > 0) {
    std::cout << "Inconsistent values for ps1 and ps2, only can be enabled at the same time\n";
    std::cout << "(ps1 = " << ps1 << " and ps2 = " << ps2 << ")" << std::endl;
    return -127;
    // if no singles trigger, try the coin trigger
  } else if (ps == -1) {
    std::cout
        << "Warning: no data with singles pre-scaler taking, using coincidence trigger instead\n";
    singles_trigger = false;
    ps              = std::max(ps5, ps6);
    if (ps5 > 0 && ps6 > 0) {
      std::cout << "Inconsistent values for ps5 and ps6, only can be enabled at the same time\n";
      std::cout << "(ps5 = " << ps5 << " and ps6 = " << ps6 << ")" << std::endl;
      return -127;
    }
    std::cout << "Selected " << ((ps5 > ps6) ? "ps5" : "ps6") << std::endl;
  } else {
    std::cout << "Selected " << ((ps1 > ps2) ? "ps1" : "ps2") << std::endl;
  }
  // Should never happen
  if (ps == -1) {
    std::cout << "ERROR: no pre-scaler was set for the SHMS, unable to proceed." << std::endl;
    return -127;
  }

  // Calculate the prescale factor
  const double ps_factor = (ps == 0) ? 1. : (std::pow(2, ps - 1) + 1);
  std::cout << "Using prescale factor " << ps_factor << std::endl;

  // ===============================================================================================
  // Get the ROOT file
  // ===============================================================================================
  std::vector<std::string> file_paths = {
      fmt::format("../full_online/{}_replay_production_{}_{}.root", mode, RunNumber, nevents),
      fmt::format("../ROOTfiles_volatile/{}_replay_production_{}_{}.root", mode, RunNumber,
                  nevents),
      fmt::format("../ROOTfiles/{}_replay_production_{}_{}.root", mode, RunNumber, nevents)};
  const std::string rootfile = get_root_file(file_paths);
  if (rootfile.size() == 0) {
    std::cout << " Error: suitable root file not found\n";
    return -127;
  }

  // ===============================================================================================
  // Dataframe
  // ===============================================================================================

  // run 12 threads in parallel
  ROOT::EnableImplicitMT(12);

  //---------------------------------------------------------------------------
  // Detector tree
  ROOT::RDataFrame d("T", rootfile);
  // SHMS Scaler tree
  ROOT::RDataFrame d_sh("TSP", rootfile);

  // Select SHMS singles only
  // event type 1 is a singles event (ps1 or ps2)
  // event type 4 is a coincidence event (ps5 or ps6)
  auto dSHMS = d.Filter(singles_trigger ? "fEvtHdr.fEvtType == 2" : "fEvtHdr.fEvtType == 4");

  // Good track cuts
  auto dGoodTrack = dSHMS.Filter(goodTrack);

  // PID cuts
  auto dHad = dGoodTrack.Filter(hadronCut);

  // Data frame index to make histogram definition a bit less cumbersome
  // here we have three versions of the data corresponding to the 3 stages of the
  // dataframe filters:
  //  1. "raw" --> after event type selection
  //  2. "Cuts: tracking" --> after the goodTrack cuts
  //  3. "Cuts: tracking+PID" --> after goodTrack && hadronCut
  std::vector<std::pair<std::string, RDFInfo>> dfs = {{"raw", {dSHMS, "SHMS"}},
                                                      {"tracked", {dGoodTrack, "Cuts: tracking"}},
                                                      {"identified", {dHad, "Cuts: tracking+PID"}}};

  // =========================================================================================
  // Histograms
  //  * store histograms in a "Map of maps":
  //     (plot name): (histogram name: histogram)
  //
  // =========================================================================================
  using Histo1DMap    = std::map<std::string, Histo1DProxy>;
  using Histo1DMapMap = std::map<std::string, Histo1DMap>;
  Histo1DMapMap h1D;
  for (auto& kval : dfs) {
    std::string name{kval.first};
    RDFInfo&    df_info{kval.second};
    // Calorimeter
    h1D["PcalEP"][name] =
        df_info.df.Histo1D({("P.cal.etottracknorm_" + name).c_str(),
                            (df_info.title + ";SHMS E/P;counts").c_str(), 200, -.5, 2.},
                           "P.cal.etottracknorm");
    // Cherenkov
    h1D["PcerNphe"][name] =
        df_info.df.Histo1D({("P.hgcer.npeSum_" + name).c_str(),
                            (df_info.title + "; SHMS HGCer #phe; counts ").c_str(), 200, -1, 40},
                           "P.hgcer.npeSum");
    // Delta
    h1D["Pdp"][name] =
        df_info.df.Histo1D({("P.gtr.dp_" + name).c_str(),
                            (df_info.title + ";#deltap [%];counts").c_str(), 200, -30, 40},
                           "P.gtr.dp");
    // theta
    h1D["Pth"][name] =
        df_info.df.Histo1D({("P.gtr.th_" + name).c_str(),
                            (df_info.title + ";#theta_{SHMS};counts ").c_str(), 200, -0.1, 0.1},
                           "P.gtr.th");
    // phi
    h1D["Pph"][name] =
        df_info.df.Histo1D({("P.gtr.ph_" + name).c_str(),
                            (df_info.title + ";#phi_{SHMS};counts ").c_str(), 200, -0.1, 0.1},
                           "P.gtr.ph");
    // y
    h1D["Py"][name] =
        df_info.df.Histo1D({("P.gtr.y_" + name).c_str(),
                            (df_info.title + ";ytar_{SHMS};counts ").c_str(), 200, -10, 10},
                           "P.gtr.y");
  }

  // =========================================================================================
  // scalers and counts
  // =========================================================================================
  auto total_charge  = d_sh.Max("P.BCM4B.scalerChargeCut");
  auto time_1MHz_cut = d_sh.Max("P.1MHz.scalerTimeCut");

  auto count_raw        = dSHMS.Count();
  auto count_tracked    = dGoodTrack.Count();
  auto count_identified = dHad.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------
  const double good_total_charge = *total_charge / 1000.0;  // mC
  const double good_time         = *time_1MHz_cut;          // s

  std::map<std::string, double> counts = {
      {"count_h", (*count_identified)},
      {"count_tracked", (*count_tracked)},
      {"count_raw", (*count_raw)},
      {"count_identified_pscorr", (*count_identified) * ps_factor},
      {"count_tracked_pscorr", (*count_tracked) * ps_factor},
      {"count_raw_pscorr", (*count_raw) * ps_factor},
      {"good_total_charge", good_total_charge},
      {"good_time", good_time}};

  // Update counts list
  json countdb;
  {
    std::ifstream input_file(db_path + "/countdb_shms.json");
    try {
      input_file >> countdb;
    } catch (json::parse_error) {
      std::cerr << "error: json file is incomplete or has broken syntax.\n";
      return -127;
    }
  }
  std::string run_str = std::to_string(RunNumber);
  std::cout << "----------------------------------------------------------" << std::endl;
  for (const auto& kv : counts) {
    if (kv.first.find("pscorr") != std::string::npos && ps_factor == 1) {
      continue;
    }
    std::cout << " " << kv.first << ": " << kv.second;
    if (kv.first.find("count") != std::string::npos) {
      if (kv.first.find("pscorr") == std::string::npos) {
        std::cout << " --> yield (counts / mC): " << kv.second / good_total_charge << " +- "
                  << sqrt(kv.second) / good_total_charge;
      } else {
        std::cout << " --> yield (counts / mC): " << kv.second / good_total_charge << " +- "
                  << (1 / sqrt(kv.second / ps_factor)) * kv.second / good_total_charge;
      }
    }
    std::cout << "\n";
    countdb[run_str][kv.first] = kv.second;
  }

  countdb[run_str]["charge bcm4b 2u cut"] = good_total_charge;
  countdb[run_str]["time 1MHz 2u cut"]    = good_time;
  countdb[run_str]["ps_factor"]           = ps_factor;

  if (update) {
    std::cout << "Updating " << db_path << "/countdb_shms.json with shms counts\n";
    std::ofstream output_file(db_path + "/countdb_shms.json");
    output_file << std::setw(4) << countdb << "\n";
  }
  // =====================================================================================
  // Create monitoring plots
  // =====================================================================================
  const std::string plot_name = fmt::format("{}/mon_shms_{}.pdf", mon_path, RunNumber);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TCanvas* c = new TCanvas("canvas", "", 800, 600);
  c->Print((plot_name + "[").c_str());
  for (auto& [name, plot] : h1D) {
    c->SetLogy();
    plot["raw"]->SetTitle(plot["raw"]->GetXaxis()->GetTitle());
    plot["raw"]->SetLineColor(kGreen + 2);
    plot["raw"]->SetLineColor(kGreen + 2);
    plot["raw"]->SetLineWidth(2);
    plot["raw"]->DrawClone();
    plot["tracked"]->SetLineColor(kMagenta + 2);
    plot["tracked"]->SetLineWidth(2);
    plot["tracked"]->DrawClone("same");
    plot["identified"]->SetLineColor(kBlue + 1);
    plot["identified"]->SetLineWidth(2);
    plot["identified"]->DrawClone("same");
    c->BuildLegend();
    c->Print(plot_name.c_str());
  }
  c->Print((plot_name + "]").c_str());

  return 0;
}
