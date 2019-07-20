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
//    (2) Calorimeter cut should be fine,
//        make sure the Cherenkov cut is strict enough
// =================================================================================
std::string goodTrack = "H.gtr.dp > -8 && H.gtr.dp < 8 &&"
                        "-0.08 < H.gtr.th && H.gtr.th < 0.06 && "
                        "-0.03 < H.gtr.ph && H.gtr.ph < 0.04"
                        "&& TMath::Abs(H.gtr.y) < 2.0";
//std::string electronCut = "H.cal.etottracknorm > 0.80 && H.cal.etottracknorm < 2.&&"
std::string electronCut = "H.hod.betanotrack>0.55 && H.hod.betanotrack<0.8";

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

int good_hms_counter(int RunNumber = 8585, int nevents = -1, const std::string& mode = "coin",
                     int update = 1) {
  // Sanity check the input
  if (mode != "coin" && mode != "hms") {
    std::cerr << "error: please specify a valid mode, either 'coin' or 'hms'" << std::endl;
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
  // was this data taken with ps3 or ps4 or coin?
  int  ps3             = -1;    // HMS ELREAL
  int  ps4             = -1;    // HMS 3/4
  int  ps5             = -1;    // COIN: SHMS 3/4 & HMS ELREAL
  int  ps6             = -1;    // COIN: SHMS 3/4 & HMS 3/4
  bool singles_trigger = true;  // Was there a HMS singles trigger for this run?
  if (db[runnum_str].find("daq") != db[runnum_str].end()) {
    ps3 = db[runnum_str]["daq"]["ps3"].get<int>();
    ps4 = db[runnum_str]["daq"]["ps4"].get<int>();
    ps5 = db[runnum_str]["daq"]["ps5"].get<int>();
    ps6 = db[runnum_str]["daq"]["ps6"].get<int>();
    std::cout << "ps3 = " << ps3 << " and ps4 = " << ps4 << "\n";
    std::cout << "ps5 = " << ps5 << " and ps6 = " << ps6 << "\n";
  } else {
    std::cout << "Error: pre-scaler unspecified in " << run_list_fname << std::endl;
    return -127;
  }

  // first try the singles trigger
  int ps = std::max(ps3, ps4);
  if (ps3 > 0 && ps4 > 0) {
    std::cout << "Inconsistent values for ps3 and ps4, only can be enabled at the same time\n";
    std::cout << "(ps3 = " << ps3 << " and ps4 = " << ps4 << ")" << std::endl;
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
    std::cout << "Selected " << ((ps3 > ps4) ? "ps3" : "ps4") << std::endl;
  }
  // Should never happen
  if (ps == -1) {
    std::cout << "ERROR: no pre-scaler was set for the HMS, unable to proceed." << std::endl;
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

  // Select HMS singles only
  // event type 2 is a singles event (ps3 or ps4)
  // event type 4 is a coincidence event (ps5 or ps6)
  auto dHMS = d.Filter(singles_trigger ? "fEvtHdr.fEvtType == 2" : "fEvtHdr.fEvtType == 4");

  // Good track cuts
  auto dGoodTrack = dHMS.Filter(goodTrack);

  // PID cuts
  auto dEl = dGoodTrack.Filter(electronCut);

  // Data frame index to make histogram definition a bit less cumbersome
  // here we have three versions of the data corresponding to the 3 stages of the
  // dataframe filters:
  //  1. "raw" --> after event type selection
  //  2. "Cuts: tracking" --> after the goodTrack cuts
  //  3. "Cuts: tracking+PID" --> after goodTrack && electronCut
  std::vector<std::pair<std::string, RDFInfo>> dfs = {{"raw", {dHMS, "HMS"}},
                                                      {"tracked", {dGoodTrack, "Cuts: tracking"}},
                                                      {"identified", {dEl, "Cuts: tracking+PID"}}};

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
    h1D["HcalEP"][name] =
        df_info.df.Histo1D({("H.cal.etottracknorm_" + name).c_str(),
                            (df_info.title + ";HMS E/P;counts").c_str(), 200, -.5, 2.},
                           "H.cal.etottracknorm");
    // Cherenkov
    h1D["HcerNphe"][name] =
        df_info.df.Histo1D({("H.cer.npeSum_" + name).c_str(),
                            (df_info.title + "; HMS Cer #phe; counts ").c_str(), 200, -1, 40},
                           "H.cer.npeSum");
    // Delta
    h1D["Hdp"][name] =
        df_info.df.Histo1D({("H.gtr.dp_" + name).c_str(),
                            (df_info.title + ";#deltap [%];counts").c_str(), 200, -30, 40},
                           "H.gtr.dp");
    // theta
    h1D["Hth"][name] =
        df_info.df.Histo1D({("H.gtr.th_" + name).c_str(),
                            (df_info.title + ";#theta_{HMS};counts ").c_str(), 200, -0.1, 0.1},
                           "H.gtr.th");
    // phi
    h1D["Hph"][name] =
        df_info.df.Histo1D({("H.gtr.ph_" + name).c_str(),
                            (df_info.title + ";#phi_{HMS};counts ").c_str(), 200, -0.1, 0.1},
                           "H.gtr.ph");
    // y
    h1D["Hy"][name] =
        df_info.df.Histo1D({("H.gtr.y_" + name).c_str(),
                            (df_info.title + ";ytar_{HMS};counts ").c_str(), 200, -10, 10},
                           "H.gtr.y");
  }

  // =========================================================================================
  // scalers and counts
  // =========================================================================================
  auto total_charge  = d_sh.Max("P.BCM4B.scalerChargeCut");
  auto time_1MHz_cut = d_sh.Max("P.1MHz.scalerTimeCut");

  auto count_raw        = dHMS.Count();
  auto count_tracked    = dGoodTrack.Count();
  auto count_identified = dEl.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------
  const double good_total_charge = *total_charge / 1000.0;  // mC
  const double good_time         = *time_1MHz_cut;          // s

  std::map<std::string, double> counts = {
      {"count_e", (*count_identified)},
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
    std::ifstream input_file(db_path + "/countdb_hms.json");
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
    std::cout << "Updating " << db_path << "/countdb_hms.json with hms counts\n";
    std::ofstream output_file(db_path + "/countdb_hms.json");
    output_file << std::setw(4) << countdb << "\n";
  }
  // =====================================================================================
  // Create monitoring plots
  // =====================================================================================
  const std::string plot_name = fmt::format("{}/mon_hms_{}.pdf", mon_path, RunNumber);

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
