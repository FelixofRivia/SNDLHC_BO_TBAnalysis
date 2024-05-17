#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "TChain.h"
#include "TClonesArray.h"
#include "SNDLHCEventHeaderConst.h"
#include "./libs/Inclusion.h"
#include "./libs/SciFiPlaneView.h"
#include "./libs/USPlaneView.h"

cfg setCfg( bool istb, int runN ) {
  cfg config;
  if (istb) {
    config.SCIFI_STATIONS = 4;
    config.SCIFI_BOARDPERPLANE = 1;
    config.SCIFI_NCHANNELS = 512;
    config.SCIFI_TIMECUT = 1;  
    config.SCIFI_DIMCLUSTER = 35;
    config.SCIFI_GAPCLUSTER = 2;
    config.SCIFI_DENSITYWINDOW = 128;
    config.SCIFI_DENSITYHITS = 36;
    config.SCIFI_F = 0.17; 

    config.US_STATIONS = 5;
    config.US_TIMECUT = 3;

    config.SCIFI_DIM = 13;

    config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/";
    config.OUTFILENAME = "output/TB_output";
  }
  else {
    config.SCIFI_STATIONS = 5;
    config.SCIFI_BOARDPERPLANE = 3;
    config.SCIFI_NCHANNELS = 512*3;
    config.SCIFI_TIMECUT = 1;  
    config.SCIFI_DIMCLUSTER = 18;
    config.SCIFI_GAPCLUSTER = 2;
    config.SCIFI_DENSITYWINDOW = 128;
    config.SCIFI_DENSITYHITS = 18;
    config.SCIFI_F = 0.17;

    config.US_STATIONS = 5;
    config.US_TIMECUT = 3;

    config.SCIFI_DIM = 13*3;
    if (runN < 5422) {
      config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/physics/2022/";
    }
    else {
      config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/physics/2023/";
    }
    config.OUTFILENAME = "output/TI18_output";
  }
  return config;
}

//*****************************************************
// Skim function:
// THE ONLY PART OF CODE ONE NEEDS TO CHANGE
//
// Here we assume:
// - that the cut is base donly on SciFi hits
// - that they come sorted by station in the TClonesArray
//*****************************************************
std::vector<SciFiPlaneView> fillSciFi(cfg configuration, TClonesArray *sf_hits){

  std::vector<SciFiPlaneView> scifi_planes;

  int begin{0};
  int count{0};

  int n_sf_hits{sf_hits->GetEntries()};

  for (int st{1}; st <= configuration.SCIFI_STATIONS; ++st) {
    begin = count;
    while (count < n_sf_hits &&
           st == static_cast<sndScifiHit *>(sf_hits->At(count))->GetStation()) {
      ++count;
    }

    auto plane = SciFiPlaneView(configuration, sf_hits, begin, count, st);
    scifi_planes.emplace_back(plane);
  }

  return scifi_planes;
}

bool vetoCut(TClonesArray *mufi_hits){

  int n_mufi_hits{mufi_hits->GetEntries()};
  int system{0};
  // skip veto/beam monitor
  for (int i{0}; i<n_mufi_hits && system<2; ++i){
    system = static_cast<MuFilterHit *>(mufi_hits->At(i))->GetSystem();
    if (system == 1) return false;
  }
  return true;
}

bool isOneVetoHit(TClonesArray *mufi_hits){

  int n_mufi_hits{mufi_hits->GetEntries()};
  int system{0};
  int n_hits{0};
  // skip veto/beam monitor
  for (int i{0}; i<n_mufi_hits && system<2; ++i){
    system = static_cast<MuFilterHit *>(mufi_hits->At(i))->GetSystem();
    if (system == 1) {
      n_hits++;
      if (n_hits>1) return false;
    }
  }
  return n_hits==1 ? true : false;
}

bool isStableBeam(SNDLHCEventHeader *header){

  return (header->GetBeamMode() == static_cast<int>(LhcBeamMode::StableBeams)) ? true : false;
}

void timeCutGuil (std::vector<SciFiPlaneView> &Scifi) {
  TH1D* times = new TH1D ("times", "times; clk cycles; entries", 1000, 0, 50);

  for (auto &plane : Scifi) {
    auto time = plane.getTime();
    for (int i{0}; i<time.x.size(); ++i) {
      if (time.x[i]>DEFAULT) {
        times->Fill(time.x[i]);
      }
      if (time.y[i]>DEFAULT) {
        times->Fill(time.y[i]);
      }
    }
  }

  double referenceTime = times->GetXaxis()->GetBinCenter(times->GetMaximumBin()); // in clk cycles

  for (auto &plane : Scifi) {
    plane.timeCut(referenceTime - 0.5, referenceTime + 0.5);
  }

  delete times;
}

int checkShower_with_density(std::vector<SciFiPlaneView> scifi_planes ) {
  //find start of shower
  for (auto &plane : scifi_planes) {
    if (plane.infoDensity(plane.getConfig().SCIFI_DENSITYWINDOW, plane.getConfig().SCIFI_DENSITYHITS)) return plane.getStation();
  }
  return -1;
}

bool skim_function(cfg configuration, SNDLHCEventHeader *header, TClonesArray *sf_hits, TClonesArray *mufi_hits) {
  // Stable beam
  if (!isStableBeam(header)) return false;
  // No hit in Veto
  if (!vetoCut(mufi_hits)) return false;
  // One hit in Veto
  // if (!isOneVetoHit(mufi_hits)) return false;
  // Shower tagged
  auto scifi_planes = fillSciFi(configuration, sf_hits);
  timeCutGuil(scifi_planes);
  if (checkShower_with_density(scifi_planes) == -1) return false;

  return true;
};

//*****************************************************
// Steering function called by run_skim.sh
//*****************************************************

void skim(std::string file_name, int run_number, std::string out_folder, bool isTB) {

  auto start = std::chrono::steady_clock::now();

  std::cout << "[skim] processing: " << file_name << std::endl;

  // ##################### Set right parameters for data type (TB/TI18)
  cfg configuration = setCfg(isTB, run_number); 

  // ##################### Read input file #####################
  auto tree = new TChain("rawConv");
  tree->Add(file_name.c_str());

  auto header = new SNDLHCEventHeader();
  auto sf_hits = new TClonesArray("sndScifiHit");
  auto mu_hits = new TClonesArray("MuFilterHit");

  if (!isTB && run_number<5422){
    tree->SetBranchAddress("EventHeader", &header);
  }
  else {
    tree->SetBranchAddress("EventHeader.", &header);
  }
  tree->SetBranchAddress("Digi_ScifiHits", &sf_hits);
  tree->SetBranchAddress("Digi_MuFilterHits", &mu_hits);

  // ##################### Create output file ##################
  TFile output_file((std::filesystem::path(out_folder) /
                     std::filesystem::path(file_name).filename())
                        .c_str(),
                    "RECREATE");
  output_file.cd();

  auto new_tree = tree->CloneTree(0);

  // ##################### Event loop ##########################
  auto n_entries{tree->GetEntries()};

  for (long long i{}; i != n_entries; ++i) {
    if (i % 100000 == 0)
      std::cout << "[skim] reading event: " << i << std::endl;

    tree->GetEntry(i);

    if (skim_function(configuration, header, sf_hits, mu_hits)) {
      new_tree->Fill();
    }
  }

  output_file.Write();

  auto finish = std::chrono::steady_clock::now();

  using namespace std::chrono_literals;
  std::cout << "[skim] efficiency: " << std::setprecision(9)
            << (100.0 * new_tree->GetEntries() / tree->GetEntries()) << " %"
            << std::endl;
  std::cout << "[skim] elapsed time: " << (finish - start) / 60s << " minutes"
            << std::endl;

  output_file.Close();
};
