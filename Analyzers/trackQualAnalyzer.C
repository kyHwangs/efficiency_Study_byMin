#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TSystem.h>

#include "MuonHLTNtupleRun3.h"

using namespace std;


void DEBUG(int i, TString str = "")
{
    cout << "DEBUG: " << i << "\t" << str << endl;
}

void printRunTime(TStopwatch timer_)
{
    Double_t cpuTime = timer_.CpuTime();
    Double_t realTime = timer_.RealTime();

    cout << endl;
    cout << "************************************************" << endl;
    cout << "Total real time: " << realTime << " (seconds)" << endl;
    cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
    cout << "  CPU time / real time = " << cpuTime / realTime << endl;
    cout << "************************************************" << endl;
}

void printRunTimeShort(TStopwatch timer_)
{
    // Double_t cpuTime = timer_.CpuTime();
    Double_t realTime = timer_.RealTime();

    cout << "Total real time: " << realTime << " (seconds)" << endl;
    // cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
}

static inline void printMemory( TString tab = "" )
{
    ifstream proc_status("/proc/self/status");
    string buffer;
    while (proc_status.peek() != EOF) {
        getline(proc_status, buffer);
        TString str = buffer;
        if(str.Contains("RSS")) {
            cout << tab << str << endl;
            break;
        }
    }
}

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if( x == n )
        cout << endl;

    if ( x % (n/r +1) != 0 ) return;

    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;

    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );

    // Show the load bar.
    for (int x=0; x<c; x++) cout << "=";

    for (int x=c; x<w; x++) cout << " ";

    // ANSI Control codes to go back to the
    // previous line and clear it.
    cout << "]\r" << flush;
}


const static int n_eta_bins = 15;
double eta_bins[n_eta_bins] = {
  -2.4, -2.1, -1.6, -1.2, -0.9,
  -0.3, -0.2,  0.0,  0.2,  0.3,
   0.9,  1.2,  1.6,  2.1,  2.4
};

const static int n_pt_bins = 16;
double pt_bins[n_pt_bins] = {
     0, 10, 12, 15, 19,
    25, 31, 39, 50, 63,
    79, 100, 140, 200, 500,
    1000
};

struct sort_by_pt
{
    inline bool operator() (const Object& a, const Object& b)
    {
        return (a.pt > b.pt);
    }
};

const double MU_MASS = 0.1056583745;

bool acceptance(Object obj)
{
    return ( fabs(obj.eta) < 2.4 );
}

bool offlineSel(Object obj)
{
    bool out = (
        acceptance(obj) &&
        obj.get("isTight") &&
        obj.get("relPFIso") < 0.15
        //obj.get("isHighPtNew") &&
        //obj.get("relTrkIso") < 0.10
    );

    return out;
}

double invMass(Object obj1, Object obj2)
{
    TLorentzVector mu0, mu1;
    mu0.SetPtEtaPhiM(obj1.pt,
                     obj1.eta,
                     obj1.phi,
                     MU_MASS);
    mu1.SetPtEtaPhiM(obj2.pt,
                     obj2.eta,
                     obj2.phi,
                     MU_MASS);
    return (mu0+mu1).M();
}

void trackQualAnalyzer(
    TString ver = "v00", TString tag = "TEST",
    vector<TString> vec_Dataset = {}, TString JobId = "",
    TString outputDir = "./",
    const bool doDimuon = false, double ZmassWindow = -1,
    Int_t maxEv = -1, bool doMem = false, int nMem = 10001, bool doBar = true  // HERE
) {
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);

    TStopwatch timer_total;
    timer_total.Start();

    // -- Input
    vector<TString> paths = vec_Dataset;
    if(tag == "TEST") {
        paths = { "./ntuple_1.root" };
    }

    // -- Output
    TString fileName = TString::Format( "hist-%s-%s", ver.Data(), tag.Data() );
    if(JobId != "")  fileName = fileName + TString::Format("--%s", JobId.Data());
    TFile *f_output = TFile::Open(outputDir+fileName+"-trkQual.root", "RECREATE");

    // -- Event chain
    TChain *_chain_Ev          = new TChain("ntupler/ntuple");
    for(size_t f = 0; f<paths.size(); ++f) {
        _chain_Ev->Add(paths[f]);
        cout << "Adding path: " << paths[f] << endl;
    }
    cout << endl;

    unsigned nEvent      = _chain_Ev->GetEntries();
    if(maxEv >= 0)  nEvent = maxEv;
    cout << "\t nEvent: " << nEvent << endl;

    vector<TString> trk_vars = {
        "inner_trkChi2",
        "inner_validFraction",
        "inner_trackerLayers",
        "inner_trackerHits",
        "inner_lostTrackerHits",
        "inner_lostTrackerHitsIn",
        "inner_lostTrackerHitsOut",
        "inner_lostPixelHits",
        "inner_lostPixelBarrelHits",
        "inner_lostPixelEndcapHits",
        "inner_lostStripHits",
        "inner_lostStripTIBHits",
        "inner_lostStripTIDHits",
        "inner_lostStripTOBHits",
        "inner_lostStripTECHits",
        "inner_pixelLayers",
        "inner_pixelHits",

        "global_muonHits",
        "global_trkChi2",
        "global_trackerLayers",
        "global_trackerHits",
    };
    vector<TString> trkonly_vars = {
        "muonHits",
        "trkChi2",
        "validFraction",
        "trackerLayers",
        "trackerHits",
        "lostTrackerHits",
        "lostTrackerHitsIn",
        "lostTrackerHitsOut",
        "pixelLayers",
        "pixelHits",
    };
    vector<TString> muon_vars = {
        "isGLB",
        "isSTA",
        "isTRK",

        "inner_trkChi2",
        "inner_validFraction",
        "inner_trackerLayers",
        "inner_trackerHits",
        "inner_lostTrackerHits",
        "inner_lostTrackerHitsIn",
        "inner_lostTrackerHitsOut",
        "inner_lostPixelHits",
        "inner_lostPixelBarrelHits",
        "inner_lostPixelEndcapHits",
        "inner_lostStripHits",
        "inner_lostStripTIBHits",
        "inner_lostStripTIDHits",
        "inner_lostStripTOBHits",
        "inner_lostStripTECHits",
        "inner_pixelLayers",
        "inner_pixelHits",

        "global_muonHits",
        "global_trkChi2",
        "global_trackerLayers",
        "global_trackerHits",

        "momentumChi2",
        "positionChi2",
        "glbKink",
        "glbTrackProbability",
        "globalDeltaEtaPhi",
        "localDistance",
        "staRelChi2",
        "tightMatch",
        "trkKink",
        "trkRelChi2",
        "segmentCompatibility",
    };

    unique_ptr<MuonHLTNtupleRun3>  nt( new MuonHLTNtupleRun3( _chain_Ev, {""} ) );

    // -- Histograms
    TH1D *h_nEvents = new TH1D("h_nEvents",  "", 3, -1, 2);
    TH1D *h_nRuns = new TH1D("h_nRuns",  "", 100, 350000, 400000);

    vector<TString> L3types = {
        "iterL3OI",
        "iterL3IOFromL2",
        "iterL3FromL2",
        "iterL3IOFromL1", //Track-only
        "muon",
        "iterL3MuonNoID",
        "iterL3Muon",
    };

    // -- To compare
    vector<TString> etas = {
      "in",  // 1.5 < eta < 1.7
      "out",
    };
    vector<TString> runs = {
      "before", // before Run 368220 (<= 367910)
      "after", // after Run 368220 (>= 368223)
    };

    vector<vector<TH1D *>> vh_L3types_in_before  = {};
    vector<vector<TH1D *>> vh_L3types_in_after  = {};
    vector<vector<TH1D *>> vh_L3types_out_before  = {};
    vector<vector<TH1D *>> vh_L3types_out_after  = {};

    vector<vector<TH2D *>> vh_L3types_before  = {};
    vector<vector<TH2D *>> vh_L3types_after  = {};

    for (unsigned i=0; i<L3types.size(); ++i) {
        vh_L3types_in_before.push_back( {} );
        vh_L3types_in_after.push_back( {} );
        vh_L3types_out_before.push_back( {} );
        vh_L3types_out_after.push_back( {} );

        vh_L3types_before.push_back( {} );
        vh_L3types_after.push_back( {} );

        vector<TString> vars;
        if (i==3) vars = trkonly_vars;
        else if(i<3) vars = trk_vars;
        else vars = muon_vars;

        for (unsigned j=0; j<vars.size(); ++j) {
            TString var = vars.at(j);
            int nbins = 101;
            double xmax = 1.01;
            if(var.Contains("Kink")) xmax = 101;
            else if(var.Contains("lost")) {nbins = 10; xmax = 10;}
            else if(var.Contains("Hits")) {nbins = 60; xmax = 60;}
            else if(var.Contains("Layers")) {nbins = 20; xmax = 20;}
            else if(var.Contains("is") || var.Contains("tightMatch")) {nbins = 2; xmax = 2;}
            else if(var.Contains("Probability") || var.Contains("position") || var.Contains("momentum")) xmax = 30.3;
            else if(var.Contains("Chi2") || var.Contains("localDistance")) xmax = 10.1;
            else if(var.Contains("globalDelta")) xmax = 0.101;

            TString name = TString::Format("h_%s_%s_in_before", L3types.at(i).Data(), var.Data());
            TH1D *h_L3type_in_before  = new TH1D(name,  "", 101, 0, xmax);
            vh_L3types_in_before.at(i).push_back( h_L3type_in_before );

            name = TString::Format("h_%s_%s_in_after", L3types.at(i).Data(), var.Data());
            TH1D *h_L3type_in_after  = new TH1D(name,  "", 101, 0, xmax);
            vh_L3types_in_after.at(i).push_back( h_L3type_in_after );

            name = TString::Format("h_%s_%s_out_before", L3types.at(i).Data(), var.Data());
            TH1D *h_L3type_out_before  = new TH1D(name,  "", 101, 0, xmax);
            vh_L3types_out_before.at(i).push_back( h_L3type_out_before );

            name = TString::Format("h_%s_%s_out_after", L3types.at(i).Data(), var.Data());
            TH1D *h_L3type_out_after  = new TH1D(name,  "", 101, 0, xmax);
            vh_L3types_out_after.at(i).push_back( h_L3type_out_after );

            name = TString::Format("h_%s_%s_2D_before", L3types.at(i).Data(), var.Data());
            TH2D *h_L3type_before  = new TH2D(name,  "", 48, -2.4, 2.4, 30, -TMath::Pi(), TMath::Pi());
            vh_L3types_before.at(i).push_back( h_L3type_before );

            name = TString::Format("h_%s_%s_2D_after", L3types.at(i).Data(), var.Data());
            TH2D *h_L3type_after  = new TH2D(name,  "", 48, -2.4, 2.4, 30, -TMath::Pi(), TMath::Pi());
            vh_L3types_after.at(i).push_back( h_L3type_after );
        }
    }

    // -- Event loop
    for(unsigned i_ev=0; i_ev<nEvent; i_ev++) {
        if(doBar)
            loadBar(i_ev+1, nEvent, 100, 100);
        else if( doMem && i_ev !=0 && i_ev % nMem == 0 )
            printMemory("\t");
        else
            printRunTimeShort(timer_total);

        nt->GetEntry( i_ev );

        double genWeight = nt->genEventWeight > 0.0 ? 1.0 : -1.0;
        if (nt->isRealData)
            genWeight = 1.;
        h_nEvents->Fill( genWeight );
        h_nRuns->Fill( nt->runNum, genWeight );

        // -- Get object collections
        vector<Object> iterL3OI = nt->get_iterL3OI();
        vector<Object> iterL3IOFromL2 = nt->get_iterL3IOFromL2();
        vector<Object> iterL3FromL2 = nt->get_iterL3FromL2();
        vector<Object> iterL3IOFromL1 = nt->get_iterL3IOFromL1();

        vector<Object> muons = nt->get_offlineMuons();
        vector<Object> iterL3MuonNoID = nt->get_iterL3MuonNoID();
        vector<Object> iterL3Muon = nt->get_iterL3Muon();

        vector<vector<Object>*> L3MuonColls {
            &iterL3OI,
            &iterL3IOFromL2,
            &iterL3FromL2,
            &iterL3IOFromL1,
            &muons,
            &iterL3MuonNoID,
            &iterL3Muon,
        };

        // -- run                                                                                                                                                                                                       
        for (unsigned irun = 0; irun < runs.size(); ++irun) {
          if (irun == 0 && nt->runNum > 368220) continue;
          if (irun == 1 && nt->runNum < 368220) continue;

          //### L3types loop ###
          for (unsigned i=0; i<L3types.size(); ++i) {
            vector<Object>* L3Coll = L3MuonColls.at(i);
            TString L3type = L3types.at(i);

            vector<TString> vars;
            if (i==3) vars = trkonly_vars;
            else if(i<3) vars = trk_vars;
            else vars = muon_vars;

            // -- L3 objects loop
            for (auto& L3: *L3Coll) {
              for (unsigned ieta = 0; ieta < etas.size(); ++ieta) {
                if (ieta == 0 && (L3.eta < 1.5 || L3.eta > 1.7)) continue; 
                if (ieta == 1 && (L3.eta > 1.5 && L3.eta < 1.7)) continue;

                for (unsigned j=0; j<vars.size(); ++j) {
                  if (irun == 0 && ieta == 0) vh_L3types_in_before.at(i).at(j)->Fill(L3.get(vars.at(j)), genWeight);
                  else if (irun == 0 && ieta == 1) vh_L3types_out_before.at(i).at(j)->Fill(L3.get(vars.at(j)), genWeight);
                  else if (irun == 1 && ieta == 0) vh_L3types_in_after.at(i).at(j)->Fill(L3.get(vars.at(j)), genWeight);
                  else if (irun == 1 && ieta == 1) vh_L3types_out_after.at(i).at(j)->Fill(L3.get(vars.at(j)), genWeight);

                  if (vars.at(j).Contains("lost")) {
                    if (irun == 0) vh_L3types_before.at(i).at(j)->Fill(L3.eta, L3.phi, L3.get(vars.at(j)));
                    else vh_L3types_after.at(i).at(j)->Fill(L3.eta, L3.phi, L3.get(vars.at(j)));
                  }
                }
              }
            }
          }
        }
    }

    // -- Save output and Clear memory
    // delete _chain_Ev;

    f_output->cd();

    h_nEvents->Write();
    h_nRuns->Write();

    for(unsigned i=0; i<L3types.size(); ++i) {
        TDirectory* dir = f_output->mkdir(L3types.at(i));
        dir->cd();

        vector<TString> vars;
        if (i==3) vars = trkonly_vars;
        else if(i<3) vars = trk_vars;
        else vars =muon_vars;

        for (unsigned j=0; j<vars.size(); ++j) {
            vh_L3types_in_before.at(i).at(j)->Write();
            vh_L3types_out_before.at(i).at(j)->Write();
            vh_L3types_in_after.at(i).at(j)->Write();
            vh_L3types_out_after.at(i).at(j)->Write();

            if (vars.at(j).Contains("lost")) {
                vh_L3types_before.at(i).at(j)->Write();
                vh_L3types_after.at(i).at(j)->Write();
            }
        }
        f_output->cd();
    }
    f_output->Close();
    printRunTime(timer_total);
}


