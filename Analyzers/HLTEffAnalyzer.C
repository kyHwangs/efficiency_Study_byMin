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

class HistContainer
{
public:
    HistContainer(
      TString _Tag,
      vector<TString> _variables = { "pt", "eta", "phi", "nvtx"},//, "pu", "lumi"},
      vector<vector<double>> _ranges = {
        { 6000, 0, 3000 },
        { 48, -2.4, 2.4 },
        { 30, -TMath::Pi(), TMath::Pi() },
        { 75, 0, 75 }//,
        //{ 75, 0, 75 },
        //{ 25, 0, 2.5 }
      }
    ) {

      if(_variables.size() != _ranges.size()) {
        cout << "HistContainer: _variables.size() != _ranges.size()" << endl;
        exit(1);
      }

      this->Tag = _Tag;
      this->variables = _variables;
      this->ranges = _ranges;
      this->nVar = variables.size();

      this->Init();
    }

  void fill_den( Object obj, double Nvtx, double PU, double Lumi, double weight = 1.0 ) {
      for(int k=0; k < nVar; ++k) {
        if( variables[k] == "nvtx" ) {
          v_den[k]->Fill( Nvtx, weight );
        }
        else if( variables[k] == "pu" ) {
          v_den[k]->Fill( PU, weight );
        }
        else if( variables[k] == "lumi" ) {
          v_den[k]->Fill( Lumi, weight );
        }
        else if( variables[k] == "mva" ) {
          v_den[k]->Fill( 1./(1.+exp(-1.*obj.get(variables[k]))), weight );
        }
        else if(obj.has(variables[k])) {
          v_den[k]->Fill( obj.get(variables[k]), weight );
        }
      }
    }

  void fill_num( Object obj, double Nvtx, double PU, double Lumi, double weight = 1.0 ) {
      for(int k=0; k < nVar; ++k) {
        if( variables[k] == "nvtx" ) {
          v_num[k]->Fill( Nvtx, weight );
        }
        else if( variables[k] == "pu" ) {
          v_num[k]->Fill( PU, weight );
        }
        else if( variables[k] == "lumi" ) {
          v_num[k]->Fill( Lumi, weight );
        }
        else if( variables[k] == "mva" ) {
          v_num[k]->Fill( 1./(1.+exp(-1.*obj.get(variables[k]))), weight );
        }
        else if(obj.has(variables[k])) {
          v_num[k]->Fill( obj.get(variables[k]), weight );
        }
      }
    }

    void fill_den( TString varname, int ivar, double value, double weight = 1.0 ) {
      if( variables[ivar] == varname ) {
        v_den[ivar]->Fill( value, weight );
      }
    }

    void fill_num( TString varname, int ivar, double value, double weight = 1.0 ) {
      if( variables[ivar] == varname ) {
        v_num[ivar]->Fill( value, weight );
      }
    }

    void Save( TFile *f_output )
    {
      f_output->cd();

      for(int k=0; k < nVar; ++k) {
        v_den[k]->Write();
        v_num[k]->Write();
        delete v_den[k];
        delete v_num[k];
      }
    }

    void Save( TDirectory *dir )
    {
      dir->cd();

      for(int k=0; k < nVar; ++k) {
        v_den[k]->SetDirectory(dir);
        v_num[k]->SetDirectory(dir);
        v_den[k]->Write();
        v_num[k]->Write();
        delete v_den[k];
        delete v_num[k];
      }
    }

    ~HistContainer() {}

private:
    TString Tag;
    int nVar;
    vector<TString> variables;
    vector<vector<double>> ranges;

    vector<TH1D*> v_den;
    vector<TH1D*> v_num;

    void Init()
    {
      TH1::SetDefaultSumw2(kTRUE);
      TH2::SetDefaultSumw2(kTRUE);
      TH1::AddDirectory(kFALSE);
      TH2::AddDirectory(kFALSE);

      TString Tag_tmp = this->Tag == "" ? "" : "_"+this->Tag;

      for(int k=0; k < nVar; ++k) {
        TString name = TString::Format("%s_%s", Tag_tmp.Data(), variables[k].Data() );

        TH1D *den = new TH1D("den"+name, "", ranges[k][0], ranges[k][1], ranges[k][2]);
        TH1D *num = new TH1D("num"+name, "", ranges[k][0], ranges[k][1], ranges[k][2]);

        v_den.push_back( den );
        v_num.push_back( num );
      }
    }
};

class HistContainer2D
{
public:
    HistContainer2D(
      TString _Tag,
      vector<vector<TString>> _variables = {
        {"eta", "pt"},
        {"eta", "phi"},
      },
      vector<vector<double>> _ranges = {
        {24, -2.4, 2.4, 40, 0, 200},
        {24, -2.4, 2.4, 30, -TMath::Pi(), TMath::Pi()},
      }
    ) {

      if(_variables.size() != _ranges.size()) {
        cout << "HistContainer2D: _variables.size() != _ranges.size()" << endl;
        exit(1);
      }

      this->Tag = _Tag;
      this->variables = _variables;
      this->ranges = _ranges;
      this->nVar = variables.size();

      this->Init();
    }

  void fill_den( Object obj, double weight = 1.0 ) {
    for(int k=0; k < nVar; ++k) {
      if(obj.has(variables.at(k).at(0)) && obj.has(variables.at(k).at(1))) {
	v_den.at(k)->Fill( obj.get(variables.at(k).at(0)), obj.get(variables.at(k).at(1)), weight );
      }
    }
  }

  void fill_num( Object obj, double weight = 1.0 ) {
    for(int k=0; k < nVar; ++k) {
      if(obj.has(variables.at(k).at(0)) && obj.has(variables.at(k).at(1))) {
	v_num.at(k)->Fill( obj.get(variables.at(k).at(0)), obj.get(variables.at(k).at(1)), weight );
      }
    }
  }

  void Save( TFile *f_output )
  {
    f_output->cd();

    for(int k=0; k < nVar; ++k) {
      v_den.at(k)->Write();
      v_num.at(k)->Write();
      delete v_den.at(k);
      delete v_num.at(k);
    }
  }

  void Save( TDirectory *dir )
  {
    dir->cd();

    for(int k=0; k < nVar; ++k) {
      v_den.at(k)->SetDirectory(dir);
      v_num.at(k)->SetDirectory(dir);
      v_den.at(k)->Write();
      v_num.at(k)->Write();
      delete v_den.at(k);
      delete v_num.at(k);
    }
  }

  ~HistContainer2D() {}

private:
    TString Tag;
    int nVar;
    vector<vector<TString>> variables;
    vector<vector<double>> ranges;

    vector<TH2D*> v_den;
    vector<TH2D*> v_num;

    void Init()
    {
      TH1::SetDefaultSumw2(kTRUE);
      TH2::SetDefaultSumw2(kTRUE);
      TH1::AddDirectory(kFALSE);
      TH2::AddDirectory(kFALSE);

      TString Tag_tmp = this->Tag == "" ? "" : "_"+this->Tag;

      for(int k=0; k < nVar; ++k) {
        TString name = TString::Format("%s_%s_%s", Tag_tmp.Data(), variables.at(k).at(0).Data(), variables.at(k).at(1).Data());

        TH2D *den = new TH2D("den"+name, "", ranges.at(k)[0], ranges.at(k)[1], ranges.at(k)[2], ranges.at(k)[3], ranges.at(k)[4], ranges.at(k)[5]);
        TH2D *num = new TH2D("num"+name, "", ranges.at(k)[0], ranges.at(k)[1], ranges.at(k)[2], ranges.at(k)[3], ranges.at(k)[4], ranges.at(k)[5]);

        v_den.push_back( den );
        v_num.push_back( num );
      }
    }
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

// echo 'gROOT->LoadMacro("HLTEffAnalyzer.C+"); gSystem->Exit(0);' | root -b -l
// root -l -b -q 'HLTEffAnalyzer.C("v00", "TEST")' >&log&



void HLTEffAnalyzer(
    TString ver = "v00", TString tag = "TEST",
    vector<TString> vec_Dataset = {}, TString JobId = "",
    TString outputDir = "./",
    const bool doDimuon = false, double ZmassWindow = -1,
    Int_t maxEv = -1, bool doMem = false, int nMem = 10001, bool doBar = true  // HERE
    // Int_t maxEv = -1, bool doMem = false, int nMem = 10001, bool doBar = false
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
        TFile *f_output = TFile::Open(outputDir+fileName+"-Eff.root", "RECREATE");

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

        vector<TString> branch_tags = {
            "genParticle",
            "muon",
            "vec_",
            "L1Muon",
            "L2Muon",

            "hltPixelTracks",
            "hltIterL3OIMuonTrackAssociated",
            "hltIter0IterL3MuonTrackAssociated",
            "hltIterL3MuonMergedAssociated",
            "hltIter0IterL3FromL1MuonTrackAssociated",
            "hltIterL3MuonAndMuonFromL1MergedAssociated",
            "iterL3MuonNoID",
            "iterL3Muon",
        };

        unique_ptr<MuonHLTNtupleRun3>  nt( new MuonHLTNtupleRun3( _chain_Ev, branch_tags ) );

    // -- Histograms
        TH1D *h_nEvents = new TH1D("h_nEvents",  "", 3, -1, 2);

        TH1D *h_gen_pt  = new TH1D("h_gen_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_eta = new TH1D("h_gen_eta", "", 60, -3, 3);
        TH1D *h_gen_phi = new TH1D("h_gen_phi", "", 64, -3.2, 3.2);

        TH1D *h_gen_acc_pt  = new TH1D("h_gen_acc_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_acc_eta = new TH1D("h_gen_acc_eta", "", 60, -3, 3);
        TH1D *h_gen_acc_phi = new TH1D("h_gen_acc_phi", "", 64, -3.2, 3.2);

        TH1D *h_gen_hard_pt  = new TH1D("h_gen_hard_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_hard_eta = new TH1D("h_gen_hard_eta", "", 60, -3, 3);
        TH1D *h_gen_hard_phi = new TH1D("h_gen_hard_phi", "", 64, -3.2, 3.2);

        TH1D *h_gen_hard_acc_pt  = new TH1D("h_gen_hard_acc_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_hard_acc_eta = new TH1D("h_gen_hard_acc_eta", "", 60, -3, 3);
        TH1D *h_gen_hard_acc_phi = new TH1D("h_gen_hard_acc_phi", "", 64, -3.2, 3.2);

        vector<TString> L3types = {
            "L1Muon",
            "L2Muon",

            "hltPixelTracks",
            "hltPixelTracksInRegionL2",
            "hltPixelTracksInRegionL1",
            "hltOI",
            "hltIter0",
            "hltL3FromL2Merged",
            "hltIter0FromL1",
            "hltL3Merged",
            "hltIterL3MuonNoID",
            "hltIterL3Muon",
            "hltIterL3MuonNoIDTrack",
            "hltIterL3MuonTrack",

            "L1sSingleMu22",
            "Mu24",
            "Mu50",
            "Mu50OrOldMu100OrTkMu100",
            "ECALIsoMu24",
            "HCALIsoMu24",
            "IsoMu24",

            "myL1sSingleMu22",
            "myMu24",
            "myMu50",
            "myMu50OrOldMu100OrTkMu100",
            "myECALIsoMu24",
            "myHCALIsoMu24",
            "myIsoMu24",
        };

        vector<TString> HLTpaths = {
            "Mu24",
            "Mu50",
            "OldMu100",
            "TkMu100",
            "Mu50OrOldMu100OrTkMu100",
            "ECALIsoMu24",
            "HCALIsoMu24",
            "IsoMu24",
            "Mu17Mu8",
            "Mu37TkMu27"
        };

        // -- Efficiency
            vector<double> Eff_genpt_mins = {
                0,
                10,
                26,
                // 30,
                53,
                // 105,
            };

            vector<vector<double>> Etas_bin = {
                {0., 2.4},
            };
            vector<TString> Etas_str = {
                "I",
            };

            vector<vector<int>> Runs_bin = {
                {-1, 999999},
                {-1, 367661-1},
                {367661, 367990-1},
                {367990, 368765},
                {368765+1, 999999},
                {367905, 367990-1}, // 5 RUNs before new align - 367905, 367906, 367907, 367908, 367910
                {368223, 368320},   // 5 RUNs after new align  - 368223, 368224, 368229, 368318, 368320
            };

            vector<vector<vector<vector<HistContainer*>>>> hc_Eff = {};  // Eff[L3 type][run bin][eta bin][gen pt min]
            vector<vector<vector<vector<HistContainer*>>>> hc_Eff_L1SQ22 = {};  // Eff[L3 type][run bin][eta bin][gen pt min]
            vector<vector<vector<vector<HistContainer*>>>> hc_Eff_L1DQ8 = {};  // Eff[L3 type][run bin][eta bin][gen pt min]

            vector<vector<vector<vector<HistContainer2D*>>>> hc2D_Eff = {};  // Eff[L3 type][run bin][eta bin][gen pt min]
            vector<vector<vector<vector<HistContainer2D*>>>> hc2D_Eff_L1SQ22 = {};  // Eff[L3 type][run bin][eta bin][gen pt min]
            vector<vector<vector<vector<HistContainer2D*>>>> hc2D_Eff_L1DQ8 = {};  // Eff[L3 type][run bin][eta bin][gen pt min]

            // -- res
            vector<vector<TH1D *>> vh_L3_qbpt_pt  = {};
            vector<vector<TH1D *>> vh_L3_qbpt_eta = {};
            vector<vector<TH1D *>> vh_L3_pt_pt    = {};
            vector<vector<TH1D *>> vh_L3_pt_eta   = {};

            int iL3type = 0;
            for(auto& L3type: L3types) {

                hc_Eff.push_back( {} );
                hc_Eff_L1SQ22.push_back( {} );
                hc_Eff_L1DQ8.push_back( {} );
                hc2D_Eff.push_back( {} );
                hc2D_Eff_L1SQ22.push_back( {} );
                hc2D_Eff_L1DQ8.push_back( {} );

                for (unsigned irun = 0; irun < Runs_bin.size(); ++irun) {
                    TString run_str = TString::Format("Run%d_%d", Runs_bin.at(irun).at(0), Runs_bin.at(irun).at(1));
                    if (Runs_bin.at(irun).at(0) < 0 && Runs_bin.at(irun).at(1) == 999999)
                        run_str = "RunAll";

                    hc_Eff.at(iL3type).push_back( {} );
                    hc_Eff_L1SQ22.at(iL3type).push_back( {} );
                    hc_Eff_L1DQ8.at(iL3type).push_back( {} );
                    hc2D_Eff.at(iL3type).push_back( {} );
                    hc2D_Eff_L1SQ22.at(iL3type).push_back( {} );
                    hc2D_Eff_L1DQ8.at(iL3type).push_back( {} );

                    for (unsigned ieta = 0; ieta < Etas_bin.size(); ++ieta) {
                        hc_Eff.at(iL3type).at(irun).push_back( {} );
                        hc_Eff_L1SQ22.at(iL3type).at(irun).push_back( {} );
                        hc_Eff_L1DQ8.at(iL3type).at(irun).push_back( {} );
                        hc2D_Eff.at(iL3type).at(irun).push_back( {} );
                        hc2D_Eff_L1SQ22.at(iL3type).at(irun).push_back( {} );
                        hc2D_Eff_L1DQ8.at(iL3type).at(irun).push_back( {} );

                        for(auto& Eff_genpt_min: Eff_genpt_mins) {
                            HistContainer* hc_tmp0   = new HistContainer( TString::Format("Eff_%s_%s_%s_genpt%.0f",      L3type.Data(), run_str.Data(), Etas_str.at(ieta).Data(), Eff_genpt_min) );
                            HistContainer* hc_tmp1_2 = new HistContainer( TString::Format("Eff_L1SQ22_%s_%s_%s_genpt%.0f", L3type.Data(), run_str.Data(), Etas_str.at(ieta).Data(), Eff_genpt_min) );
                            HistContainer* hc_tmp2_1 = new HistContainer( TString::Format("Eff_L1DQ8_%s_%s_%s_genpt%.0f", L3type.Data(), run_str.Data(), Etas_str.at(ieta).Data(), Eff_genpt_min) );
                            HistContainer2D* hc2D_tmp0 = new HistContainer2D( TString::Format("2DEff_%s_%s_%s_genpt%.0f", L3type.Data(), run_str.Data(), Etas_str.at(ieta).Data(), Eff_genpt_min) );
                            HistContainer2D* hc2D_tmp1 = new HistContainer2D( TString::Format("2DEff_L1SQ22_%s_%s_%s_genpt%.0f", L3type.Data(), run_str.Data(), Etas_str.at(ieta).Data(), Eff_genpt_min) );
                            HistContainer2D* hc2D_tmp2 = new HistContainer2D( TString::Format("2DEff_L1DQ8_%s_%s_%s_genpt%.0f", L3type.Data(), run_str.Data(), Etas_str.at(ieta).Data(), Eff_genpt_min) );
                            hc_Eff.at(iL3type).at(irun).at(ieta).push_back( hc_tmp0 );
                            hc_Eff_L1SQ22.at(iL3type).at(irun).at(ieta).push_back( hc_tmp1_2 );
                            hc_Eff_L1DQ8.at(iL3type).at(irun).at(ieta).push_back( hc_tmp2_1 );
                            hc2D_Eff.at(iL3type).at(irun).at(ieta).push_back( hc2D_tmp0 );
                            hc2D_Eff_L1SQ22.at(iL3type).at(irun).at(ieta).push_back( hc2D_tmp1 );
                            hc2D_Eff_L1DQ8.at(iL3type).at(irun).at(ieta).push_back( hc2D_tmp2 );
                        }
                    }
                }

                // -- res
                vh_L3_qbpt_pt.push_back( {} );
                vh_L3_qbpt_eta.push_back( {} );
                vh_L3_pt_pt.push_back( {} );
                vh_L3_pt_eta.push_back( {} );

                for(int ipt=0; ipt<n_pt_bins; ++ipt) {
                    TString name = TString::Format("h_%s_qbpt_pt_%d", L3types.at(iL3type).Data(), ipt);
                    TH1D *h_L3_qbpt  = new TH1D(name,  "", 400, -1, 1);
                    vh_L3_qbpt_pt.at(iL3type).push_back( h_L3_qbpt );

                    name = TString::Format("h_%s_pt_pt_%d", L3types.at(iL3type).Data(), ipt);
                    TH1D *h_L3_pt  = new TH1D(name,  "", 400, -1, 1);
                    vh_L3_pt_pt.at(iL3type).push_back( h_L3_pt );
                }

                for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
                    TString name = TString::Format("h_%s_qbpt_eta_%d", L3types.at(iL3type).Data(), ieta);
                    TH1D *h_L3_qbpt  = new TH1D(name,  "", 400, -1, 1);
                    vh_L3_qbpt_eta.at(iL3type).push_back( h_L3_qbpt );

                    name = TString::Format("h_%s_pt_eta_%d", L3types.at(iL3type).Data(), ieta);
                    TH1D *h_L3_pt  = new TH1D(name,  "", 400, -1, 1);
                    vh_L3_pt_eta.at(iL3type).push_back( h_L3_pt );
                }

                iL3type += 1;
            }

    // -- Event loop
    for(unsigned i_ev=0; i_ev<nEvent; i_ev++) {
        //cout<<i_ev<<endl;
        if(doBar)
            loadBar(i_ev+1, nEvent, 100, 100);
        else if( doMem && i_ev !=0 && i_ev % nMem == 0 )
            printMemory("\t");
        else
            printRunTimeShort(timer_total);
            // cout << "\nEvent: " << i_ev << endl;  // HERE

        nt->GetEntry( i_ev );

        double genWeight = nt->genEventWeight > 0.0 ? 1.0 : -1.0;
        if (nt->isRealData)
            genWeight = 1.;
        h_nEvents->Fill( genWeight );

        /*
        vector<Object> GenParticles = nt->get_GenParticles();

        bool isDimuon = false;
        vector<Object> GenMuonsFromHardProcess = {};
        vector<Object> GenMuonsFromHPInAcc = {};
        // -- Dimuon skim for DY
            bool found0 = false;
            bool found1 = false;
            for(auto& genP: GenParticles) {

                if( fabs(genP.get("ID")) == 13 && genP.get("fromHardProcessFinalState") == 1 ) {
                    GenMuonsFromHardProcess.push_back( genP );

                    h_gen_hard_pt->Fill( genP.pt, genWeight );
                    h_gen_hard_eta->Fill( genP.eta, genWeight );
                    h_gen_hard_phi->Fill( genP.phi, genWeight );

                    if( acceptance( genP ) ) {
                        GenMuonsFromHPInAcc.push_back(genP);

                        h_gen_hard_acc_pt->Fill( genP.pt, genWeight );
                        h_gen_hard_acc_eta->Fill( genP.eta, genWeight );
                        h_gen_hard_acc_phi->Fill( genP.phi, genWeight );
                    }
                }

                if( fabs(genP.get("ID")) == 13 && genP.get("status") == 1 ) {
                    h_gen_pt->Fill( genP.pt, genWeight );
                    h_gen_eta->Fill( genP.eta, genWeight );
                    h_gen_phi->Fill( genP.phi, genWeight );

                    if( acceptance( genP ) ) {
                        h_gen_acc_pt->Fill( genP.pt, genWeight );
                        h_gen_acc_eta->Fill( genP.eta, genWeight );
                        h_gen_acc_phi->Fill( genP.phi, genWeight );
                    }
                }

                if( genP.get("ID") == 13 && genP.get("isHardProcess") == 1 )
                    found0 = true;
                if( genP.get("ID") == -13 && genP.get("isHardProcess") == 1 )
                    found1 = true;
            }
            isDimuon = (found0 && found1);
            if (doDimuon && !isDimuon) {
                continue;
            }

            std::sort(GenMuonsFromHardProcess.begin(), GenMuonsFromHardProcess.end(), sort_by_pt());
            std::sort(GenMuonsFromHPInAcc.begin(), GenMuonsFromHPInAcc.end(), sort_by_pt());

            if(doDimuon && isDimuon && ZmassWindow > 0.) {
                TLorentzVector mu0, mu1;
                mu0.SetPtEtaPhiM(GenMuonsFromHardProcess.at(0).pt,
                                 GenMuonsFromHardProcess.at(0).eta,
                                 GenMuonsFromHardProcess.at(0).phi,
                                 MU_MASS);
                mu1.SetPtEtaPhiM(GenMuonsFromHardProcess.at(1).pt,
                                 GenMuonsFromHardProcess.at(1).eta,
                                 GenMuonsFromHardProcess.at(1).phi,
                                 MU_MASS);
                double Zmass = (mu0+mu1).M();
                if (abs(Zmass-91) > ZmassWindow) {
                    continue;
                }
            }
        */

        // -- Get object collections
            vector<Object> L1Muons = nt->get_L1Muons();
            vector<Object> L2Muons = nt->get_L2Muons();

            vector<Object> hltPixelTracksAssociated = nt->get_hltPixelTracksAssociated();
            vector<Object> hltPixelTracksInRegionL2Associated = nt->get_hltPixelTracksInRegionL2Associated();
            vector<Object> hltPixelTracksInRegionL1Associated = nt->get_hltPixelTracksInRegionL1Associated();
            vector<Object> hltIterL3OIMuonTrackAssociated = nt->get_hltIterL3OIMuonTrackAssociated();
            vector<Object> hltIter0IterL3MuonTrackAssociated = nt->get_hltIter0IterL3MuonTrackAssociated();
            vector<Object> hltIterL3MuonMergedAssociated = nt->get_hltIterL3MuonMergedAssociated();
            vector<Object> hltIter0IterL3FromL1MuonTrackAssociated = nt->get_hltIter0IterL3FromL1MuonTrackAssociated();
            vector<Object> hltIterL3MuonAndMuonFromL1MergedAssociated = nt->get_hltIterL3MuonAndMuonFromL1MergedAssociated();
            vector<Object> iterL3MuonNoID = nt->get_iterL3MuonNoID();
            vector<Object> iterL3Muon = nt->get_iterL3Muon();
            vector<Object> iterL3MuonNoIDTrackAssociated = nt->get_iterL3MuonNoIDTrackAssociated();
            vector<Object> iterL3MuonTrackAssociated = nt->get_iterL3MuonTrackAssociated();

            vector<Object> L1sSingleMu22_HLT = nt->get_HLTObjects("hltL1fL1sMu22L1Filtered0");
            vector<Object> Mu24_HLT = nt->get_HLTObjects("hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q");
            vector<Object> Mu50_HLT = nt->get_HLTObjects("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q");
            vector<Object> TkMu100_HLT = nt->get_HLTObjects("hltL3fL1sMu25f0TkFiltered100Q");
            vector<Object> OldMu100_HLT = nt->get_HLTObjects("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q");
            vector<Object> ECALIsoMu24_HLT = nt->get_HLTObjects("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfecalIsoRhoFiltered");
            vector<Object> HCALIsoMu24_HLT = nt->get_HLTObjects("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfhcalIsoRhoFiltered");
            vector<Object> IsoMu24_HLT = nt->get_HLTObjects("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered");

            vector<Object> L1sSingleMu22_MYHLT = nt->get_myHLTObjects("hltL1fL1sMu22L1Filtered0");
            vector<Object> Mu24_MYHLT = nt->get_myHLTObjects("hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q");
            vector<Object> Mu50_MYHLT = nt->get_myHLTObjects("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q");
            vector<Object> TkMu100_MYHLT = nt->get_myHLTObjects("hltL3fL1sMu25f0TkFiltered100Q");
            vector<Object> OldMu100_MYHLT = nt->get_myHLTObjects("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q");
            vector<Object> ECALIsoMu24_MYHLT = nt->get_myHLTObjects("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfecalIsoRhoFiltered");
            vector<Object> HCALIsoMu24_MYHLT = nt->get_myHLTObjects("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfhcalIsoRhoFiltered");
            vector<Object> IsoMu24_MYHLT = nt->get_myHLTObjects("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered");

        // -- TnP selection
            vector<Object> muons = nt->get_offlineMuons();
            vector<Object> probes = {};
            // tag
            for (auto & i_mu : muons) {
                if (i_mu.pt < 27.)
                    continue;
                if (!offlineSel(i_mu))
                    continue;
                if (!i_mu.matched(IsoMu24_HLT, 0.1))
                    continue;

                // probe
                for (auto & j_mu : muons) {
                    if (!offlineSel(j_mu))
                        continue;
                    if (i_mu.get("charge") * j_mu.get("charge") > 0)
                        continue;
                    double pair_mass = invMass(i_mu, j_mu);
                    if (pair_mass < 81.)
                        continue;
                    if (pair_mass > 101.)
                        continue;

                    j_mu.addVar("pair_mass", pair_mass);

                    probes.push_back(j_mu);
                }
            }

            vector<vector<Object>*> L3MuonColls {
                &probes,  // for L1 muon eff
                &L2Muons,

                &hltPixelTracksAssociated,
                &hltPixelTracksInRegionL2Associated,
                &hltPixelTracksInRegionL1Associated,
                &hltIterL3OIMuonTrackAssociated,
                &hltIter0IterL3MuonTrackAssociated,
                &hltIterL3MuonMergedAssociated,
                &hltIter0IterL3FromL1MuonTrackAssociated,
                &hltIterL3MuonAndMuonFromL1MergedAssociated,
                &iterL3MuonNoID,
                &iterL3Muon,
                &iterL3MuonNoIDTrackAssociated,
                &iterL3MuonTrackAssociated,

                &L1sSingleMu22_HLT,
                &Mu24_HLT,
                &Mu50_HLT,
                &Mu50_HLT,// Mu50OrOldMu100OrTkMu100
                &ECALIsoMu24_HLT,
                &HCALIsoMu24_HLT,
                &IsoMu24_HLT,

                &L1sSingleMu22_MYHLT,
                &Mu24_MYHLT,
                &Mu50_MYHLT,
                &Mu50_MYHLT, // Mu50OrOldMu100OrTkMu100
                &ECALIsoMu24_MYHLT,
                &HCALIsoMu24_MYHLT,
                &IsoMu24_MYHLT
            };

            if (L3types.size() != L3MuonColls.size()) {
                cout << "ERROR: L3types.size() != L3MuonColls.size()" << endl;
                return;
            }

        //### L3types loop ###
        for(unsigned i=0; i<L3types.size(); ++i) {
            vector<Object>* L3Coll = L3MuonColls.at(i);
            TString L3type = L3types.at(i);
            L3type.ReplaceAll("my","");
            bool looseMatch = L3type.Contains("L2Muon");

            // -- Efficiency
            for (unsigned irun = 0; irun < Runs_bin.size(); ++irun) {
                if (nt->runNum < Runs_bin.at(irun).at(0))
                    continue;
                if (nt->runNum > Runs_bin.at(irun).at(1))
                    continue;
                for (unsigned ieta = 0; ieta < Etas_bin.size(); ++ieta) {
                    //### offlineMuons loop ###
                    int iprobe = -1;
                    for(auto& probemu: probes) {
                        L3Coll = L3MuonColls.at(i);

                        if( !acceptance( probemu ) )
                            continue;
                        if (Etas_bin.at(ieta).at(0) > abs(probemu.eta))
                            continue;
                        if (Etas_bin.at(ieta).at(1) < abs(probemu.eta))
                            continue;

                        // New L1s matching
                        /*
                        //if (probemu.pt < 26)
                        //    continue;
                        iprobe++;
                        if(probemu.get("l1ptByQ") < 1.0){
                          cout<<"\n\nOffline mu("<<i<<")        : pt = "<<probemu.pt<<", eta = "<<probemu.eta<<", phi = "<<probemu.phi<<", charge = "<<probemu.get("charge")<<endl;
                          cout<<"recol1MatchByQ       : pt = "<<probemu.get("l1ptByQ")<<", eta = "<<probemu.get("l1etaByQ")<<", phi = "<<probemu.get("l1phiByQ")<<", qual = "<<probemu.get("l1qByQ")<<", charge = "<<probemu.get("l1chargeByQ")<<endl;
                            for(unsigned il1=0; il1<L1Muons.size(); il1++){
                              cout<<"L1 Collection's "<<il1<<"th  : pt = "<<L1Muons.at(il1).pt<<", eta = "<<L1Muons.at(il1).eta<<", phi = "<<L1Muons.at(il1).phi<<", dR = "<<probemu.deltaR(L1Muons.at(il1))<<", charge = "<<L1Muons.at(il1).get("charge")<<", qual = "<<L1Muons.at(il1).get("quality")<<endl;
                            }
                            cout<<"nl1t = "<<probemu.get("nl1t")<<endl;
                            cout<<"l1t size = "<<probemu.getvec("l1tpt").size()<<endl;
                            for(unsigned il1t=0; il1t<probemu.get("nl1t"); il1t++){
                              cout<<"L1T Collection's "<<il1t<<"th : pt = "<<probemu.getvec("l1tpt").at(il1t)<<", eta = "<<probemu.getvec("l1teta").at(il1t)<<", phi = "<<probemu.getvec("l1tphi").at(il1t)<<", dR = "<<probemu.getvec("l1tdr").at(il1t)<<", charge = "<<probemu.getvec("l1tcharge").at(il1t)<<", qual = "<<probemu.getvec("l1tq").at(il1t)<<endl;
                            }
                        }
                        */
                        bool matched_L1SQ0 = (
                            probemu.get("l1ptByQ") > -1.0 &&
                            probemu.get("l1drByQ") < 0.3 &&
                            probemu.get("l1qByQ") > 11
                        );
                        bool matched_L1SQ8 = (
                            probemu.get("l1ptByQ") > 8.0 &&
                            probemu.get("l1drByQ") < 0.3 &&
                            probemu.get("l1qByQ") > 11
                        );
                        bool matched_L1SQ22 = (
                            probemu.get("l1ptByQ") > 22.0 &&
                            probemu.get("l1drByQ") < 0.3 &&
                            probemu.get("l1qByQ") > 11
                        );

                        bool matched_L1DQ0 = (
                            probemu.get("l1ptByQ") > -1.0 &&
                            probemu.get("l1drByQ") < 0.3 &&
                            probemu.get("l1qByQ") > 7
                        );
                        bool matched_L1DQ8 = (
                            probemu.get("l1ptByQ") > 8.0 &&
                            probemu.get("l1drByQ") < 0.3 &&
                            probemu.get("l1qByQ") > 7
                        );
                        bool matched_L1DQ22 = (
                            probemu.get("l1ptByQ") > 22.0 &&
                            probemu.get("l1drByQ") < 0.3 &&
                            probemu.get("l1qByQ") > 7
                        );

                        bool l1matched_L1SQ0 = probemu.l1matched(0., 11);
                        bool l1matched_L1SQ8 = probemu.l1matched(8., 11);
                        bool l1matched_L1SQ22 = probemu.l1matched(22., 11);
                        bool l1matched_L1DQ0 = probemu.l1matched(0., 7);
                        bool l1matched_L1DQ8 = probemu.l1matched(8., 7);
                        bool l1matched_L1DQ22 = probemu.l1matched(22., 7);

                        //cout<<"matched = "<<(matched_L1SQ22?"true":"false")<<", l1matched = "<<(m1atched_L1SQ22?"true":"false")<<endl;
                        int matched_idx = -1e6;

                        // HERE !!!
                        vector<int> L3map(L3Coll->size(), -1);
                        if (L3type.Contains("L1Muon")) {
                            matched_idx = -1e6;
                        }
                        else if (
                            L3type.Contains("OI") ||
                            L3type.Contains("L3") ||
                            L3type.Contains("GlbTrkMuon") ||
                            (std::find(HLTpaths.begin(), HLTpaths.end(), L3type) != HLTpaths.end())
                        ) {
                            matched_idx = probemu.matched( *L3Coll, L3map, 0.1 );
                        }
                        else if (
                            L3type.Contains("hltPixelTracks")
                        ) {
                            matched_idx = probemu.matched( *L3Coll, L3map, 0.01 );  // Pixel tracks - TRK guys are using dR 0.01 btw Offline, HLT
                        }
                        else if (
                            L3type.Contains("L1sSingleMu22")
                        ) {
                            matched_idx = probemu.matched( *L3Coll, L3map, 0.3 );
                        }
                        else {
                            matched_idx = looseMatch ? probemu.matched( *L3Coll, L3map, 0.3 ) :  // L2 muon
                                                       probemu.matched( *L3Coll, L3map, 0.1, 0.5 );  // IO tracks
                        }

                        // Mu50OrOldMu100OrTkMu100
                        if (L3type.Contains("Mu50OrOldMu100OrTkMu100") &&
                            matched_idx < 0) {
                            vector<int> TkMu100map(L3types.at(i).Contains("my")? TkMu100_MYHLT.size() : TkMu100_HLT.size(), -1);
                            matched_idx = probemu.matched(L3types.at(i).Contains("my")? TkMu100_MYHLT : TkMu100_HLT, TkMu100map, 0.1);
                            if (matched_idx >= 0) L3Coll = L3types.at(i).Contains("my")? &TkMu100_MYHLT : &TkMu100_HLT;
                            if (matched_idx < 0) {
                                vector<int> OldMu100map(L3types.at(i).Contains("my")? OldMu100_MYHLT.size() : OldMu100_HLT.size(), -1);
                                matched_idx = probemu.matched(L3types.at(i).Contains("my")? OldMu100_MYHLT : OldMu100_HLT, OldMu100map, 0.1);
                                if (matched_idx >= 0) L3Coll = L3types.at(i).Contains("my")? &OldMu100_MYHLT : &OldMu100_HLT;
                            }
                        }

                        if (matched_idx < 0) {
                            probemu.addVar("dR", -1.);
                        } else {
                            probemu.addVar("dR", probemu.deltaR(L3Coll->at(matched_idx)));
                        }

                        int matched_idx_res = -1e6;
                        vector<int> L3map2(L3Coll->size(), -1);
                        matched_idx_res = looseMatch ? probemu.matched( *L3Coll, L3map2, 0.3 ) :  // L2 muon
                                                       probemu.matched( *L3Coll, L3map2, 0.1 );

                        // --  Efficiency / Gen or L1
                        for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                            if( probemu.pt > Eff_genpt_mins.at(j) ) {
                                hc_Eff.at(i).at(irun).at(ieta).at(j)->fill_den( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                hc2D_Eff.at(i).at(irun).at(ieta).at(j)->fill_den( probemu, genWeight );
                                if( matched_idx > -1 ) {
                                    hc_Eff.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                    hc2D_Eff.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, genWeight );
                                }

                                if(l1matched_L1DQ8) {
                                    hc_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->fill_den( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                    hc2D_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->fill_den( probemu, genWeight );

                                    if(L3type.Contains("L1Muon")) {
                                        hc_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                        hc2D_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, genWeight );
                                    }
                                    else {
                                        if( matched_idx > -1 ) {
                                            hc_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                            hc2D_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, genWeight );
                                        }
                                    }
                                }

                                if(l1matched_L1SQ22) {
                                    hc_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->fill_den( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                    hc2D_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->fill_den( probemu, genWeight );

                                    if(L3type.Contains("L1Muon")) {
                                        hc_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                        hc2D_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, genWeight );
                                    }
                                    else {
                                        if( matched_idx > -1 ) {
                                            hc_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, nt->nVertex, nt->DataPU, nt->InstLumi, genWeight );
                                            hc2D_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->fill_num( probemu, genWeight );
                                        }
                                    }
                                }
                            }
                        }

                        // -- res
                        if (matched_idx_res > -1 &&
                            (std::find(HLTpaths.begin(), HLTpaths.end(), L3type) == HLTpaths.end()) &&
                            Runs_bin.at(irun).at(0) < 0 &&
                            Etas_bin.at(ieta).at(0) == 0. &&
                            Etas_bin.at(ieta).at(1) == 2.4
                        ) {
                            if (L3type.Contains("L1Muon"))
                                continue;

                            if (!L3Coll->at(matched_idx_res).has("charge")) {
                                if (L3Coll->at(matched_idx_res).has("inner_charge")) {
                                    L3Coll->at(matched_idx_res).addVar("charge", L3Coll->at(matched_idx_res).get("inner_charge"));
                                }
                                else {
                                    L3Coll->at(matched_idx_res).addVar("charge", 0.);
                                }
                            }

                            double L3charge  = L3Coll->at(matched_idx_res).get("charge");
                            double gencharge = probemu.get("charge");
                            if (L3charge == 0.) {
                                L3charge = 1.;
                                gencharge = 1.;
                            }
                            double L3pt      = L3Coll->at(matched_idx_res).pt;
                            double L3qbpt    = L3charge / L3Coll->at(matched_idx_res).pt;
                            double genpt     = probemu.pt;
                            double genqbpt   = gencharge / probemu.pt;
                            double res_qbpt  = (genqbpt - L3qbpt) / genqbpt;
                            double res_pt    = (L3pt - genpt) / genpt;

                            for(int ipt=0; ipt<n_pt_bins-1; ++ipt) {
                                if( probemu.pt >= pt_bins[ipt] && probemu.pt < pt_bins[ipt+1] ) {
                                    vh_L3_qbpt_pt.at(i)[ipt]->Fill( res_qbpt, genWeight );
                                    vh_L3_pt_pt.at(i)[ipt]->Fill(   res_pt,   genWeight );
                                    break;
                                }
                            }
                            if (probemu.pt > 26.0) {
                                vh_L3_qbpt_pt.at(i)[n_pt_bins-1]->Fill( res_qbpt, genWeight );
                                vh_L3_pt_pt.at(i)[n_pt_bins-1]->Fill(   res_pt,   genWeight );
                            }
                            for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
                                if( probemu.eta >= eta_bins[ieta] && probemu.eta < eta_bins[ieta+1] && probemu.pt > 26.0 ) {
                                    vh_L3_qbpt_eta.at(i)[ieta]->Fill( res_qbpt, genWeight );
                                    vh_L3_pt_eta.at(i)[ieta]->Fill(   res_pt,   genWeight );
                                    break;
                                }
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

        h_gen_pt->Write();
        h_gen_eta->Write();
        h_gen_phi->Write();

        h_gen_acc_pt->Write();
        h_gen_acc_eta->Write();
        h_gen_acc_phi->Write();

        h_gen_hard_pt->Write();
        h_gen_hard_eta->Write();
        h_gen_hard_phi->Write();

        h_gen_hard_acc_pt->Write();
        h_gen_hard_acc_eta->Write();
        h_gen_hard_acc_phi->Write();

        TDirectory* dir0 = f_output->mkdir("Eff");
        dir0->cd();

        for(unsigned i=0; i<L3types.size(); ++i) {
            TDirectory* dir1 = dir0->mkdir(L3types.at(i));
            dir1->cd();

            for (unsigned irun = 0; irun < Runs_bin.size(); ++irun) {
                for (unsigned ieta = 0; ieta < Etas_bin.size(); ++ieta) {
                    for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                        hc_Eff.at(i).at(irun).at(ieta).at(j)->Save(dir1);
                        hc_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->Save(dir1);
                        hc_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->Save(dir1);
                        hc2D_Eff.at(i).at(irun).at(ieta).at(j)->Save(dir1);
                        hc2D_Eff_L1SQ22.at(i).at(irun).at(ieta).at(j)->Save(dir1);
                        hc2D_Eff_L1DQ8.at(i).at(irun).at(ieta).at(j)->Save(dir1);
                    }
                }
            }
            dir0->cd();
        }

        TDirectory* dir1 = f_output->mkdir("Res");
        dir1->cd();

        for(unsigned i=0; i<L3types.size(); ++i) {
            for(int ipt=0; ipt<n_pt_bins; ++ipt) {
                vh_L3_qbpt_pt.at(i)[ipt]->Write();
                vh_L3_pt_pt.at(i)[ipt]->Write();
            }
            for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
                vh_L3_qbpt_eta.at(i)[ieta]->Write();
                vh_L3_pt_eta.at(i)[ieta]->Write();
            }
        }

        f_output->cd();
        f_output->Close();

    printRunTime(timer_total);
}


