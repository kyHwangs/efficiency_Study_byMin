#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

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
      vector<TString> _variables = { "pt", "eta", "phi", "pu" },
      vector<vector<double>> _ranges = {
        { 1000, 0, 1000 },
        { 48, -2.4, 2.4 },
        { 60, -TMath::Pi(), TMath::Pi() },
        { 250, 0, 250}
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

    void fill_den( Object obj, double PU, double weight = 1.0 ) {
      for(int k=0; k < nVar; ++k) {
        if( variables[k] == "pu" ) {
          v_den[k]->Fill( PU, weight );
        }
        else if(obj.has(variables[k])) {
          v_den[k]->Fill( obj.get(variables[k]), weight );
        }
      }
    }

    void fill_num( Object obj, double PU, double weight = 1.0 ) {
      for(int k=0; k < nVar; ++k) {
        if( variables[k] == "pu" ) {
          v_num[k]->Fill( PU, weight );
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

struct sort_by_pt
{
    inline bool operator() (const Object& a, const Object& b)
    {
        return (a.pt > b.pt);
    }
};

bool acceptance( Object obj )
{
    return ( fabs(obj.eta) < 2.4 );
}

// echo 'gROOT->LoadMacro("HLTEffAnalyzer.C+"); gSystem->Exit(0);' | root -b -l
// root -l -b -q 'HLTEffAnalyzer.C("v00", "TEST")' >&log&

const double MU_MASS = 0.1056583745;

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
            paths = { "./ntuple*.root" };
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
            "vec_my",
            "L1Muon",
            "L2Muon",
            "iterL3OI",
            "iterL3IOFromL2",
            "iterL3FromL2",
            "iterL3IOFromL1",
            "iterL3MuonNoID",
            "iterL3Muon",
            "hltIter"
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
            "L3OI",
            "L3IOFromL2",
            "L3IOFromL1",
            "L3MuonNoId",
            "L3Muon",

            "hltOI",
            "hltIter0",
            "hltIter2",
            "hltIOFromL2Merged",
            "hltL3FromL2Merged",
            "hltIter0FromL1",
            "hltIter2FromL1",
            "hltFromL1Merged",
            "hltL3Merged",
            "hltIterL3MuonNoID",
            "hltIterL3Muon"
        };

        // -- Efficiency
            vector<double> Eff_genpt_mins = {
                0,
                10,
                26  //,
                // 53
            };
            vector<double> Eff_L3pt_mins = {
                // 10,
                // 22,
                // 24,
                // 50
            };

            vector<vector<HistContainer*>> hc_Eff = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1SQ0 = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1SQ8 = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1SQ22 = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1DQ0 = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1DQ8 = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1DQ22 = {};  // Eff[L3 type][gen pt min]

            vector<vector<HistContainer*>> hc_EffTO = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1SQ0 = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1SQ8 = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1SQ22 = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1DQ0 = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1DQ8 = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1DQ22 = {};  // Eff[L3 type][L3 pt min]

            // -- res
            vector<vector<TH1D *>> vh_L3_qbpt_pt  = {};
            vector<vector<TH1D *>> vh_L3_qbpt_eta = {};
            vector<vector<TH1D *>> vh_L3_pt_pt    = {};
            vector<vector<TH1D *>> vh_L3_pt_eta   = {};

            int iL3type = 0;
            for(auto& L3type: L3types) {

                hc_Eff.push_back( {} );
                hc_Eff_L1SQ0.push_back( {} );
                hc_Eff_L1SQ8.push_back( {} );
                hc_Eff_L1SQ22.push_back( {} );
                hc_Eff_L1DQ0.push_back( {} );
                hc_Eff_L1DQ8.push_back( {} );
                hc_Eff_L1DQ22.push_back( {} );
                hc_EffTO.push_back( {} );
                hc_EffTO_L1SQ0.push_back( {} );
                hc_EffTO_L1SQ8.push_back( {} );
                hc_EffTO_L1SQ22.push_back( {} );
                hc_EffTO_L1DQ0.push_back( {} );
                hc_EffTO_L1DQ8.push_back( {} );
                hc_EffTO_L1DQ22.push_back( {} );

                for(auto& Eff_genpt_min: Eff_genpt_mins) {
                    HistContainer* hc_tmp0   = new HistContainer( TString::Format("Eff_%s_genpt%.0f",      L3type.Data(), Eff_genpt_min) );
                    HistContainer* hc_tmp1_0 = new HistContainer( TString::Format("Eff_L1SQ0_%s_genpt%.0f", L3type.Data(), Eff_genpt_min) );
                    HistContainer* hc_tmp1_1 = new HistContainer( TString::Format("Eff_L1SQ8_%s_genpt%.0f", L3type.Data(), Eff_genpt_min) );
                    HistContainer* hc_tmp1_2 = new HistContainer( TString::Format("Eff_L1SQ22_%s_genpt%.0f", L3type.Data(), Eff_genpt_min) );
                    HistContainer* hc_tmp2_0 = new HistContainer( TString::Format("Eff_L1DQ0_%s_genpt%.0f", L3type.Data(), Eff_genpt_min) );
                    HistContainer* hc_tmp2_1 = new HistContainer( TString::Format("Eff_L1DQ8_%s_genpt%.0f", L3type.Data(), Eff_genpt_min) );
                    HistContainer* hc_tmp2_2 = new HistContainer( TString::Format("Eff_L1DQ22_%s_genpt%.0f", L3type.Data(), Eff_genpt_min) );
                    hc_Eff.at(iL3type).push_back( hc_tmp0 );
                    hc_Eff_L1SQ0.at(iL3type).push_back( hc_tmp1_0 );
                    hc_Eff_L1SQ8.at(iL3type).push_back( hc_tmp1_1 );
                    hc_Eff_L1SQ22.at(iL3type).push_back( hc_tmp1_2 );
                    hc_Eff_L1DQ0.at(iL3type).push_back( hc_tmp2_0 );
                    hc_Eff_L1DQ8.at(iL3type).push_back( hc_tmp2_1 );
                    hc_Eff_L1DQ22.at(iL3type).push_back( hc_tmp2_2 );
                }

                for(auto& Eff_L3pt_min: Eff_L3pt_mins) {
                    HistContainer* hc_tmp0   = new HistContainer( TString::Format("EffTO_%s_L3pt%.0f",      L3type.Data(), Eff_L3pt_min) );
                    HistContainer* hc_tmp1_0 = new HistContainer( TString::Format("EffTO_L1SQ0_%s_L3pt%.0f", L3type.Data(), Eff_L3pt_min) );
                    HistContainer* hc_tmp1_1 = new HistContainer( TString::Format("EffTO_L1SQ8_%s_L3pt%.0f", L3type.Data(), Eff_L3pt_min) );
                    HistContainer* hc_tmp1_2 = new HistContainer( TString::Format("EffTO_L1SQ22_%s_L3pt%.0f", L3type.Data(), Eff_L3pt_min) );
                    HistContainer* hc_tmp2_0 = new HistContainer( TString::Format("EffTO_L1DQ0_%s_L3pt%.0f", L3type.Data(), Eff_L3pt_min) );
                    HistContainer* hc_tmp2_1 = new HistContainer( TString::Format("EffTO_L1DQ8_%s_L3pt%.0f", L3type.Data(), Eff_L3pt_min) );
                    HistContainer* hc_tmp2_2 = new HistContainer( TString::Format("EffTO_L1DQ22_%s_L3pt%.0f", L3type.Data(), Eff_L3pt_min) );
                    hc_EffTO.at(iL3type).push_back( hc_tmp0 );
                    hc_EffTO_L1SQ0.at(iL3type).push_back( hc_tmp1_0 );
                    hc_EffTO_L1SQ8.at(iL3type).push_back( hc_tmp1_1 );
                    hc_EffTO_L1SQ22.at(iL3type).push_back( hc_tmp1_2 );
                    hc_EffTO_L1DQ0.at(iL3type).push_back( hc_tmp2_0 );
                    hc_EffTO_L1DQ8.at(iL3type).push_back( hc_tmp2_1 );
                    hc_EffTO_L1DQ22.at(iL3type).push_back( hc_tmp2_2 );
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
        if(doBar)
            loadBar(i_ev+1, nEvent, 100, 100);
        else if( doMem && i_ev !=0 && i_ev % nMem == 0 )
            printMemory("\t");
        else
            printRunTimeShort(timer_total);
            // cout << "\nEvent: " << i_ev << endl;  // HERE

        nt->GetEntry( i_ev );

        double genWeight = nt->genEventWeight > 0.0 ? 1.0 : -1.0;
        h_nEvents->Fill( genWeight );

        vector<Object> GenParticles = nt->get_GenParticles();

        bool isDimuon = false;
        vector<Object> GenMuonsFromHardProcess = {};
        vector<Object> GenMuonsInAcc = {};
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

        // -- Get object collections
            // vector<Object> L1Muons              = nt->get_L1Muons();
            vector<Object> L2Muons              = nt->get_L2Muons();

            vector<Object> iterL3OI             = nt->get_iterL3OI();
            vector<Object> iterL3IOFromL2       = nt->get_iterL3IOFromL2();
            vector<Object> iterL3IOFromL1       = nt->get_iterL3IOFromL1();
            vector<Object> iterL3MuonNoID       = nt->get_iterL3MuonNoID();
            vector<Object> iterL3Muon           = nt->get_iterL3Muon();

            vector<Object> hltOI                = nt->get_hltIterL3OIMuonTrackAssociated();
            vector<Object> hltIter0             = nt->get_hltIter0IterL3MuonTrackAssociated();
            vector<Object> hltIter2             = nt->get_hltIter2IterL3MuonTrackAssociated();
            vector<Object> hltIOFromL2Merged    = nt->get_hltIter2IterL3MuonMergedAssociated();
            vector<Object> hltL3FromL2Merged    = nt->get_hltIterL3MuonMergedAssociated();
            vector<Object> hltIter0FromL1       = nt->get_hltIter0IterL3FromL1MuonTrackAssociated();
            vector<Object> hltIter2FromL1       = nt->get_hltIter2IterL3FromL1MuonTrackAssociated();
            vector<Object> hltFromL1Merged      = nt->get_hltIter2IterL3FromL1MuonMergedAssociated();
            vector<Object> hltL3Merged          = nt->get_hltIterL3MuonAndMuonFromL1MergedAssociated();
            vector<Object> hltIterL3MuonNoID    = nt->get_iterL3MuonNoIDTrackAssociated();
            vector<Object> hltIterL3Muon        = nt->get_iterL3MuonTrackAssociated();

            vector<vector<Object>*> L3MuonColls {
                &GenMuonsFromHardProcess,  // for L1 muon eff
                &L2Muons,
                &iterL3OI,
                &iterL3IOFromL2,
                &iterL3IOFromL1,
                &iterL3MuonNoID,
                &iterL3Muon,

                &hltOI,
                &hltIter0,
                &hltIter2,
                &hltIOFromL2Merged,
                &hltL3FromL2Merged,
                &hltIter0FromL1,
                &hltIter2FromL1,
                &hltFromL1Merged,
                &hltL3Merged,
                &hltIterL3MuonNoID,
                &hltIterL3Muon
            };

        for(unsigned i=0; i<L3types.size(); ++i) {
            vector<Object>* L3Coll = L3MuonColls.at(i);

            bool looseMatch = L3types.at(i).Contains("L2Muon");

            // -- Efficiency
            if( !doDimuon || (doDimuon && isDimuon) ) {

                for(auto& genmu: GenParticles) {
                    if( fabs(genmu.get("ID")) != 13 )
                        continue;

                    if( !acceptance( genmu ) )
                        continue;

                    if( genmu.get("status") != 1 )
                        continue;

                    bool fromHardProcess = genmu.matched(GenMuonsFromHardProcess, 0.001);
                    if( doDimuon && !fromHardProcess )
                        continue;

                    bool matched_L1SQ0 = (
                        genmu.get("l1ptByQ") > -1.0 &&
                        genmu.get("l1drByQ") < 0.3 &&
                        genmu.get("l1qByQ") > 11
                    );
                    bool matched_L1SQ8 = (
                        genmu.get("l1ptByQ") > 8.0 &&
                        genmu.get("l1drByQ") < 0.3 &&
                        genmu.get("l1qByQ") > 11
                    );
                    bool matched_L1SQ22 = (
                        genmu.get("l1ptByQ") > 22.0 &&
                        genmu.get("l1drByQ") < 0.3 &&
                        genmu.get("l1qByQ") > 11
                    );

                    bool matched_L1DQ0 = (
                        genmu.get("l1ptByQ") > -1.0 &&
                        genmu.get("l1drByQ") < 0.3 &&
                        genmu.get("l1qByQ") > 7
                    );
                    bool matched_L1DQ8 = (
                        genmu.get("l1ptByQ") > 8.0 &&
                        genmu.get("l1drByQ") < 0.3 &&
                        genmu.get("l1qByQ") > 7
                    );
                    bool matched_L1DQ22 = (
                        genmu.get("l1ptByQ") > 22.0 &&
                        genmu.get("l1drByQ") < 0.3 &&
                        genmu.get("l1qByQ") > 7
                    );

                    int matched_idx = -1e6;

                    // HERE !!!
                    vector<int> L3map(L3Coll->size(), -1);
                    if (L3types.at(i).Contains("OI") || L3types.at(i).Contains("L3Muon")) {
                        matched_idx = genmu.matched( *L3Coll, L3map, 0.1 );
                    }
                    else {
                        matched_idx = looseMatch ? genmu.matched( *L3Coll, L3map, 0.3 ) :  // L2 muon
                                                   genmu.matched( *L3Coll, L3map, 0.1, 0.5 );  // IO tracks
                    }

                    int matched_idx_res = -1e6;
                    vector<int> L3map2(L3Coll->size(), -1);
                    matched_idx_res = looseMatch ? genmu.matched( *L3Coll, L3map2, 0.3 ) :  // L2 muon
                                                   genmu.matched( *L3Coll, L3map2, 0.1 );

                    // --  Efficiency / Gen or L1
                    for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                        if( genmu.pt > Eff_genpt_mins.at(j) ) {
                            hc_Eff.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if(matched_L1SQ0) {
                                hc_Eff_L1SQ0.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                                if(L3types.at(i).Contains("L1Muon")) {
                                    hc_Eff_L1SQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                                else {
                                    if( matched_idx > -1 ) {
                                        hc_Eff_L1SQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                    }
                                }
                            }
                            if(matched_L1DQ0) {
                                hc_Eff_L1DQ0.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                                if(L3types.at(i).Contains("L1Muon")) {
                                    hc_Eff_L1DQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                                else {
                                    if( matched_idx > -1 ) {
                                        hc_Eff_L1DQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                    }
                                }
                            }

                            if(matched_L1SQ8) {
                                hc_Eff_L1SQ8.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                                if(L3types.at(i).Contains("L1Muon")) {
                                    hc_Eff_L1SQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                                else {
                                    if( matched_idx > -1 ) {
                                        hc_Eff_L1SQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                    }
                                }
                            }
                            if(matched_L1DQ8) {
                                hc_Eff_L1DQ8.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                                if(L3types.at(i).Contains("L1Muon")) {
                                    hc_Eff_L1DQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                                else {
                                    if( matched_idx > -1 ) {
                                        hc_Eff_L1DQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                    }
                                }
                            }

                            if(matched_L1SQ22) {
                                hc_Eff_L1SQ22.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                                if(L3types.at(i).Contains("L1Muon")) {
                                    hc_Eff_L1SQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                                else {
                                    if( matched_idx > -1 ) {
                                        hc_Eff_L1SQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                    }
                                }
                            }
                            if(matched_L1DQ22) {
                                hc_Eff_L1DQ22.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                                if(L3types.at(i).Contains("L1Muon")) {
                                    hc_Eff_L1DQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                                else {
                                    if( matched_idx > -1 ) {
                                        hc_Eff_L1DQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                    }
                                }
                            }
                        }
                    }

                    // --  Efficiency turn-on / Gen or L1
                    for(unsigned j=0; j<Eff_L3pt_mins.size(); ++j) {
                        hc_EffTO.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                        if( matched_L1SQ0 ) {
                            hc_EffTO_L1SQ0.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if(L3types.at(i).Contains("L1Muon")) {
                                if(genmu.get("l1ptByQ") > Eff_L3pt_mins.at(j)) {
                                    hc_EffTO_L1SQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                            else {
                                if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                    hc_EffTO_L1SQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }
                        if( matched_L1DQ0 ) {
                            hc_EffTO_L1DQ0.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if(L3types.at(i).Contains("L1Muon")) {
                                if(genmu.get("l1ptByQ") > Eff_L3pt_mins.at(j)) {
                                    hc_EffTO_L1DQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                            else {
                                if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                    hc_EffTO_L1DQ0.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }

                        if( matched_L1SQ8 ) {
                            hc_EffTO_L1SQ8.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if(L3types.at(i).Contains("L1Muon")) {
                                if(genmu.get("l1ptByQ") > Eff_L3pt_mins.at(j)) {
                                    hc_EffTO_L1SQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                            else {
                                if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                    hc_EffTO_L1SQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }
                        if( matched_L1DQ8 ) {
                            hc_EffTO_L1DQ8.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if(L3types.at(i).Contains("L1Muon")) {
                                if(genmu.get("l1ptByQ") > Eff_L3pt_mins.at(j)) {
                                    hc_EffTO_L1DQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                            else {
                                if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                    hc_EffTO_L1DQ8.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }

                        if( matched_L1SQ22 ) {
                            hc_EffTO_L1SQ22.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if(L3types.at(i).Contains("L1Muon")) {
                                if(genmu.get("l1ptByQ") > Eff_L3pt_mins.at(j)) {
                                    hc_EffTO_L1SQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                            else {
                                if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                    hc_EffTO_L1SQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }
                        if( matched_L1DQ22 ) {
                            hc_EffTO_L1DQ22.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if(L3types.at(i).Contains("L1Muon")) {
                                if(genmu.get("l1ptByQ") > Eff_L3pt_mins.at(j)) {
                                    hc_EffTO_L1DQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                            else {
                                if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                    hc_EffTO_L1DQ22.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }
                    }

                    // -- res
                    if (matched_idx_res > -1) {
                        if (L3types.at(i).Contains("L1Muon"))
                            continue;

                        if (!L3Coll->at(matched_idx_res).has("charge")) {
                            if (!L3Coll->at(matched_idx_res).has("inner_charge")) {
                                L3Coll->at(matched_idx_res).addVar("charge", L3Coll->at(matched_idx_res).get("inner_charge"));
                            }
                            else {
                                L3Coll->at(matched_idx_res).addVar("charge", 0.);
                            }
                        }

                        double genpt   = genmu.pt;
                        double genqbpt = genmu.get("charge") / genmu.pt;
                        double L3pt    = L3Coll->at(matched_idx_res).pt;
                        double L3qbpt  = L3Coll->at(matched_idx_res).get("charge") / L3Coll->at(matched_idx_res).pt;
                        double res_qbpt = (genqbpt - L3qbpt) / genqbpt;
                        double res_pt   = (L3pt - genpt) / genpt;

                        for(int ipt=0; ipt<n_pt_bins-1; ++ipt) {
                            if( genmu.pt >= pt_bins[ipt] && genmu.pt < pt_bins[ipt+1] ) {
                                vh_L3_qbpt_pt.at(i)[ipt]->Fill( res_qbpt, genWeight );
                                vh_L3_pt_pt.at(i)[ipt]->Fill(   res_pt,   genWeight );
                                break;
                            }
                        }
                        if (genmu.pt > 26.0) {
                            vh_L3_qbpt_pt.at(i)[n_pt_bins-1]->Fill( res_qbpt, genWeight );
                            vh_L3_pt_pt.at(i)[n_pt_bins-1]->Fill(   res_pt,   genWeight );
                        }
                        for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
                            if( genmu.eta >= eta_bins[ieta] && genmu.eta < eta_bins[ieta+1] && genmu.pt > 26.0 ) {
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
            for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                hc_Eff.at(i).at(j)->Save( dir0 );
                hc_Eff_L1SQ0.at(i).at(j)->Save( dir0 );
                hc_Eff_L1DQ0.at(i).at(j)->Save( dir0 );
                hc_Eff_L1SQ8.at(i).at(j)->Save( dir0 );
                hc_Eff_L1DQ8.at(i).at(j)->Save( dir0 );
                hc_Eff_L1SQ22.at(i).at(j)->Save( dir0 );
                hc_Eff_L1DQ22.at(i).at(j)->Save( dir0 );
                // delete hc_Eff.at(i).at(j);
                // delete hc_Eff_L1SQ0.at(i).at(j);
                // delete hc_Eff_L1DQ0.at(i).at(j);
            }

            for(unsigned j=0; j<Eff_L3pt_mins.size(); ++j) {
                hc_EffTO.at(i).at(j)->Save( dir0 );
                hc_EffTO_L1SQ0.at(i).at(j)->Save( dir0 );
                hc_EffTO_L1DQ0.at(i).at(j)->Save( dir0 );
                hc_EffTO_L1SQ8.at(i).at(j)->Save( dir0 );
                hc_EffTO_L1DQ8.at(i).at(j)->Save( dir0 );
                hc_EffTO_L1SQ22.at(i).at(j)->Save( dir0 );
                hc_EffTO_L1DQ22.at(i).at(j)->Save( dir0 );
                // delete hc_EffTO.at(i).at(j);
                // delete hc_EffTO_L1SQ0.at(i).at(j);
                // delete hc_EffTO_L1DQ0.at(i).at(j);
            }
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

    delete f_output;

    printRunTime(timer_total);
}


