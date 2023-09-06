
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TSystem.h>

#define DEBUG (0)
#include "PlotTools.h"

using namespace std;

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

void drawtnp2DEffL3wrtL1(
  TString efftag = "hltIterL3Muon", TString ver = "vRun3_01", TString SAMPLE = "Run2023", TString tag = "Muon",
  TString L1tag = "L1SQ22", TString L1str = "L1 qual > 11, p_{T}^{L1} > 22 GeV",
  //TString L1tag = "L1DQ8", TString L1str = "L1 qual > 7, p_{T}^{L1} > 8 GeV",
  bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Dir = "./plots_"+ver+"_PPD/"+tag+"/Eff_"+efftag+"/"+L1tag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);

  vector<TString> files = {
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023D_DCS-Eff.root",
  };
  vector<TString> types = {
    //"Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_Run-1_367660",
    //"Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_Run367661_367989",
    //"Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_Run367990_368765",
    //"Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_Run368766_999999",
    //"Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_Run367905_367989",
    //"Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_Run368223_368320",
    "Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_RunAll",
    "Eff/"+efftag+"/num_2DEff_"+L1tag+"_"+efftag+"_RunAll",
  };
  vector<TString> types_den = {
    //TString("Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_RunAll").ReplaceAll("my", ""),

    //"Eff/"+efftag+"den_2DEff_"+L1tag+"_"+efftag+"_Run-1_367660",
    //"Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_Run367661_367989",
    //"Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_Run367990_368765",
    //"Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_Run368766_999999",
    //"Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_Run367905_367989",
    //"Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_Run368223_368320",
    "Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_RunAll",
    "Eff/"+efftag+"/den_2DEff_"+L1tag+"_"+efftag+"_RunAll",

  };
  vector<TString> types_str = {
    //efftag+" : Run2022G Data",
    //efftag+" : Run3Winer23 DY",
    efftag+" : Run2023BC Data",
    //efftag+" : Run2023 Data, HLTv1.0",
    //efftag+" : Run2023 Data, HLTv1.1 & before new align",
    //efftag+" : Run2023 Data, HLTv1.1 & after new align",
    //efftag+" : Run2023 Data, and so on",
    //efftag+" : 5 RUNs before new align",
    //efftag+" : 5 RUNs after new align",
    efftag+" : Run2023D Data (DCSOnly json)",
  };

  TH2F* prev_eff;
  for(int i = 0; i<(int)files.size(); ++i) {
    TString canvasName = TString::Format("2DEff_%s_%s_I_genpt26_eta_phi",
                                         types_str.at(i).Data(),
                                         L1tag.Data());
    canvasName.ReplaceAll(".","p").ReplaceAll("-","_").ReplaceAll("my", "").ReplaceAll(" : ","_").ReplaceAll(" ","_");
    TCanvas *c = new TCanvas("eff", "eff", 1500, 900);
    c->cd();

    TString fileName = files.at(i);

    TString the_type_num = types[i];
    TString the_type_den = types_den[i];

    TString den_name = TString::Format("%s_I_genpt26_eta_phi", the_type_den.Data());
    TString num_name = TString::Format("%s_I_genpt26_eta_phi", the_type_num.Data());

    TH2F* den = Get_Hist_2D( fileName, den_name );
    TH2F* num = Get_Hist_2D( fileName, num_name );

    TH2F* eff = (TH2F*)num->Clone();
    eff->Divide(den);
    //eff->Draw("colz text");
    eff->Draw("colz");
    eff->SetTitle(types_str.at(i).Data());
    eff->GetXaxis()->SetTitle(GetTitleX("eta_reco"));
    eff->GetYaxis()->SetTitle(GetTitleX("phi_reco"));
    eff->SetMaximum(1.0); eff->SetMinimum(0.0);

    c->Modified();  c->Update();  c->RedrawAxis();
    gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
    //c->SaveAs(Dir+canvasName+".pdf","pdf");
    c->SaveAs(Dir+canvasName+".png","png");
    gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

    if(i==0) prev_eff = eff;
    else {
      TCanvas *c2 = new TCanvas("ratio_eff", "ratio_eff", 1500, 900);
      c2->cd();

      TH2F* tmp_prev = (TH2F*)eff->Clone();
      eff->Divide(prev_eff);
      eff->Draw("colz");
      eff->SetTitle("Efficiency Ratio : "+types_str.at(i));
      eff->SetMaximum(1.05); eff->SetMinimum(0.85);
      c2->SaveAs(Dir+canvasName+"_ratio.png","png");
      c2->Close();

      prev_eff = tmp_prev;
    }
    c->Close();
  }
  printRunTime(timer_total);
}
