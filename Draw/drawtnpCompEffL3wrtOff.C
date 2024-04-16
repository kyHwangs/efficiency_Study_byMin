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

void drawtnpCompEffL3wrtOff(
  //TString efftag = "IsoMu24", bool gen = false, TString ver = "vRun3_05", TString SAMPLE = "Run2023Dv1", TString tag = "Muon0",
  //TString efftag = "IsoMu24", bool gen = false, TString ver = "vRun3_05", TString SAMPLE = "Drell-Yan Simulation", TString tag = "DYToLL_M50_133X",
  //TString efftag = "IsoMu24", bool gen = true, TString ver = "vRun3_0Jpsi", TString SAMPLE = "JPsi Simulation", TString tag = "JPsi_133X",
  TString efftag = "IsoMu24", bool gen = true, TString ver = "vRun3_0Bs", TString SAMPLE = "Bs Simulation", TString tag = "Bs_126X",
  bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Eff = gen? "Effgen" : "Eff";
  TString muon = gen? "gen" : "offline";
  TString run = gen? "" : "_RunAll";
  TString Dir = "./plots_"+ver+"/"+tag+"/"+Eff+"/"+efftag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);

  vector<TString> v_var = {"pt_zoom", "pt", "eta", "phi", "nvtx"};//, "pu"};
  vector< vector<double> > range = {
    {1, 0, 200},  // pt
    {1, 0, 200},  // pt
    {1, -2.4, 2.4},  // eta
    {1, -TMath::Pi(), TMath::Pi()},
    {1, 10, 85}  // PU
  };
  if (tag.Contains( "JPsi") || tag.Contains("Bs")) {
    range.at(0) = {1, 0, 40};
    range.at(1) = {1, 0, 40};
  }

  int n_pt_bins = 29-1;
  double pt_bins[29] = {
    0, 1, 2, 3, 4,
    5, 6, 7, 8, 9,
    10, 12, 15, 18, 20,
    21, 22, 23, 26, 30,
    40, 60, 90, 130, 200,
    300, 450, 700, 1000
  };

  if (efftag.Contains("Mu50")) {
    range = {
      {1, 0, 1000},  // pt
      {1, 0, 1000},  // pt
      {1, -2.4, 2.4},  // eta
      {1, -TMath::Pi(), TMath::Pi()},
      {1, 10, 85}  // PU
    };
  }
  int n_highpt_bins = 21-1;
  double highpt_bins[21] = {
    0, 10, 20, 30, 40,
    45, 46, 47, 48, 49,
    50, 53, 55, 60, 75,
    100, 200, 300, 500, 1000,
    3000
  };

  int n_eta_bins = 23-1;
  double eta_bins[23] = {
    -2.4, -2.1, -1.9, -1.7, -1.6, -1.5, -1.3, -1.2, -0.9,
    -0.3, -0.2,  0.0,  0.2,  0.3,
     0.9,  1.2, 1.3, 1.5, 1.6, 1.7, 1.9, 2.1,  2.4
  };
  vector<TString> etas_str = {"I"};//, "BB", "BE", "EB", "EE"};
  vector<TString> etas_str_long = {"|#eta^{("+muon+")}| < 2.4"};//, "|#eta^{("+muon+")}| < 0.9", "0.9 < |#eta^{("+muon+")}| < 1.2", "1.2 < |#eta^{("+muon+")}| < 2.1", "2.1 < |#eta^{("+muon+")}| < 2.4"};

  vector<Color_t> v_color = {
    kBlack,
    kBlue,
    kRed,
    //kOrange,
    kGreen+2,
    //kCyan+2,
    //kPink+4,
    //kGray+2,
    //kMagenta,
  };
  vector<int> v_marker = {
    22,
    25,
    26,
    23,
    //22,
    //26,
    //23,
    //32,
  };
  vector<TString> files = {
    "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-"+Eff+".root",
    //"./Outputs_"+ver+"/hist-"+ver+"-chaining-"+tag+"-"+Eff+".root",
    "./Outputs_"+ver+"/hist-"+ver+"-noPixOI-"+tag+"-"+Eff+".root",
  };
  vector<TString> types = {
    //TString("Eff/"+efftag+"/num_Eff_"+efftag+run).ReplaceAll("my", ""),
    "Eff/"+efftag+"/num_Eff_"+efftag+run,
    "Eff/"+efftag+"/num_Eff_"+efftag+run,
    "Eff/"+efftag+"/num_Eff_"+efftag+run,
    "Eff/"+efftag+"/num_Eff_"+efftag+run,
    "Eff/"+efftag+"/num_Eff_"+efftag+run,

    //"Eff/"+efftag+"/num_Eff_L1SQ22_"+efftag+run,
    //"Eff/"+efftag+"/den_Eff_L1SQ22_"+efftag+run,
  };
  vector<TString> types_den = {
    //TString("Eff/"+efftag+"/den_Eff_"+efftag+run).ReplaceAll("my", ""),
    "Eff/"+efftag+"/den_Eff_"+efftag+run,
    "Eff/"+efftag+"/den_Eff_"+efftag+run,
    "Eff/"+efftag+"/den_Eff_"+efftag+run,
    "Eff/"+efftag+"/den_Eff_"+efftag+run,
    "Eff/"+efftag+"/den_Eff_"+efftag+run,

    //"Eff/"+efftag+"/den_Eff_L1SQ22_"+efftag+run,
    //"Eff/"+efftag+"/num_Eff_"+efftag+run,
  };
  vector<TString> types_str = {
    "Default Menu",
    //"With IO chaining",
    "With OI change - no Pixetl hits Cut",
  };

  vector<TString> v_pts = {
    "genpt0",
    "genpt10",
    "genpt26",
    "genpt53",
  };

  vector<TString> v_pts_str = {
    "",
    "p_{T}^{("+muon+")} > 10 GeV",
    "p_{T}^{("+muon+")} > 26 GeV",
    "p_{T}^{("+muon+")} > 53 GeV",
  };

  for(unsigned i_eta=0; i_eta<etas_str.size(); i_eta++){
    for(int ipt=0; ipt<(int)v_pts.size(); ++ipt) {
      for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

        double xmin = range[ivar][1];
        double xmax = range[ivar][2];
        double ymin = 0.0;
        double ymax = 1.6;

        if(!v_var[ivar].Contains("pt") || v_var[ivar] == "pt_zoom") {
          ymin = 0.6;//0.7;//0.5;//0.6;//0.85;
          ymax = 1.2;//1.15;//1.25;//1.2;//1.1;
        }

        TString canvasName = TString::Format("%s_%s_%s_%s_%s",
                                             Eff.Data(),
                                             efftag.Data(),
                                             etas_str.at(i_eta).Data(),
                                             v_pts[ipt].Data(),
                                             v_var[ivar].Data());
        canvasName.ReplaceAll(".","p").ReplaceAll("-","_").ReplaceAll("my", "");
        TCanvas *c;
        SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
        c->cd();
        if(isLogy) c->SetLogy();
        if(tag == "Zprime" && v_var[ivar].Contains("pt")) c->SetLogx();

        TLegend *legend;
        SetLegend( legend, 0.17, 0.72, 0.94, 0.85, -1);

        bool isFirst = true;
        for(int i = 0; i<(int)files.size(); ++i) {
          TString fileName = files.at(i);

          TString the_type_num = types[i];
          TString the_type_den = types_den[i];
          TString the_type_str = types_str[i].ReplaceAll("my","");

          TString hist_var = v_var[ivar];
          hist_var.ReplaceAll("_zoom", "");

          TString titleX = GetTitleX(hist_var+"_reco");
          if(gen) titleX = GetTitleX(hist_var+"_gen");
          TString titleY = "L1+HLT Efficiency";
          if(efftag.Contains("L1sSingleMu22") || efftag == "L1Muon") titleY.ReplaceAll("+HLT", "");

          TString den_name = TString::Format("%s_%s_%s_%s", the_type_den.Data(), etas_str.at(i_eta).Data(), v_pts[ipt].Data(), hist_var.Data());
          TString num_name = TString::Format("%s_%s_%s_%s", the_type_num.Data(), etas_str.at(i_eta).Data(), v_pts[ipt].Data(), hist_var.Data());

          TH1F* den = Get_Hist( fileName, den_name );
          TH1F* num = Get_Hist( fileName, num_name );

          if(v_var[ivar].Contains("pt")) {
            if (efftag.Contains("Mu50")) {
              den = (TH1F*)den->Rebin(n_highpt_bins, den_name+"_rb", highpt_bins);
              num = (TH1F*)num->Rebin(n_highpt_bins, num_name+"_rb", highpt_bins);
            } else {
              den = (TH1F*)den->Rebin(n_pt_bins, den_name+"_rb", pt_bins);
              num = (TH1F*)num->Rebin(n_pt_bins, num_name+"_rb", pt_bins);
            }
          }
          else if(v_var[ivar] == "eta") {
            den = (TH1F*)den->Rebin(n_eta_bins, den_name+"_rb", eta_bins);
            num = (TH1F*)num->Rebin(n_eta_bins, num_name+"_rb", eta_bins);
          }
          else{
            den = (TH1F*)den->Rebin(range[ivar][0]);
            num = (TH1F*)num->Rebin(range[ivar][0]);
          }

          int nbins = den->GetNbinsX();
          TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
          g->Divide(num, den, "n e0");
          //g->Divide(num, den, "pois");

          for(int ip=0; ip<nbins; ++ip) {
            if(g->GetPointY(ip) == 0.)  g->SetPointEYhigh(ip, 0.0);
          }

          g->SetTitle("");
          g->SetMarkerSize(1.5);
          g->SetMarkerStyle(v_marker[i]);
          g->SetMarkerColor(v_color[i]);
          g->SetLineColor(  v_color[i]);
          g->SetLineWidth(1);

          g->GetXaxis()->SetLimits( xmin, xmax );
          g->GetXaxis()->SetRangeUser( xmin, xmax );
          g->GetYaxis()->SetRangeUser( ymin, ymax );

          SetAxis_SinglePad( g->GetXaxis(), g->GetYaxis(), titleX, titleY );

          if(isFirst) {
            g->Draw("APE");
            isFirst = false;
          }
          else {
            g->Draw("PE same");
          }

          legend->AddEntry( g, TString::Format("%s", the_type_str.Data()), "lep" );
        }

        TLine eff1p0(xmin,1.0, xmax,1.0);
        eff1p0.SetLineColor(kGray);
        eff1p0.SetLineWidth(1);
        eff1p0.Draw("same");
        legend->Draw();

        TString L3str = "";
        if(efftag == "L2Muon") L3str = "L2 Muon";
        else if(efftag == "hltOI") L3str = "Outside-in L3 MuonTrack";
        else if(efftag == "hltPixelTracksInRegionL2") L3str = "PixelTrack near L2";
        else if(efftag == "hltPixelTracksInRegionL1") L3str = "PixelTrack near L1";
        else if(efftag == "hltPixelTracks") L3str = "PixelTrack";
        else if(efftag == "hltIter0FromL1") L3str = "Inside-out L3 MuonTrack from L1";
        else if(efftag == "hltL3FromL2Merged") L3str = "L3 MuonTrack from L2";
        else if(efftag == "hltL3Merged") L3str = "L3 MuonTrack";
        else if(efftag.Contains("hltIterL3MuonNoID")) L3str = "L3 Muon";
        else if(efftag == "hltIterL3Muon") L3str = "L3 Muon after Trigger ID";
        else if(efftag.Contains("L1sSingleMu22")) L3str = "Good quality L1 muon with p_{T}^{L1} > 22 GeV";
        else if(efftag.Contains("IsoMu24")) L3str = "Isolated muon with p_{T} > 24 GeV";
        else if(efftag.Contains("Mu24")) L3str = "Non-isolated muon with p_{T} > 24 GeV";
        else if(efftag.Contains("Mu50OrOldMu100OrTkMu100")) L3str = "Non-isolated muon with p_{T} > 50 GeV";
        else if(efftag.Contains("Mu50")) L3str = "Non-isolated muon with p_{T} > 50 GeV";

        TLatex latex;
        Latex_Preliminary_13p6TeV( latex );
        latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
        latex.DrawLatexNDC(0.16, 0.90, "#font[42]{#scale[0.6]{"+L3str+"}}");
        latex.DrawLatexNDC((i_eta==2?0.66:0.70), 0.89, "#font[42]{#scale[0.8]{"+etas_str_long.at(i_eta)+"}}");
        if(v_var[ivar] != "pt" ) latex.DrawLatexNDC(0.68, 0.84, "#font[42]{#scale[0.8]{"+v_pts_str.at(ipt)+"}}");

        TString logy_tag = isLogy ? "_log" : "";
        // CMS_lumi(c, 98, 11);
        c->Modified();  c->Update();  c->RedrawAxis();
        gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
        //c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
        c->SaveAs(Dir+canvasName+logy_tag+".png","png");
        gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

        c->Close();
      }
    }
  }
  printRunTime(timer_total);
}
