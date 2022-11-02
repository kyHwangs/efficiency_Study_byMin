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

void drawtnpEffL3wrtL1(
  TString ver = "vRun3_01", TString SAMPLE = "Run 2022D", TString tag = "Muon_Run2022D",
  TString L1tag = "L1SQ0", TString L1str = "L1: SQ, no p_{T} cut",
  bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString efftag = "L3_L1";
  TString Dir = "./plots_"+ver+"/"+tag+"/Eff_"+efftag+"/"+L1tag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);

  vector<Color_t> v_color = {
    kWhite,
    kWhite,
    kBlack,
    kBlack,
    kBlue,
    kBlue,
    kRed,
    kRed,
    kMagenta,
    kMagenta,
    kGreen+2,
    kGreen+2,

    kYellow
  };

  vector<int> v_marker = {
    1,
    1,
    24,
    20,
    26,
    22,
    32,
    23,
    27,
    33,
    28,
    34,

    21,
    24,
    25
  };


  vector<TString> v_var = {"pt_zoom", "pt", "eta", "phi", "pu"};
  vector< vector<double> > range = {
    {1, 0, 200},  // pt
    {1, 0, 200},  // pt
    {1, -2.4, 2.4},  // eta
    {1, -TMath::Pi(), TMath::Pi()},
    {1, 10, 60}  // PU
  };

  int n_pt_bins = 29-1;
  double pt_bins[29] = {
    0, 1, 2, 3, 4,
    5, 6, 7, 8, 9,
    10, 12, 15, 18, 20,
    21, 22, 23, 26, 30,
    40, 60, 90, 130, 200,
    300, 450, 700, 1000
  };

  int n_eta_bins = 23-1;
  double eta_bins[23] = {
    -2.4, -2.1, -1.9, -1.7, -1.6, -1.5, -1.3, -1.2, -0.9,
    -0.3, -0.2,  0.0,  0.2,  0.3,
     0.9,  1.2, 1.3, 1.5, 1.6, 1.7, 1.9, 2.1,  2.4
  };

  // vector<TString> files = {
  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",
  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",

  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",
  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",

  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",
  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",

  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",
  //   "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root",
  // };

  TString file_name = "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"-Eff.root";

  vector<TString> types = {
    "Eff/hltOI/num_Eff_"+L1tag+"_hltOI_RunAll_I",
    "Eff/hltOI/num_Eff_"+L1tag+"_hltOI_Run357610_357900_I",

    "Eff/hltOI/num_Eff_"+L1tag+"_hltOI_RunAll_I",
    "Eff/hltOI/num_Eff_"+L1tag+"_hltOI_Run357610_357900_I",

    "Eff/hltL3FromL2Merged/num_Eff_"+L1tag+"_hltL3FromL2Merged_RunAll_I",
    "Eff/hltL3FromL2Merged/num_Eff_"+L1tag+"_hltL3FromL2Merged_Run357610_357900_I",

    "Eff/hltIterL3MuonNoID/num_Eff_"+L1tag+"_hltIterL3MuonNoID_RunAll_I",
    "Eff/hltIterL3MuonNoID/num_Eff_"+L1tag+"_hltIterL3MuonNoID_Run357610_357900_I",

    "Eff/hltIterL3Muon/num_Eff_"+L1tag+"_hltIterL3Muon_RunAll_I",
    "Eff/hltIterL3Muon/num_Eff_"+L1tag+"_hltIterL3Muon_Run357610_357900_I",

    //"Eff/IsoMu24/num_Eff_"+L1tag+"_IsoMu24_RunAll_I",
    //"Eff/IsoMu24/num_Eff_"+L1tag+"_IsoMu24_Run357610_357900_I",
  };

  vector<TString> types_den = {
    "Eff/hltOI/den_Eff_"+L1tag+"_hltOI_RunAll_I",
    "Eff/hltOI/den_Eff_"+L1tag+"_hltOI_Run357610_357900_I",

    "Eff/hltOI/den_Eff_"+L1tag+"_hltOI_RunAll_I",
    "Eff/hltOI/den_Eff_"+L1tag+"_hltOI_Run357610_357900_I",

    "Eff/hltL3FromL2Merged/den_Eff_"+L1tag+"_hltL3FromL2Merged_RunAll_I",
    "Eff/hltL3FromL2Merged/den_Eff_"+L1tag+"_hltL3FromL2Merged_Run357610_357900_I",

    "Eff/hltIterL3MuonNoID/den_Eff_"+L1tag+"_hltIterL3MuonNoID_RunAll_I",
    "Eff/hltIterL3MuonNoID/den_Eff_"+L1tag+"_hltIterL3MuonNoID_Run357610_357900_I",

    "Eff/hltIterL3Muon/den_Eff_"+L1tag+"_hltIterL3Muon_RunAll_I",
    "Eff/hltIterL3Muon/den_Eff_"+L1tag+"_hltIterL3Muon_Run357610_357900_I",

    //"Eff/IsoMu24/den_Eff_"+L1tag+"_IsoMu24_RunAll_I",
    //"Eff/IsoMu24/den_Eff_"+L1tag+"_IsoMu24_Run357610_357900_I",
  };

  vector<TString> types_str = {
    "Run < 357610",
    "Run #geq 357610",

    "OI",
    "OI",

    "OI + IO (L2)",
    "OI + IO (L2)",

    "L3 muon",
    "L3 muon",

    "L3 muon + ID",
    "L3 muon + ID",

    //"IsoMu24",
    //"IsoMu24",
  };

  vector<TString> v_pts = {
    "genpt0",
    "genpt26"
  };

  vector<TString> v_pts_str = {
    "",
    "p_{T}^{reco} > 26 GeV"
  };

  for(int ipt=0; ipt<(int)v_pts.size(); ++ipt) {

    for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

      double xmin = range[ivar][1];
      double xmax = range[ivar][2];
      double ymin = 0.0;
      double ymax = 1.6;

      if(!v_var[ivar].Contains("pt") || v_var[ivar] == "pt_zoom") {
        ymin = 0.85;
        ymax = 1.1;
      }

      TString canvasName = TString::Format("Eff_%s_%s_%s_%s_%s",
                                           efftag.Data(),
                                           // type_name.Data(),
                                           tag.Data(),
                                           L1tag.Data(),
                                           v_pts[ipt].Data(),
                                           v_var[ivar].Data());
      canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
      TCanvas *c;
      SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
      c->cd();
      if(isLogy)
        c->SetLogy();
      if(tag == "MuGunPU" && v_var[ivar].Contains("pt"))
        c->SetLogx();

      TLegend *legend;
      SetLegend( legend, 0.15, 0.65, 0.90, 0.83, -1);
      legend->SetNColumns(2);

      bool isFirst = true;
      for(int i = 0; i<(int)types.size(); ++i) {
        TString fileName = file_name;  // files.at(i);

        TString the_type_num = types[i];
        TString the_type_den = types_den[i];
        TString the_type_str = types_str[i];

        TString hist_var = v_var[ivar];
        hist_var.ReplaceAll("_zoom", "");

        TString titleX = GetTitleX(hist_var+"_reco");
        TString titleY = "L3/L1 efficiency";

        TString den_name = TString::Format("%s_%s_%s", the_type_den.Data(),
                                                       v_pts[ipt].Data(),
                                                       hist_var.Data());
        TString num_name = TString::Format("%s_%s_%s", the_type_num.Data(),
                                                       v_pts[ipt].Data(),
                                                       hist_var.Data());

        TH1F* den = Get_Hist( fileName, den_name );
        TH1F* num = Get_Hist( fileName, num_name );

        if(v_var[ivar] == "pt" || v_var[ivar] == "pt_zoom") {
          den = (TH1F*)den->Rebin(n_pt_bins, den_name+"_rb", pt_bins);
          num = (TH1F*)num->Rebin(n_pt_bins, num_name+"_rb", pt_bins);
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

        c->cd();

        TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
        g->Divide(num, den, "n e0");

        for(int ip=0; ip<nbins; ++ip) {
          if(g->GetPointY(ip) == 0.)  g->SetPointEYhigh(ip, 0.0);
        }

        g->SetTitle("");
        g->SetMarkerSize(1.5);
        g->SetMarkerStyle(v_marker[i]);
        g->SetMarkerColor(v_color[i]);
        g->SetLineColor(  v_color[i]);
        g->SetLineWidth(1);
        if(i%2==0) {
          g->SetLineStyle(2);
        }

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

      legend->Draw();

      TLatex latex;
      Latex_Preliminary_13p6TeV( latex );
      latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
      latex.DrawLatexNDC(0.20, 0.87, "#font[42]{#scale[0.8]{"+L1str+"}}");
      if(v_var[ivar] != "pt" )
        latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{"+v_pts_str.at(ipt)+"}}");


      TString logy_tag = isLogy ? "_log" : "";
      c->Modified();  c->Update();  c->RedrawAxis();
      gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
      c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
      gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

      c->Close();
    }
  }


  printRunTime(timer_total);
}
