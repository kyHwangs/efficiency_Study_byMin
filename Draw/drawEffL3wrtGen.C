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

// echo 'gROOT->LoadMacro("drawEffL3wrtGen.C+"); gSystem->Exit(0);' | root -b -l
// rootbq 'drawEffL3wrtGen.C("v01", "DY Winter21", "DYToLL_M50", "L1Tk")'
// rootbq 'drawEffL3wrtGen.C("v01", "DY Winter21", "DYToLL_M50", "")'

void drawEffL3wrtGen(
  TString ver = "v02", TString SAMPLE = "DY", TString tag = "DYToLL_M50",
  TString L1tag = "L1SQ0", TString L1str = "L1: SQ, no p_{T} cut",
  bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);


  TString efftag = "L3_Gen";
  TString Dir = "./plots_"+ver+"/"+tag+"/Eff_"+efftag+"/"+L1tag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);

  vector<Color_t> v_color = {
    kWhite,
    kWhite,
    kGray+2,
    kGray+2,
    kBlack,
    kBlack,
    kBlue,
    kBlue,
    kRed,
    kRed,

    kMagenta,
    kGreen+2,
    kYellow
  };

  vector<int> v_marker = {
    1,
    1,
    25,
    21,
    24,
    20,
    26,
    22,
    32,
    23,

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
    {1, 30, 81}  // PU
  };
  if (tag == "Jpsi") {
    range.at(0) = {1, 0, 60};
    range.at(1) = {1, 0, 60};
  }
  if (tag == "MuGunPU") {
    range.at(0) = {1, 1., 1000};
    range.at(1) = {1, 1., 1000};
  }
  if (tag == "Zprime_M6000") {
    range.at(0) = {1, 30, 1000};
    range.at(1) = {1, 30, 1000};
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

  int n_eta_bins = 23-1;
  double eta_bins[23] = {
    -2.4, -2.1, -1.9, -1.7, -1.6, -1.5, -1.3, -1.2, -0.9,
    -0.3, -0.2,  0.0,  0.2,  0.3,
     0.9,  1.2, 1.3, 1.5, 1.6, 1.7, 1.9, 2.1,  2.4
  };

  TString file_tag = (ver == "v33") ? "Sim" : "";
  vector<TString> files = {
    "./Outputs_"+ver+"/hist-"+ver+"-wp00-"+tag+"-Eff"+file_tag+".root",
    "./Outputs_"+ver+"/hist-"+ver+"-wp01L2L1-"+tag+"-Eff"+file_tag+".root",

    "./Outputs_"+ver+"/hist-"+ver+"-wp00-"+tag+"-Eff"+file_tag+".root",
    "./Outputs_"+ver+"/hist-"+ver+"-wp01L2L1-"+tag+"-Eff"+file_tag+".root",

    "./Outputs_"+ver+"/hist-"+ver+"-wp00-"+tag+"-Eff"+file_tag+".root",
    "./Outputs_"+ver+"/hist-"+ver+"-wp01L2L1-"+tag+"-Eff"+file_tag+".root",

    "./Outputs_"+ver+"/hist-"+ver+"-wp00-"+tag+"-Eff"+file_tag+".root",
    "./Outputs_"+ver+"/hist-"+ver+"-wp01L2L1-"+tag+"-Eff"+file_tag+".root",

    "./Outputs_"+ver+"/hist-"+ver+"-wp00-"+tag+"-Eff"+file_tag+".root",
    "./Outputs_"+ver+"/hist-"+ver+"-wp01L2L1-"+tag+"-Eff"+file_tag+".root",
  };

  vector<TString> types = {
    "Eff/den_Eff_"+L1tag+"_hltOI",  // L1
    "Eff/den_Eff_"+L1tag+"_hltOI",  // L1

    "Eff/den_Eff_"+L1tag+"_hltOI",  // L1
    "Eff/den_Eff_"+L1tag+"_hltOI",  // L1

    "Eff/num_Eff_"+L1tag+"_hltOI",
    "Eff/num_Eff_"+L1tag+"_hltOI",

    "Eff/num_Eff_"+L1tag+"_hltL3FromL2Merged",
    "Eff/num_Eff_"+L1tag+"_hltL3FromL2Merged",

    "Eff/num_Eff_"+L1tag+"_hltL3Merged",
    "Eff/num_Eff_"+L1tag+"_hltL3Merged",
  };

  vector<TString> types_den = {
    "Eff/den_Eff_L1Muon",
    "Eff/den_Eff_L1Muon",

    "Eff/den_Eff_L1Muon",
    "Eff/den_Eff_L1Muon",

    "Eff/den_Eff_hltOI",
    "Eff/den_Eff_hltOI",

    "Eff/den_Eff_hltL3FromL2Merged",
    "Eff/den_Eff_hltL3FromL2Merged",

    "Eff/den_Eff_hltL3Merged",
    "Eff/den_Eff_hltL3Merged",
  };

  vector<TString> types_str = {
    "No IO MVA cut",
    "IO MVA_{L2, L1} > 0.01",

    "L1 muons: "+L1str,
    "L1 muons: "+L1str,

    "OI",
    "OI",

    "OI + IO (L2)",
    "OI + IO (L2)",

    "OI + IO (L2) + IO (L1)",
    "OI + IO (L2) + IO (L1)",
  };

  vector<TString> v_pts = {
    "genpt0",
    "genpt10",
    "genpt26"
  };

  vector<TString> v_pts_str = {
    "",
    "p_{T}^{gen} > 10 GeV",
    "p_{T}^{gen} > 26 GeV"
  };

  for(int ipt=0; ipt<(int)v_pts.size(); ++ipt) {

    for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

      double xmin = range[ivar][1];
      double xmax = range[ivar][2];
      double ymin = 0.0;
      double ymax = 1.6;

      if(!v_var[ivar].Contains("pt") || v_var[ivar] == "pt_zoom") {
        ymin = 0.7;
        ymax = 1.2;
      }

      TString canvasName = TString::Format("Eff_%s_%s_%s_%s_%s",
                                           efftag.Data(),
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
        TString fileName = files.at(i);

        TString the_type_num = types[i];
        TString the_type_den = types_den[i];
        TString the_type_str = types_str[i];

        TString hist_var = v_var[ivar];
        hist_var.ReplaceAll("_zoom", "");

        TString titleX = GetTitleX(hist_var+"_gen");
        TString titleY = "L1 + L3 efficiency";

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
      Latex_Simulation_14TeV( latex );
      latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
      latex.DrawLatexNDC(0.20, 0.87, "#font[42]{#scale[0.8]{"+L1str+"}}");
      if(v_var[ivar] != "pt" )
        latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{"+v_pts_str.at(ipt)+"}}");


      TString logy_tag = isLogy ? "_log" : "";
      // CMS_lumi(c, 98, 11);
      c->Modified();  c->Update();  c->RedrawAxis();
      gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
      c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
      // c->SaveAs(Dir+canvasName+logy_tag+".C","C");
      // c->SaveAs(Dir+canvasName+logy_tag+".root","root");
      gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

      c->Close();
    }
  }


  printRunTime(timer_total);
}
