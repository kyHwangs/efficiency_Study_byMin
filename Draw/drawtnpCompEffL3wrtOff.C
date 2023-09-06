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
  TString efftag = "IsoMu24", TString ver = "vRun3_01", TString SAMPLE = "Run2023", TString tag = "Muon",
  bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Dir = "./plots_"+ver+"_PPD/"+tag+"/Eff_"+efftag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);

  vector<TString> v_var = {"pt_zoom", "pt", "eta", "phi", "nvtx"};//, "pu"};
  vector< vector<double> > range = {
    {1, 0, 200},  // pt
    {1, 0, 200},  // pt
    {1, -2.4, 2.4},  // eta
    {1, -TMath::Pi(), TMath::Pi()},
    {1, 10, 75}  // PU
  };
  if (tag == "JPsi" || tag == "Bs") {
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

  int n_eta_bins = 23-1;
  double eta_bins[23] = {
    -2.4, -2.1, -1.9, -1.7, -1.6, -1.5, -1.3, -1.2, -0.9,
    -0.3, -0.2,  0.0,  0.2,  0.3,
     0.9,  1.2, 1.3, 1.5, 1.6, 1.7, 1.9, 2.1,  2.4
  };
  vector<TString> etas_str = {"I"};//, "BB", "BE", "EB", "EE"};
  vector<TString> etas_str_long = {"|#eta^{offline}| < 2.4"};//, "|#eta^{offline}| < 0.9", "0.9 < |#eta^{offline}| < 1.2", "1.2 < |#eta^{offline}| < 2.1", "2.1 < |#eta^{offline}| < 2.4"};

  vector<Color_t> v_color = {
    kBlack,
    //kBlue,
    //kRed,
    //kOrange,
    kGreen+2,
    //kCyan+2,
    //kPink+4,
    //kGray+2,
    //kMagenta,
  };
  vector<int> v_marker = {
    20,
    //25,
    //26,
    23,
    //22,
    //26,
    //23,
    //32,
  };
  vector<TString> files = {
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2022G-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-DYToLL_M50_126X-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023B-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023C-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023BC-Eff.root",
    "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023D_DCS-Eff.root",

    //"./Outputs_"+ver+"/hist-"+ver+"-new_align-"+tag+"_Run2023Cv4-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-new_align-"+tag+"_Run2023D-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-old_align-"+tag+"_Run2023Cv4-Eff.root",
    //"./Outputs_"+ver+"/hist-"+ver+"-old_align-"+tag+"_Run2023D-Eff.root",

  };
  vector<TString> types = {
    //TString("Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll").ReplaceAll("my", ""),
    "Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run-1_367660",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run367661_367989",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run367990_368765",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run368766_999999",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run367905_367989",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run368223_368320",
    "Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run367990_368765",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run367905_367989",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_Run368223_368320",

    //"Eff/"+efftag+"/num_Eff_L1SQ22_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/den_Eff_L1SQ22_"+efftag+"_RunAll",
  };
  vector<TString> types_den = {
    //TString("Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll").ReplaceAll("my", ""),
    "Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run-1_367660",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run367661_367989",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run367990_368765",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run368766_999999",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run367905_367989",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run368223_368320",
    "Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run367990_368765",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run367905_367989",
    //"Eff/"+efftag+"/den_Eff_"+efftag+"_Run368223_368320",

    //"Eff/"+efftag+"/den_Eff_L1SQ22_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
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
    //efftag+" : Run2023Cv4 Data, HLTv1.1 & after new align - OLD ALIGN",
    //efftag+" : Run2023D Data - OLD ALIGN",
    //efftag+" : 5 RUNs before new align - OLD ALIGN",
    //efftag+" : 5 RUNs after new align - OLD ALIGN",
  };

  vector<TString> v_pts = {
    "genpt0",
    //"genpt10",
    "genpt26",
    "genpt53",
  };

  vector<TString> v_pts_str = {
    "",
    //"p_{T}^{reco} > 10 GeV",
    "p_{T}^{reco} > 26 GeV",
    "p_{T}^{reco} > 53 GeV",
  };

  for(unsigned i_eta=0; i_eta<etas_str.size(); i_eta++){
    for(int ipt=0; ipt<(int)v_pts.size(); ++ipt) {
      for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

        double xmin = range[ivar][1];
        double xmax = range[ivar][2];
        double ymin = 0.0;
        double ymax = 1.6;

	if(!v_var[ivar].Contains("pt") || v_var[ivar] == "pt_zoom") {
	  ymin = 0.6;//0.5;//0.7;//0.6;//0.85;
	  ymax = 1.2;//1.25;//1.15;//1.2;//1.1;
	}

        TString canvasName = TString::Format("Eff_%s_%s_%s_%s",
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
        //SetLegend( legend, 0.14, 0.71, 0.94, 0.84, -1);
        SetLegend( legend, 0.14, 0.67, 0.94, 0.8, -1);

        bool isFirst = true;
        for(int i = 0; i<(int)files.size(); ++i) {
          TString fileName = files.at(i);

          TString the_type_num = types[i];
          TString the_type_den = types_den[i];
          TString the_type_str = types_str[i].ReplaceAll("my","");

          TString hist_var = v_var[ivar];
          hist_var.ReplaceAll("_zoom", "");

          TString titleX = GetTitleX(hist_var+"_reco");
          TString titleY = "L1+HLT Efficiency";
          if(efftag.Contains("L1sSingleMu22") || efftag == "L1Muon") titleY.ReplaceAll("+HLT", "");

          TString den_name = TString::Format("%s_%s_%s_%s", the_type_den.Data(), etas_str.at(i_eta).Data(), v_pts[ipt].Data(), hist_var.Data());
          TString num_name = TString::Format("%s_%s_%s_%s", the_type_num.Data(), etas_str.at(i_eta).Data(), v_pts[ipt].Data(), hist_var.Data());

          TH1F* den = Get_Hist( fileName, den_name );
          TH1F* num = Get_Hist( fileName, num_name );

          if(v_var[ivar].Contains("pt")) {
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

        legend->Draw();

        TLatex latex;
	Latex_Preliminary_13p6TeV( latex );
        latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
        latex.DrawLatexNDC((i_eta==2?0.66:0.70), 0.89, "#font[42]{#scale[0.8]{"+etas_str_long.at(i_eta)+"}}");
        if(v_var[ivar] != "pt" ) latex.DrawLatexNDC(0.70, 0.84, "#font[42]{#scale[0.8]{"+v_pts_str.at(ipt)+"}}");

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
