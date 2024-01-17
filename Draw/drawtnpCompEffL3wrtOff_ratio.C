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

void canvas_margin(TCanvas *c, TPad *c_up, TPad *c_down);
TGraphAsymmErrors* MakeRatioGraph(TGraphAsymmErrors *_g_Type1, TGraphAsymmErrors *_g_Type2, TString errorPropagation = "AoverB");
void SetAxis_Pad_Top( TAxis *X_axis, TAxis *Y_axis, TString YTitle );
void SetAxis_Pad_Bottom( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle );
Double_t ReturnLargerValue(Double_t a, Double_t b)
{
  if( a > b )
    return a;
  else
    return b;
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

void drawtnpCompEffL3wrtOff_ratio(
  TString efftag = "IsoMu24", TString ver = "vRun3_04", TString SAMPLE = "Run2022, 2023", TString tag = "Muon",
  bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Dir = "./plots_"+ver+"/"+tag+"/Eff_"+efftag+"/";
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

  if (efftag.Contains("Mu50")) {
    range = {
      {1, 0, 500},  // pt
      {1, 0, 500},  // pt
      {1, -2.4, 2.4},  // eta
      {1, -TMath::Pi(), TMath::Pi()},
      {1, 10, 75}  // PU
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
  vector<TString> etas_str_long = {"|#eta^{offline}| < 2.4"};//, "|#eta^{offline}| < 0.9", "0.9 < |#eta^{offline}| < 1.2", "1.2 < |#eta^{offline}| < 2.1", "2.1 < |#eta^{offline}| < 2.4"};

  vector<Color_t> v_color = {
    //kBlue,
    kRed,
    kBlack,
    //kOrange,
    //kGreen+2,
    //kCyan+2,
    //kPink+4,
    //kGray+2,
    //kMagenta,
  };
  vector<int> v_marker = {
    //21,
    23,
    20,
    //22,
    //25,
    //26,
    //23,
    //22,
    //26,
    //23,
    //32,
  };
  vector<TString> files = {
    //"./Outputs_"+ver+"/hist-"+ver+"-DYToLL_M50_126X-Eff_1326.root",
    "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2022-Eff_1326.root",
    "./Outputs_"+ver+"/hist-"+ver+"-"+tag+"_Run2023-Eff_1326.root",
  };
  vector<TString> types = {
    //TString("Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll").ReplaceAll("my", ""),
    "Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",

    //"Eff/"+efftag+"/num_Eff_L1SQ22_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/den_Eff_L1SQ22_"+efftag+"_RunAll",
  };
  vector<TString> types_den = {
    //TString("Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll").ReplaceAll("my", ""),
    "Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",
    "Eff/"+efftag+"/den_Eff_"+efftag+"_RunAll",

    //"Eff/"+efftag+"/den_Eff_L1SQ22_"+efftag+"_RunAll",
    //"Eff/"+efftag+"/num_Eff_"+efftag+"_RunAll",
  };
  vector<TString> types_str = {
    //"Drell-Yan Simulation",
    "Run2022 Data (35 fb^{-1} )",
    "Run2023 Data (27 fb^{-1} )",
  };

  vector<TString> v_pts = {
    "genpt0",
    //"genpt10",
    "genpt26",
    "genpt53",
  };

  vector<TString> v_pts_str = {
    "",
    //"p_{T}^{offline} > 10 GeV",
    "p_{T}^{offline} > 26 GeV",
    "p_{T}^{offline} > 53 GeV",
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

        TString canvasName = TString::Format("Eff_%s_%s_%s_%s",
                                             efftag.Data(),
                                             etas_str.at(i_eta).Data(),
                                             v_pts[ipt].Data(),
                                             v_var[ivar].Data());
        canvasName.ReplaceAll(".","p").ReplaceAll("-","_").ReplaceAll("my", "");
        TCanvas *c = new TCanvas("c_eff", "", 900, 900);
        TPad *c_up = new TPad("c_up", "", 0.01,0.25, 0.99,0.99);
        TPad *c_down = new TPad("c_down", "", 0.01,0.01, 0.99,0.25);
        canvas_margin(c, c_up, c_down);

        c_up->cd();
        if(isLogy) c_up->SetLogy();
        if(tag == "Zprime" && v_var[ivar].Contains("pt")) c_up->SetLogx();

        TLegend *legend;
        SetLegend( legend, 0.15, 0.685, 0.65, 0.775, -1);

        bool isFirst = true;
        TGraphAsymmErrors* g_ref;
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

          SetAxis_Pad_Top( g->GetXaxis(), g->GetYaxis(), titleY );
          if(isFirst) {
            g->Draw("APE");
            g_ref = g;
            isFirst = false;
          }
          else {
            g->Draw("PE same");
            c_down->cd();
            TGraphAsymmErrors* g_ratio = MakeRatioGraph(g, g_ref);
            g_ratio->GetXaxis()->SetLimits( xmin, xmax );
            g_ratio->GetXaxis()->SetRangeUser( xmin, xmax );
            g_ratio->GetYaxis()->SetRangeUser( 0.89, 1.11 );
            if(i == 1) g_ratio->Draw("APE");
            else g_ratio->Draw("PE same");
            SetAxis_Pad_Bottom( g_ratio->GetXaxis(), g_ratio->GetYaxis(), titleX, titleY );
          }

          c_up->cd();
          legend->AddEntry( g, TString::Format("%s", the_type_str.Data()), "lep" );
        }

        TLine eff1p0(xmin,1.0, xmax,1.0);
        eff1p0.SetLineColor(kGray);
        eff1p0.SetLineWidth(1);
        eff1p0.Draw("same");
        legend->Draw();

	c->cd();
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
        //latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
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

void canvas_margin(TCanvas *c, TPad *c_up, TPad *c_down){
  c_up->SetTopMargin( 0.06 );
  c_up->SetBottomMargin( 0.02 );
  c_up->SetLeftMargin( 0.11 );
  c_up->SetRightMargin( 0.02 );
  c_up->Draw();

  c_down->SetTopMargin( 0.01 );
  c_down->SetBottomMargin( 0.26 );
  c_down->SetLeftMargin( 0.11 );
  c_down->SetRightMargin( 0.02 );
  c_down->SetGridx();
  c_down->SetGridy();
  c_down->Draw();

  c->SetTopMargin( 0.05 );   //0.05
  c->SetBottomMargin( 0.13 );//0.13
  c->SetRightMargin( 0.05 ); //0.05
  c->SetLeftMargin( 0.16 );  //0.16

  gStyle->SetOptStat(0);
}

void SetAxis_Pad_Top( TAxis *X_axis, TAxis *Y_axis, TString YTitle ){
  X_axis->SetLabelSize(0.000);
  X_axis->SetTitleSize(0.000);

  Y_axis->SetTitle( YTitle );
  Y_axis->SetTitleSize(0.05);
  Y_axis->SetTitleOffset(1.0);
  Y_axis->SetLabelSize(0.05);
}
void SetAxis_Pad_Bottom( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle ){
  X_axis->SetTitle( XTitle );
  X_axis->SetLabelSize(0.1);
  X_axis->SetTitleOffset(1.0);
  X_axis->SetTitleSize(0.12);
  X_axis->SetNoExponent();
  X_axis->SetMoreLogLabels();

  Y_axis->SetTitle( "2023/2022" );
  Y_axis->SetTitleSize(0.1);
  Y_axis->SetTitleOffset(0.5);
  Y_axis->SetLabelSize(0.1);
}

TGraphAsymmErrors* MakeRatioGraph(TGraphAsymmErrors *_g_Type1, TGraphAsymmErrors *_g_Type2, TString errorPropagation = "AoverB"){
  TGraphAsymmErrors* g_Type1 = (TGraphAsymmErrors*)_g_Type1->Clone();
  TGraphAsymmErrors* g_Type2 = (TGraphAsymmErrors*)_g_Type2->Clone();

  TGraphAsymmErrors*  g_ratio = (TGraphAsymmErrors*)g_Type1->Clone();
  g_ratio->Set(0); // --Remove all points (reset) -- //
  
  Int_t nPoint = g_Type1->GetN();
  Int_t nPoint_2 = g_Type2->GetN();
  if( nPoint != nPoint_2 ) {
    printf("# points is different bewteen two graph... return NULL\n");
    cout << "\tnPoint Num = " << nPoint << endl;
    cout << "\tnPoint Den = " << nPoint_2 << endl;
    return NULL;
  }

  for(Int_t i_p=0; i_p<nPoint; i_p++)
    {
      Double_t x_Type1, y_Type1;
      g_Type1->GetPoint(i_p, x_Type1, y_Type1);
      Double_t error_Type1 = ReturnLargerValue( g_Type1->GetErrorYhigh(i_p), g_Type1->GetErrorYlow(i_p) );

      //Get Type2 point
      Double_t x_Type2, y_Type2;
      g_Type2->GetPoint(i_p, x_Type2, y_Type2);
      Double_t error_Type2 = ReturnLargerValue( g_Type2->GetErrorYhigh(i_p), g_Type2->GetErrorYlow(i_p) );

      Double_t ratio;
      Double_t ratio_error = 999.;
      if(y_Type2 != 0)
        {
          ratio = y_Type1 / y_Type2;
          if(errorPropagation == "AoverB")
            ratio_error = GetUncorrelatedError(y_Type1, error_Type1, y_Type2, error_Type2);
        }
        else
          {
            ratio = 0;
            ratio_error = 0;
        }

      //Set Central value
      g_ratio->SetPoint(i_p, x_Type1, ratio);

      //Set the error
      Double_t error_XLow = g_Type1->GetErrorXlow(i_p);
      Double_t error_Xhigh = g_Type1->GetErrorXhigh(i_p);
      g_ratio->SetPointError(i_p, error_XLow, error_Xhigh, ratio_error, ratio_error);
    }
  return g_ratio;
}
