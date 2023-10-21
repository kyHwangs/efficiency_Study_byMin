void canvas_margin(TCanvas *c1, TPad *c1_up, TPad *c1_down);
void hist_axis(TH1D *hist, TH1D *hist_compare);

void drawTrackQuals(TString filePath = "./Outputs_vRun3_02/hist-vRun3_02-trackQual.root"){

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
    "lostPixelHits",
    "lostPixelBarrelHits",
    "lostPixelEndcapHits",
    "lostStripHits",
    "lostStripTIBHits",
    "lostStripTIDHits",
    "lostStripTOBHits",
    "lostStripTECHits",
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

  vector<TString> L3types = {
    "iterL3OI",
    "iterL3IOFromL2",
    "iterL3FromL2",
    "iterL3IOFromL1", //Track-only

    "muon",
    "iterL3MuonNoID",
    "iterL3Muon",
  };

  TFile *hists = new TFile(filePath);

  for(unsigned i=0; i<L3types.size(); ++i){
    TString Dir = "./plots_trackQual/"+L3types.at(i)+"/";
    if (gSystem->mkdir(Dir,kTRUE) != -1) gSystem->mkdir(Dir,kTRUE);

    vector<TString> vars;
    if (i==3) vars = trkonly_vars;
    else if(i<3) vars = trk_vars;
    else vars = muon_vars;

    for (unsigned j=0; j<vars.size(); ++j) {
      TString name = L3types.at(i)+"/h_"+L3types.at(i)+"_"+vars.at(j);

      TH1D *out_before = (TH1D*)hists->Get(name+"_out_before");
      TH1D *out_after = (TH1D*)hists->Get(name+"_out_after");
      TH1D *in_before = (TH1D*)hists->Get(name+"_in_before");
      TH1D *in_after = (TH1D*)hists->Get(name+"_in_after");

      out_before->Scale(1./out_before->Integral(1, out_before->GetXaxis()->GetNbins()));
      out_after->Scale(1./out_after->Integral(1, out_after->GetXaxis()->GetNbins()));
      in_before->Scale(1./in_before->Integral(1, in_before->GetXaxis()->GetNbins()));
      in_after->Scale(1./in_after->Integral(1, in_after->GetXaxis()->GetNbins()));

      TCanvas *c_var = new TCanvas("c_var", "", 800, 800);
      c_var->Draw();
      TPad *c1_up = new TPad("c1_up", "", 0, 0.25, 1, 1);
      TPad *c1_down = new TPad("c1_down", "", 0, 0, 1, 0.25);
      canvas_margin(c_var,c1_up,c1_down);
      c1_down->SetGridx();
      c1_down->SetGridy();
      c1_up->Draw();
      c1_down->Draw();

      c1_up->cd();
      out_after->SetLineColor(kBlack);
      out_after->SetMarkerColor(kBlack);
      out_after->SetMarkerStyle(25);
      out_after->SetMarkerSize(1);
      out_before->SetLineColor(kBlue);
      out_before->SetMarkerColor(kBlue);
      out_before->SetMarkerStyle(25);
      out_before->SetMarkerSize(1);

      in_after->SetLineColor(kGreen+2);
      in_after->SetMarkerColor(kGreen+2);
      in_after->SetMarkerStyle(20);
      in_after->SetMarkerSize(1);
      in_before->SetLineColor(kRed);
      in_before->SetMarkerColor(kRed);
      in_before->SetMarkerStyle(20);
      in_before->SetMarkerSize(1);

      TH1D *ratio_out = (TH1D*)out_after->Clone(); ratio_out->Divide(out_before);
      TH1D *ratio_in = (TH1D*)in_after->Clone(); ratio_in->Divide(in_before);

      out_after->Draw("p");
      out_before->Draw("psame");
      in_after->Draw("psame");
      in_before->Draw("psame");

      TLegend *lg = new TLegend(0.6, 0.75, 0.95, 0.95);
      lg->SetFillStyle(0);
      lg->SetBorderSize(0);
      lg->AddEntry(out_after,  "Outside the dip, After Issue", "lp");
      lg->AddEntry(out_before, "Outside the dip, Before Issue", "lp");
      lg->AddEntry(in_after,   "Inside the dip, After Issue", "lp");
      lg->AddEntry(in_before,  "Inside the dip, Before Issue", "lp");
      lg->Draw();

      out_after->SetTitle("");
      out_after->GetYaxis()->SetTitle("Normalized");
      //out_after->GetYaxis()->SetRangeUser(0.6, 1.1);
      c1_down->cd();

      ratio_out->Draw("p");
      ratio_in->Draw("psame");
      ratio_out->GetXaxis()->SetTitle(vars.at(j));
      ratio_out->GetYaxis()->SetTitle("After/Before");
      ratio_out->GetYaxis()->SetRangeUser(0, 3);

      hist_axis(out_after, ratio_out);

      c_var->cd();

      TLatex latex_CMSPriliminary, latex_selection;
      latex_CMSPriliminary.SetNDC();
      latex_selection.SetNDC();

      latex_CMSPriliminary.SetTextSize(0.045);
      latex_CMSPriliminary.DrawLatex(0.115, 0.965, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
      latex_selection.SetTextSize(0.033); //original size 0.03 and position (0.125,0.92)
      latex_selection.DrawLatex(0.75, 0.965, L3types.at(i));

      c_var->SaveAs(Dir+vars.at(j)+".png");
      delete c_var;

      // Draw 2D lostHits vs eta, phi
      if(vars.at(j).Contains("lost")){
        TH2D *before = (TH2D*)hists->Get(name+"_2D_before");
        before->SetTitle(L3types.at(i)+"_"+vars.at(j)+"_Before");
        TH2D *after = (TH2D*)hists->Get(name+"_2D_after");
        after->SetTitle(L3types.at(i)+"_"+vars.at(j)+"_After");

        TCanvas *c_before = new TCanvas("c_before", "", 1200, 800);
        TCanvas *c_after = new TCanvas("c_after", "", 1200, 800);

        c_before->cd();
        c_before->Draw();
        before->Draw("colz");
        c_before->SaveAs(Dir+vars.at(j)+"_2D_before.png");
        delete c_before;

        c_after->cd();
        c_after->Draw();
        after->Draw("colz");
        c_after->SaveAs(Dir+vars.at(j)+"_2D_after.png");
        delete c_after;
      }
    }
  }
}

void canvas_margin(TCanvas *c1, TPad *c1_up, TPad *c1_down){
  c1_up->SetTopMargin( 0.06 );     //0.07 (JaeSung's original number set)
  c1_up->SetBottomMargin( 0.02 );  //0.02
  c1_up->SetLeftMargin( 0.11 );    //0.15
  c1_up->SetRightMargin( 0.02 );   //0.03

  c1_down->SetTopMargin( 0.03 );   //0.03
  c1_down->SetBottomMargin( 0.26 );//0.4
  c1_down->SetLeftMargin( 0.11 );  //0.15
  c1_down->SetRightMargin( 0.02 ); //0.03
  //c1_down->SetGridx();
  //c1_down->SetGridy();
  
  c1->SetTopMargin( 0.05 );        //0.05
  c1->SetBottomMargin( 0.13 );     //0.13
  c1->SetRightMargin( 0.05 );      //0.05
  c1->SetLeftMargin( 0.16 );       //0.16

  gStyle->SetOptStat(0);
}
void hist_axis(TH1D *hist, TH1D *hist_compare){

  hist->SetTitle("");

  //==== top plot
  hist->GetYaxis()->SetLabelSize(0.05);           //0.05
  hist->GetYaxis()->SetTitleSize(0.07);           //0.07
  hist->GetYaxis()->SetTitleOffset(0.75);         //1.02
  //==== hide x-axis for top plot
  hist->GetXaxis()->SetLabelSize(0);

  //==== bottom plot
  hist_compare->SetTitle("");
  hist_compare->GetXaxis()->SetLabelSize(0.10);   //0.10
  hist_compare->GetXaxis()->SetTitleSize(0.13);   //0.15
  hist_compare->GetXaxis()->SetTitleOffset(0.85);
  hist_compare->GetYaxis()->SetLabelSize(0.08);   //0.08
  hist_compare->GetYaxis()->SetTitleSize(0.13);   //0.12
  hist_compare->GetYaxis()->SetTitleOffset(0.3);  //0.5
  hist_compare->SetFillColorAlpha(45,0.35);

}
