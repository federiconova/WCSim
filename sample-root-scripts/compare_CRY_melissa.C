#include <algorithm>
#include <iostream>
#include <string>

int compare_CRY_melissa(){

  TFile input_cry("cosmic_muons_e61_detector/overburden_4m/muon_hit_times.root","READ");
  TFile input_melissa("cosmic_muons_e61_detector/overburden_0m_melissa/muon_hit_times.root","READ");

  TH1F * h_n_hits_cry = (TH1F*)input_cry.Get("h_n_hits");
  h_n_hits_cry->SetLineColor(kBlack);
  h_n_hits_cry->SetLineWidth(2);
  h_n_hits_cry->SetTitle("CRY");
  TH1F * h_n_hits_melissa = (TH1F*)input_melissa.Get("h_n_hits");
  h_n_hits_melissa->SetLineColor(kRed);
  h_n_hits_melissa->SetLineWidth(2);
  h_n_hits_melissa->SetTitle("CORSIKA");

  double n_hits_cry = h_n_hits_cry->GetEntries();
  double n_hits_melissa = h_n_hits_melissa->GetEntries();

  double n_hits_min = std::min(n_hits_cry, n_hits_melissa);

  h_n_hits_cry->Scale(n_hits_min/n_hits_cry);
  h_n_hits_melissa->Scale(n_hits_min/n_hits_melissa);


  TH1F * h_times_cry = (TH1F*)input_cry.Get("h_times");
  h_times_cry->SetLineColor(kBlack);
  h_times_cry->SetLineWidth(2);
  h_times_cry->SetTitle("CRY");
  TH1F * h_times_melissa = (TH1F*)input_melissa.Get("h_times");
  h_times_melissa->SetLineColor(kRed);
  h_times_melissa->SetLineWidth(2);
  h_times_melissa->SetTitle("CORSIKA");


  double mean_time_cry = h_times_cry->GetMean();
  double mean_time_melissa = h_times_melissa->GetMean();

  int nbins = 140;
  double min = -50.;
  double max = 100.;

  TH1F * h_times_cry_translated = new TH1F("h_times_cry_translated",h_times_cry->GetTitle(),nbins,min,max);
  h_times_cry_translated->SetLineColor(h_times_cry->GetLineColor());
  h_times_cry_translated->SetLineWidth(h_times_cry->GetLineWidth());
  for(int i=0; i<=h_times_cry->GetXaxis()->GetNbins(); i++){
    h_times_cry_translated->Fill(h_times_cry->GetXaxis()->GetBinCenter(i) - mean_time_cry,h_times_cry->GetBinContent(i));
  }
  TH1F * h_times_melissa_translated = new TH1F("h_times_melissa_translated",h_times_melissa->GetTitle(),nbins,min,max);
  h_times_melissa_translated->SetLineColor(h_times_melissa->GetLineColor());
  h_times_melissa_translated->SetLineWidth(h_times_melissa->GetLineWidth());
  for(int i=0; i<=h_times_melissa->GetXaxis()->GetNbins(); i++){
    h_times_melissa_translated->Fill(h_times_melissa->GetXaxis()->GetBinCenter(i) - mean_time_melissa,h_times_melissa->GetBinContent(i));
    //    h_times_melissa_translated->SetBinContent(i,h_times_melissa->GetBinContent(i));
  }

  h_times_cry_translated->Rebin(2);
  h_times_melissa_translated->Rebin(2);

  double n_times_cry = h_times_cry->GetEntries();
  double n_times_melissa = h_times_melissa->GetEntries();

  double n_times_min = std::min(n_times_cry, n_times_melissa);

  h_times_cry_translated->Scale(n_times_min/n_times_cry);
  h_times_melissa_translated->Scale(n_times_min/n_times_melissa);


  TH1F * h_tube_id_cry = (TH1F*)input_cry.Get("h_tube_id");
  h_tube_id_cry->SetLineColor(kBlack);
  h_tube_id_cry->SetLineWidth(2);
  h_tube_id_cry->SetTitle("CRY");
  TH1F * h_tube_id_melissa = (TH1F*)input_melissa.Get("h_tube_id");
  h_tube_id_melissa->SetLineColor(kRed);
  h_tube_id_melissa->SetLineWidth(2);
  h_tube_id_melissa->SetTitle("CORSIKA");

  double n_tube_id_cry = h_tube_id_cry->GetEntries();
  double n_tube_id_melissa = h_tube_id_melissa->GetEntries();

  double n_tube_id_min = std::min(n_tube_id_cry, n_tube_id_melissa);

  h_tube_id_cry->Scale(n_tube_id_min/n_tube_id_cry);
  h_tube_id_melissa->Scale(n_tube_id_min/n_tube_id_melissa);




  TH1F * h_charge_max_hit_cry = (TH1F*)input_cry.Get("h_charge_max_hit");
  h_charge_max_hit_cry->SetLineColor(kBlack);
  h_charge_max_hit_cry->SetLineWidth(2);
  h_charge_max_hit_cry->SetTitle("CRY");
  TH1F * h_charge_max_hit_melissa = (TH1F*)input_melissa.Get("h_charge_max_hit");
  h_charge_max_hit_melissa->SetLineColor(kRed);
  h_charge_max_hit_melissa->SetLineWidth(2);
  h_charge_max_hit_melissa->SetTitle("CORSIKA");

  double n_charge_max_hit_cry = h_charge_max_hit_cry->GetEntries();
  double n_charge_max_hit_melissa = h_charge_max_hit_melissa->GetEntries();

  double n_charge_max_hit_min = std::min(n_charge_max_hit_cry, n_charge_max_hit_melissa);

  h_charge_max_hit_cry->Scale(n_charge_max_hit_min/n_charge_max_hit_cry);
  h_charge_max_hit_melissa->Scale(n_charge_max_hit_min/n_charge_max_hit_melissa);


  

  TFile output("cosmic_rays_compare_CRY_melissa.root","RECREATE");

  h_n_hits_cry->Write();
  h_n_hits_melissa->Write();

  h_times_cry->Write();
  h_times_melissa->Write();

  h_times_cry_translated->Write();
  h_times_melissa_translated->Write();

  h_tube_id_cry->Write();
  h_tube_id_melissa->Write();

  h_charge_max_hit_cry->Write();
  h_charge_max_hit_melissa->Write();


  exit(0);

  return 1;
}
