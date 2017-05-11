#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <math.h>       /* log10 */
#include <TF1.h>
#include <TRandom.h>
#include <TProfile.h>


int main(){

  // open input file
  TFile *f = new TFile("select_muons.root","READ");

  TTree  *select_tree = (TTree*)f->Get("select_tree");
  Float_t start_x, start_y, start_z;
  Float_t stop_x, stop_y, stop_z;
  Float_t p_x, p_y, p_z;
  Float_t costheta, mom;
  Int_t n_digitized_hits;
  Int_t selected;
  select_tree->SetBranchAddress("start_x",&start_x);
  select_tree->SetBranchAddress("start_y",&start_y);
  select_tree->SetBranchAddress("start_z",&start_z);
  select_tree->SetBranchAddress("stop_x",&stop_x);
  select_tree->SetBranchAddress("stop_y",&stop_y);
  select_tree->SetBranchAddress("stop_z",&stop_z);
  select_tree->SetBranchAddress("costheta",&costheta);
  select_tree->SetBranchAddress("p_x",&p_x);
  select_tree->SetBranchAddress("p_y",&p_y);
  select_tree->SetBranchAddress("p_z",&p_z);
  select_tree->SetBranchAddress("mom",&mom);
  select_tree->SetBranchAddress("n_digitized_hits",&n_digitized_hits);
  select_tree->SetBranchAddress("selected",&selected);

  double energy;
  double mom_min = select_tree->GetMinimum("mom");
  double mom_max = select_tree->GetMaximum("mom");

  Int_t n_min = select_tree->GetMinimum("n_digitized_hits");
  Int_t n_max = select_tree->GetMaximum("n_digitized_hits");


  double muMeV = 105.66;

  double Emin = sqrt(pow(mom_min,2) + pow(muMeV,2)) - muMeV;
  double Emax = sqrt(pow(mom_max,2) + pow(muMeV,2)) - muMeV;

  std::clog << " mom_min " << mom_min << " Emin " << Emin << " mom_max " << mom_max << " Emax " << Emax << std::endl;

  TFile *of = new TFile("contained_muons.root", "RECREATE");
  TH2F * h_contained_muons_all = new TH2F("h_contained_muons_all","h_contained_muons_all; K [MeV]; cos(#theta)", 50, Emin, Emax, 50, -1., 1.);
  TH2F * h_contained_muons_selected = new TH2F("h_contained_muons_selected","h_contained_muons_selected; K [MeV]; cos(#theta)", 50, Emin, Emax, 50, -1., 1.);

  TH2F * h_n_hits_vs_energy = new TH2F("h_n_hits_vs_energy","h_n_hits_vs_energy; K [MeV]; n digitized hits", 50, Emin, Emax, 50, n_min, n_max);
  TH2F * h_n_hits_vs_energy_selected = new TH2F("h_n_hits_vs_energy_selected","h_n_hits_vs_energy_selected; K [MeV]; n digitized hits", 50, Emin, Emax, 50, n_min, n_max);


  for(size_t ientry=0; ientry<select_tree->GetEntries(); ientry++){
    select_tree->GetEntry(ientry);
    energy = sqrt(pow(mom,2) + pow(muMeV,2)) - muMeV;

    h_contained_muons_all->Fill(energy, costheta);

    h_n_hits_vs_energy->Fill(energy, n_digitized_hits);

    if( selected ){
      h_contained_muons_selected->Fill(energy, costheta);
      h_n_hits_vs_energy_selected->Fill(energy, n_digitized_hits);
    }

  }

  TH2F * h_contained_muons_fraction = (TH2F*)h_contained_muons_selected->Clone("h_contained_muons_fraction");
  h_contained_muons_fraction->Divide(h_contained_muons_all);

  TProfile* p_n_hits_vs_energy = h_n_hits_vs_energy->ProfileX(); 
  p_n_hits_vs_energy->GetXaxis()->SetTitle(h_n_hits_vs_energy->GetXaxis()->GetTitle());
  p_n_hits_vs_energy->GetYaxis()->SetTitle(h_n_hits_vs_energy->GetYaxis()->GetTitle());

  TProfile* p_n_hits_vs_energy_selected = h_n_hits_vs_energy_selected->ProfileX(); 
  p_n_hits_vs_energy_selected->GetXaxis()->SetTitle(h_n_hits_vs_energy_selected->GetXaxis()->GetTitle());
  p_n_hits_vs_energy_selected->GetYaxis()->SetTitle(h_n_hits_vs_energy_selected->GetYaxis()->GetTitle());

  of->Write();

  return 1;

}


