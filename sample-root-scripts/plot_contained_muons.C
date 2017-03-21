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


int main(){

  // open input file
  TFile *f = new TFile("select.root","READ");

  TTree  *select_tree = (TTree*)f->Get("select_tree");
  Float_t start_x, start_y, start_z;
  Float_t stop_x, stop_y, stop_z;
  Float_t p_x, p_y, p_z;
  Float_t costheta, mom;
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
  select_tree->SetBranchAddress("selected",&selected);

  double energy;
  double mom_min = select_tree->GetMinimum("mom");
  double mom_max = select_tree->GetMaximum("mom");


  double muMeV = 105.66;

  double Emin = sqrt(pow(mom_min,2) + pow(muMeV,2));
  double Emax = sqrt(pow(mom_max,2) + pow(muMeV,2));

  std::clog << " mom_min " << mom_min << " Emin " << Emin << " mom_max " << mom_max << " Emax " << Emax << std::endl;

  TFile *of = new TFile("contained_muons.root", "RECREATE");
  TH2F * h_contained_muons_all = new TH2F("h_contained_muons_all","h_contained_muons_all; E [MeV]; cos(#theta)", 100, Emin, Emax, 100, -1., 1.);
  TH2F * h_contained_muons_selected = new TH2F("h_contained_muons_selected","h_contained_muons_selected; E [MeV]; cos(#theta)", 100, Emin, Emax, 100, -1., 1.);


  for(size_t ientry=0; ientry<select_tree->GetEntries(); ientry++){
    select_tree->GetEntry(ientry);
    energy = sqrt(pow(mom,2) + pow(muMeV,2));

    h_contained_muons_all->Fill(energy, costheta);
    if( selected )
      h_contained_muons_selected->Fill(energy, costheta);

  }

  TH2F * h_contained_muons = (TH2F*)h_contained_muons_selected->Clone("h_contained_muons");
  h_contained_muons->Divide(h_contained_muons_all);

  of->Write();

  return 1;

}


