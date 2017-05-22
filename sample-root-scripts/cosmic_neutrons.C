#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

int main(){

  // open input file
  TFile *f = new TFile("output_cosmic_muons.root","READ");
  
  // create output file
  TFile *of = new TFile("cosmic_neutrons.root", "RECREATE");


  // geometry tree
  TTree * geom_tree = (TTree*)f->Get("geom_tree");
  Float_t detector_length, detector_radius, pmt_radius;
  Int_t number_of_pmts;
  geom_tree->SetBranchAddress("detector_length",&detector_length);
  geom_tree->SetBranchAddress("detector_radius",&detector_radius);
  geom_tree->SetBranchAddress("pmt_radius",&pmt_radius);
  geom_tree->SetBranchAddress("number_of_pmts",&number_of_pmts);
  geom_tree->GetEntry(0);
  std::clog << " detector_length " << detector_length << " detector_radius " << detector_radius << " pmt_radius " << pmt_radius << " number_of_pmts " << number_of_pmts << std::endl;

  int nbins = 100;
  double xmin = -10000.;
  double xmax = 10000.;
  double zmin = -12000.;
  double zmax = 12000.;

  // all pmts tree
  TTree * all_pmts_tree = (TTree*)f->Get("all_pmts_tree");
  Int_t pmt_number, pmt_location;
  Float_t pmt_ux, pmt_uy, pmt_uz, pmt_x, pmt_y, pmt_z;
  all_pmts_tree->SetBranchAddress("pmt_number",&pmt_number);
  all_pmts_tree->SetBranchAddress("pmt_location",&pmt_location);
  all_pmts_tree->SetBranchAddress("pmt_ux",&pmt_ux);
  all_pmts_tree->SetBranchAddress("pmt_uy",&pmt_uy);
  all_pmts_tree->SetBranchAddress("pmt_uz",&pmt_uz);
  all_pmts_tree->SetBranchAddress("pmt_x",&pmt_x);
  all_pmts_tree->SetBranchAddress("pmt_y",&pmt_y);
  all_pmts_tree->SetBranchAddress("pmt_z",&pmt_z);
  for(int ipmt = 0; ipmt < all_pmts_tree->GetEntries(); ipmt++){
    all_pmts_tree->GetEntry(ipmt);
    // std::clog << " pmt_number " << pmt_number << " pmt_location " << pmt_location << " pmt dir: (" << pmt_ux << ", " << pmt_uy << ", " << pmt_uz << "), pos: (" << pmt_x << ", " << pmt_y << ", " << pmt_z << ")" << std::endl;
  }


  // primary events tree
  TTree * primary_events_tree = (TTree*)f->Get("primary_events_tree");

  std::vector<Int_t> *trigger_number = 0; std::vector<Int_t> * trigger_date = 0; std::vector<Int_t> * trigger_mode = 0; std::vector<Int_t> * trigger_vec_rec_number = 0; std::vector<Int_t> * trigger_jmu = 0; std::vector<Int_t> * trigger_jp = 0; std::vector<Int_t> * trigger_npar = 0; std::vector<Int_t> * trigger_ntrack = 0; std::vector<Int_t> * trigger_number_raw_hits = 0; std::vector<Int_t> * trigger_number_digitized_hits = 0; std::vector<Int_t> * trigger_number_times = 0; std::vector<Int_t> * trigger_nvertex = 0;
  std::vector<std::vector<Int_t> > * trigger_vtxvol = 0; 

  std::vector<Float_t> * trigger_sum_q = 0;
  std::vector<std::vector<Float_t> > * trigger_vtx_x = 0; std::vector<std::vector<Float_t> > * trigger_vtx_y = 0; std::vector<std::vector<Float_t> > * trigger_vtx_z = 0;

  std::vector<std::vector<Int_t> > * track_ipnu = 0; std::vector<std::vector<Int_t> > *   track_parent_type = 0; std::vector<std::vector<Int_t> > *   track_flag = 0; std::vector<std::vector<Int_t> > *   track_start_volume = 0; std::vector<std::vector<Int_t> > *   track_stop_volume = 0; std::vector<std::vector<Int_t> > *   track_id = 0;

  std::vector<std::vector<Float_t> > * track_ux = 0; std::vector<std::vector<Float_t> > *   track_uy = 0; std::vector<std::vector<Float_t> > *   track_uz = 0; std::vector<std::vector<Float_t> > *   track_M = 0; std::vector<std::vector<Float_t> > *   track_P = 0; std::vector<std::vector<Float_t> > *   track_E = 0; std::vector<std::vector<Float_t> > *   track_px = 0; std::vector<std::vector<Float_t> > *   track_py = 0; std::vector<std::vector<Float_t> > *   track_pz = 0; std::vector<std::vector<Float_t> > *   track_stop_x = 0; std::vector<std::vector<Float_t> > *   track_stop_y = 0; std::vector<std::vector<Float_t> > *   track_stop_z = 0; std::vector<std::vector<Float_t> > *   track_start_x = 0; std::vector<std::vector<Float_t> > *   track_start_y = 0; std::vector<std::vector<Float_t> > *   track_start_z = 0; std::vector<std::vector<Float_t> > *   track_time = 0;

  std::vector<std::vector<Int_t> > * raw_hit_tube_id = 0; std::vector<std::vector<Int_t> > * raw_hit_tube_times_indexes = 0; std::vector<std::vector<Int_t> > * raw_hit_tube_pe = 0;
  std::vector<std::vector<std::vector<Float_t> > > * raw_hit_times = 0;
  std::vector<std::vector<std::vector<Int_t> > > * raw_hit_parent_ids = 0;

  std::vector<std::vector<Int_t> > * digitized_hit_tube_id = 0;
  std::vector<std::vector<Float_t> > * digitized_hit_Q = 0; std::vector<std::vector<Float_t> > * digitized_hit_time = 0;
  std::vector<std::vector<std::vector<Int_t> > > * digitized_hit_photon_ids = 0;

  primary_events_tree->SetBranchAddress("trigger_number",&trigger_number);
  primary_events_tree->SetBranchAddress("trigger_date",&trigger_date);
  primary_events_tree->SetBranchAddress("trigger_mode",&trigger_mode);
  primary_events_tree->SetBranchAddress("trigger_vtxvol",&trigger_vtxvol);
  primary_events_tree->SetBranchAddress("trigger_vtx_x",&trigger_vtx_x); 
  primary_events_tree->SetBranchAddress("trigger_vtx_y",&trigger_vtx_y);
  primary_events_tree->SetBranchAddress("trigger_vtx_z",&trigger_vtx_z);
  primary_events_tree->SetBranchAddress("trigger_vec_rec_number",&trigger_vec_rec_number);
  primary_events_tree->SetBranchAddress("trigger_jmu",&trigger_jmu);
  primary_events_tree->SetBranchAddress("trigger_jp",&trigger_jp); 
  primary_events_tree->SetBranchAddress("trigger_npar",&trigger_npar);
  primary_events_tree->SetBranchAddress("trigger_ntrack",&trigger_ntrack);
  primary_events_tree->SetBranchAddress("trigger_nvertex",&trigger_nvertex);
  primary_events_tree->SetBranchAddress("trigger_number_raw_hits",&trigger_number_raw_hits);
  primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits);
  primary_events_tree->SetBranchAddress("trigger_sum_q",&trigger_sum_q);
  primary_events_tree->SetBranchAddress("trigger_number_times",&trigger_number_times);

  primary_events_tree->SetBranchAddress("track_ipnu",&track_ipnu); 
  primary_events_tree->SetBranchAddress("track_parent_type",&track_parent_type); 
  primary_events_tree->SetBranchAddress("track_ux",&track_ux); 
  primary_events_tree->SetBranchAddress("track_uy",&track_uy);
  primary_events_tree->SetBranchAddress("track_uz",&track_uz);
  primary_events_tree->SetBranchAddress("track_px",&track_px); 
  primary_events_tree->SetBranchAddress("track_py",&track_py);
  primary_events_tree->SetBranchAddress("track_pz",&track_pz);
  primary_events_tree->SetBranchAddress("track_flag",&track_flag); 
  primary_events_tree->SetBranchAddress("track_M",&track_M); 
  primary_events_tree->SetBranchAddress("track_P",&track_P); 
  primary_events_tree->SetBranchAddress("track_E",&track_E); 
  primary_events_tree->SetBranchAddress("track_start_volume",&track_start_volume); 
  primary_events_tree->SetBranchAddress("track_stop_volume",&track_stop_volume); 
  primary_events_tree->SetBranchAddress("track_stop_x",&track_stop_x); 
  primary_events_tree->SetBranchAddress("track_stop_y",&track_stop_y);
  primary_events_tree->SetBranchAddress("track_stop_z",&track_stop_z);
  primary_events_tree->SetBranchAddress("track_start_x",&track_start_x); 
  primary_events_tree->SetBranchAddress("track_start_y",&track_start_y);
  primary_events_tree->SetBranchAddress("track_start_z",&track_start_z);
  primary_events_tree->SetBranchAddress("track_time",&track_time); 
  primary_events_tree->SetBranchAddress("track_id",&track_id);

  primary_events_tree->SetBranchAddress("raw_hit_tube_id",&raw_hit_tube_id);
  primary_events_tree->SetBranchAddress("raw_hit_tube_times_indexes",&raw_hit_tube_times_indexes);
  primary_events_tree->SetBranchAddress("raw_hit_tube_pe",&raw_hit_tube_pe);
  primary_events_tree->SetBranchAddress("raw_hit_times",&raw_hit_times);
  primary_events_tree->SetBranchAddress("raw_hit_parent_ids",&raw_hit_parent_ids);

  primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  primary_events_tree->SetBranchAddress("digitized_hit_Q",&digitized_hit_Q);
  primary_events_tree->SetBranchAddress("digitized_hit_time",&digitized_hit_time);
  primary_events_tree->SetBranchAddress("digitized_hit_photon_ids",&digitized_hit_photon_ids);

  Int_t ncapturecount, neutroncount;
  const Int_t MAXNCAPTURES = 25;
  Double_t capt_x[MAXNCAPTURES], capt_y[MAXNCAPTURES], capt_z[MAXNCAPTURES], capt_t0[MAXNCAPTURES], capt_E[MAXNCAPTURES];
  Int_t capt_num[MAXNCAPTURES], capt_pid[MAXNCAPTURES], capt_nucleus[MAXNCAPTURES], capt_nphot[MAXNCAPTURES], capt_ngamma[MAXNCAPTURES];
  primary_events_tree->SetBranchAddress("ncapturecount",&ncapturecount);
  primary_events_tree->SetBranchAddress("neutroncount",&neutroncount);
  primary_events_tree->SetBranchAddress("capt_x",capt_x);
  primary_events_tree->SetBranchAddress("capt_y",capt_y);
  primary_events_tree->SetBranchAddress("capt_z",capt_z);
  primary_events_tree->SetBranchAddress("capt_t0",capt_t0);
  primary_events_tree->SetBranchAddress("capt_E",capt_E);
  primary_events_tree->SetBranchAddress("capt_num",capt_num);
  primary_events_tree->SetBranchAddress("capt_pid",capt_pid);
  primary_events_tree->SetBranchAddress("capt_nucleus",capt_nucleus);
  primary_events_tree->SetBranchAddress("capt_nphot",capt_nphot);
  primary_events_tree->SetBranchAddress("capt_ngamma",capt_ngamma);
  Int_t nneutrons;
  primary_events_tree->SetBranchAddress("nneutrons",&nneutrons);

  double max_time = 2.e3;


  int tube_id;
  int parent_id;
  int ipnu;
  double n_of_muons_that_make_hit = 0.;
  double n_of_captured_neutrons = 0.;
  double n_of_captured_neutrons_without_muons_that_make_hit = 0.;
  bool muon_makes_hit;


  for(int ievent=0; ievent<primary_events_tree->GetEntries(); ievent++){
    // loop on primary events
    primary_events_tree->GetEvent(ievent); 
    muon_makes_hit = false;

    for(size_t itrigger=0; itrigger<trigger_ntrack->size(); itrigger++){
      // loop on triggers in the event

      for(int itrack=0; itrack<trigger_ntrack->at(itrigger); itrack++){
	ipnu = (track_ipnu->at(itrigger)).at(itrack);
	if( abs(ipnu) == 13 ){
	  muon_makes_hit = true;
	  break;
	}
      }

      if( muon_makes_hit )
	n_of_muons_that_make_hit ++;
      
      if( ncapturecount ){
	n_of_captured_neutrons ++;
	if( !muon_makes_hit ){
	  n_of_captured_neutrons_without_muons_that_make_hit ++;
	}
      }

    }
  }

  std::clog << " n_of_muons_that_make_hit " << n_of_muons_that_make_hit << std::endl;
  std::clog << " n_of_captured_neutrons " << n_of_captured_neutrons << std::endl;
  std::clog << " n_of_captured_neutrons_without_muons_that_make_hit " << n_of_captured_neutrons_without_muons_that_make_hit << " ( = " << n_of_captured_neutrons_without_muons_that_make_hit/n_of_captured_neutrons << " )" << std::endl;

  of->cd();


  of->Write();

  delete trigger_number ; delete  trigger_date ; delete  trigger_mode ; delete  trigger_vtxvol ; delete  trigger_vec_rec_number ; delete  trigger_jmu ; delete  trigger_jp ; delete  trigger_npar ; delete  trigger_ntrack ; delete  trigger_number_raw_hits ; delete  trigger_number_digitized_hits ; delete  trigger_number_times ;

  delete  trigger_vtx_x ; delete  trigger_vtx_y ; delete  trigger_vtx_z ; delete  trigger_sum_q ;

  delete  track_ipnu ; delete    track_parent_type ; delete    track_flag ; delete    track_start_volume ; delete    track_stop_volume ; delete    track_id ;

  delete  track_ux ; delete    track_uy ; delete    track_uz ; delete    track_M ; delete    track_P ; delete    track_E ; delete    track_px ; delete    track_py ; delete    track_pz ; delete    track_stop_x ; delete    track_stop_y ; delete    track_stop_z ; delete    track_start_x ; delete    track_start_y ; delete    track_start_z ; delete    track_time ;

  delete raw_hit_tube_id; delete raw_hit_tube_times_indexes; delete raw_hit_tube_pe; delete raw_hit_times; delete raw_hit_parent_ids;

  delete digitized_hit_tube_id;  delete digitized_hit_Q;  delete digitized_hit_time;   delete digitized_hit_photon_ids;


  return 1;

}


