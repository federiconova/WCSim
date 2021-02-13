#include <iostream>
#include <TH1F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TPaletteAxis.h>
#include <cfloat>
#include <limits>
#include <algorithm>
#include <math.h>       /* log10 */
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TInterpreter.h>

std::string number_with_digits(int i, int n);

int n_digits = 5;
double time_bin_size = 10.; // ns
double _time_min=150.;
double _time_max=750.;
bool use_time_range = true;

int main(){

#if 1
    gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<vector<float> > >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
    gInterpreter->GenerateDictionary("vector<vector<vector<int> > >", "vector");
#endif

  int ievent = 0; 
  int itrigger = 0;

  // open input and output files
  TFile *f = new TFile("output.root","READ");
  TFile *of = new TFile("plot.root","RECREATE");

  // geometry tree
  TTree * geom_tree = (TTree*)f->Get("geom_tree");
  Float_t detector_length, detector_radius, pmt_radius;
  Int_t number_of_pmts, number_of_pmts_OD;
  geom_tree->SetBranchAddress("detector_length",&detector_length);
  geom_tree->SetBranchAddress("detector_radius",&detector_radius);
  geom_tree->SetBranchAddress("pmt_radius",&pmt_radius);
  geom_tree->SetBranchAddress("number_of_pmts",&number_of_pmts);
  geom_tree->SetBranchAddress("number_of_pmts_OD",&number_of_pmts_OD);
  geom_tree->GetEntry(0);
  detector_radius = 3800.;
  std::clog << " detector_length " << detector_length << " detector_radius " << detector_radius << " pmt_radius " << pmt_radius << " number_of_pmts " << number_of_pmts << " number_of_pmts_OD " << number_of_pmts_OD << std::endl;

  double r_limit=detector_radius;
  double z_limit=detector_length/2.;
  double z_limit_with_source = 4400.;
  double r_limit_with_source = 5200.;

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


  // all pmts tree OD
  TTree * all_pmts_tree_OD = (TTree*)f->Get("all_pmts_OD_tree");
  Int_t pmt_number_OD, pmt_location_OD;
  Float_t pmt_ux_OD, pmt_uy_OD, pmt_uz_OD, pmt_x_OD, pmt_y_OD, pmt_z_OD;
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_number",&pmt_number_OD);
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_location",&pmt_location_OD);
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_ux",&pmt_ux_OD);
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_uy",&pmt_uy_OD);
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_uz",&pmt_uz_OD);
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_x",&pmt_x_OD);
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_y",&pmt_y_OD);
  all_pmts_tree_OD->SetBranchAddress("pmt_OD_z",&pmt_z_OD);

  // primary events tree
  TTree * primary_events_tree = (TTree*)f->Get("primary_events_tree");

  std::vector<Int_t> *trigger_number = 0; std::vector<Int_t> * trigger_date = 0; std::vector<Int_t> * trigger_mode = 0; std::vector<Int_t> * trigger_vec_rec_number = 0; std::vector<Int_t> * trigger_jmu = 0; std::vector<Int_t> * trigger_jp = 0; std::vector<Int_t> * trigger_npar = 0; std::vector<Int_t> * trigger_ntrack = 0; std::vector<Int_t> * trigger_number_raw_hits = 0; std::vector<Int_t> * trigger_number_raw_hits_OD = 0; std::vector<Int_t> * trigger_number_digitized_hits = 0; std::vector<Int_t> * trigger_number_digitized_hits_OD = 0; std::vector<Int_t> * trigger_number_times = 0; std::vector<Int_t> * trigger_nvertex = 0;
  std::vector<std::vector<Int_t> > * trigger_vtxvol = 0; 

  std::vector<Float_t> * trigger_sum_q = 0;
  std::vector<std::vector<Float_t> > * trigger_vtx_x = 0; std::vector<std::vector<Float_t> > * trigger_vtx_y = 0; std::vector<std::vector<Float_t> > * trigger_vtx_z = 0;

  std::vector<std::vector<Int_t> > * track_ipnu = 0; std::vector<std::vector<Int_t> > *   track_parent_type = 0; std::vector<std::vector<Int_t> > *   track_flag = 0; std::vector<std::vector<Int_t> > *   track_start_volume = 0; std::vector<std::vector<Int_t> > *   track_stop_volume = 0; std::vector<std::vector<Int_t> > *   track_id = 0;

  std::vector<std::vector<Float_t> > * track_ux = 0; std::vector<std::vector<Float_t> > *   track_uy = 0; std::vector<std::vector<Float_t> > *   track_uz = 0; std::vector<std::vector<Float_t> > *   track_M = 0; std::vector<std::vector<Float_t> > *   track_P = 0; std::vector<std::vector<Float_t> > *   track_E = 0; std::vector<std::vector<Float_t> > *   track_px = 0; std::vector<std::vector<Float_t> > *   track_py = 0; std::vector<std::vector<Float_t> > *   track_pz = 0; std::vector<std::vector<Float_t> > *   track_stop_x = 0; std::vector<std::vector<Float_t> > *   track_stop_y = 0; std::vector<std::vector<Float_t> > *   track_stop_z = 0; std::vector<std::vector<Float_t> > *   track_start_x = 0; std::vector<std::vector<Float_t> > *   track_start_y = 0; std::vector<std::vector<Float_t> > *   track_start_z = 0; std::vector<std::vector<Float_t> > *   track_time = 0;

  std::vector<std::vector<Int_t> > * raw_hit_tube_id = 0; std::vector<std::vector<Int_t> > * raw_hit_tube_times_indexes = 0; std::vector<std::vector<Int_t> > * raw_hit_tube_pe = 0;
  std::vector<std::vector<std::vector<Float_t> > > * raw_hit_times = 0;
  std::vector<std::vector<std::vector<Int_t> > > * raw_hit_parent_ids = 0;

  std::vector<std::vector<Int_t> > * raw_hit_OD_tube_id = 0; std::vector<std::vector<Int_t> > * raw_hit_OD_tube_times_indexes = 0; std::vector<std::vector<Int_t> > * raw_hit_OD_tube_pe = 0;
  std::vector<std::vector<std::vector<Float_t> > > * raw_hit_OD_times = 0;
  std::vector<std::vector<std::vector<Int_t> > > * raw_hit_OD_parent_ids = 0;

  std::vector<std::vector<Int_t> > * digitized_hit_tube_id = 0;
  std::vector<std::vector<Float_t> > * digitized_hit_Q = 0; std::vector<std::vector<Float_t> > * digitized_hit_time = 0;
  std::vector<std::vector<std::vector<Int_t> > > * digitized_hit_photon_ids = 0;

  std::vector<std::vector<Int_t> > * digitized_hit_OD_tube_id = 0;
  std::vector<std::vector<Float_t> > * digitized_hit_OD_Q = 0; std::vector<std::vector<Float_t> > * digitized_hit_OD_time = 0;
  std::vector<std::vector<std::vector<Int_t> > > * digitized_hit_OD_photon_ids = 0;


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
  primary_events_tree->SetBranchAddress("trigger_number_raw_hits_OD",&trigger_number_raw_hits_OD);
  primary_events_tree->SetBranchAddress("trigger_number_digitized_hits_OD",&trigger_number_digitized_hits_OD);
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

  primary_events_tree->SetBranchAddress("raw_hit_OD_tube_id",&raw_hit_OD_tube_id);
  primary_events_tree->SetBranchAddress("raw_hit_OD_tube_times_indexes",&raw_hit_OD_tube_times_indexes);
  primary_events_tree->SetBranchAddress("raw_hit_OD_tube_pe",&raw_hit_OD_tube_pe);
  primary_events_tree->SetBranchAddress("raw_hit_OD_times",&raw_hit_OD_times);
  primary_events_tree->SetBranchAddress("raw_hit_OD_parent_ids",&raw_hit_OD_parent_ids);

  primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  primary_events_tree->SetBranchAddress("digitized_hit_Q",&digitized_hit_Q);
  primary_events_tree->SetBranchAddress("digitized_hit_time",&digitized_hit_time);
  primary_events_tree->SetBranchAddress("digitized_hit_photon_ids",&digitized_hit_photon_ids);

  primary_events_tree->SetBranchAddress("digitized_hit_OD_tube_id",&digitized_hit_OD_tube_id);
  primary_events_tree->SetBranchAddress("digitized_hit_OD_Q",&digitized_hit_OD_Q);
  primary_events_tree->SetBranchAddress("digitized_hit_OD_time",&digitized_hit_OD_time);
  primary_events_tree->SetBranchAddress("digitized_hit_OD_photon_ids",&digitized_hit_OD_photon_ids);

  // get primary event
  primary_events_tree->GetEvent(ievent); 

  // define space binning so that different pmts are sure to end up in different bins
  int nbins_x = (int)(3.*detector_radius/pmt_radius);
  int nbins_y = nbins_x;
  int nbins_z = (int)(3.*detector_length/(2.*pmt_radius));
  int nbins_phi  = (int)(3.*acos(-1.)*detector_radius/pmt_radius);


  // get number of hits for trigger in the event
  int number_of_raw_cherenkov_hits = trigger_number_raw_hits->at(itrigger);
  int number_of_digitized_cherenkov_hits = (digitized_hit_tube_id->at(itrigger)).size();
  int number_of_raw_cherenkov_hits_OD = trigger_number_raw_hits_OD->at(itrigger);
  int number_of_digitized_cherenkov_hits_OD = (digitized_hit_OD_tube_id->at(itrigger)).size();


  TH3F PMT_x_y_z("PMT_x_y_z","PMT_x_y_z; x [cm]; y [cm]; z [cm]",
		 nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
		 nbins_y,-(r_limit*1.1),(r_limit*1.1),
		 nbins_z,-(z_limit*1.1),z_limit);
  PMT_x_y_z.SetFillColor(kBlack);
  PMT_x_y_z.SetMarkerColor(kBlack);
  PMT_x_y_z.SetMarkerStyle(1);
  for(int ipmt = 0; ipmt < all_pmts_tree->GetEntries(); ipmt++){
    all_pmts_tree->GetEntry(ipmt);
    PMT_x_y_z.Fill(pmt_x, pmt_y, pmt_z);
  }


  TH3F PMT_OD_x_y_z("PMT_OD_x_y_z","OD pmt; x [cm]; y [cm]; z [cm]",
		 nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
		 nbins_y,-(r_limit*1.1),(r_limit*1.1),
		 nbins_z,-(z_limit*1.1),z_limit);
  PMT_OD_x_y_z.SetFillColor(kBlack);
  PMT_OD_x_y_z.SetMarkerColor(kBlack);
  PMT_OD_x_y_z.SetMarkerStyle(1);
  TH2F PMT_OD_x_y("PMT_OD_x_y","OD pmt; x [cm]; y [cm]",
		 nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
		 nbins_y,-(r_limit*1.1),(r_limit*1.1));
  PMT_OD_x_y.SetFillColor(kBlack);
  PMT_OD_x_y.SetMarkerColor(kBlack);
  PMT_OD_x_y.SetMarkerStyle(1);

  TH1F PMT_OD_z("PMT_OD_z","OD pmt; z [cm]",
				       nbins_z,-(z_limit*1.1),(z_limit*1.2));
  PMT_OD_z.SetLineColor(kBlack);
  PMT_OD_z.SetLineWidth(2);

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    PMT_OD_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
    PMT_OD_x_y.Fill(pmt_x_OD, pmt_y_OD);
    PMT_OD_z.Fill(pmt_z_OD);
  }

  TH3F *h_trigger_vtx_x_y_z = new TH3F("h_trigger_vtx_x_y_z","muon starting point; x [cm]; y [cm]; z [cm]",
				       		 nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
				       nbins_y,-(r_limit*1.1),(r_limit*1.1),
				       nbins_z,-(z_limit*1.1),z_limit);
  h_trigger_vtx_x_y_z->SetFillColor(kBlack);
  h_trigger_vtx_x_y_z->SetMarkerColor(kBlue);
  h_trigger_vtx_x_y_z->SetMarkerStyle(2);

  TH2F *h_trigger_vtx_x_y = new TH2F("h_trigger_vtx_x_y","muon starting point; x [cm]; y [cm]",
				     nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
				       nbins_y,-(r_limit*1.1),(r_limit*1.1));
  h_trigger_vtx_x_y->SetFillColor(kBlack);
  h_trigger_vtx_x_y->SetMarkerColor(kBlue);
  h_trigger_vtx_x_y->SetMarkerStyle(1);

  TH1F *h_trigger_vtx_z = new TH1F("h_trigger_vtx_z","muon starting point; z [cm]",
				       nbins_z,-(z_limit*1.1),(z_limit*1.2));
  h_trigger_vtx_z->SetLineColor(kBlue);
  h_trigger_vtx_z->SetLineWidth(2);

  TH3F *h_hit_x_y_z = new TH3F("h_hit_x_y_z","h_hit_x_y_z",
			       nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
			       nbins_y,-(r_limit*1.1),(r_limit*1.1),
			       nbins_z,-(z_limit*1.1),z_limit);
  h_hit_x_y_z->SetFillColor(kBlack);
  h_hit_x_y_z->SetMarkerColor(kRed);
  h_hit_x_y_z->SetMarkerStyle(2);

  TH2F *h_hit_x_y = new TH2F("h_hit_x_y","h_hit_x_y",
			       nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
			       nbins_y,-(r_limit*1.1),(r_limit*1.1));
  h_hit_x_y->SetFillColor(kBlack);
  h_hit_x_y->SetMarkerColor(kRed);
  h_hit_x_y->SetMarkerStyle(2);

  TH1F h_digitized_hit_OD_Q("h_digitized_hit_OD_Q","h_digitized_hit_OD_Q; hit charge [p.e.]; n of events",100,-1.,300.);

  TH1F h_digitized_hit_OD_Q_top_tubes("h_digitized_hit_OD_Q_top_tubes","top; hit charge [p.e.]; n of events",100,-1.,300.);
  h_digitized_hit_OD_Q_top_tubes.SetLineWidth(2);
  h_digitized_hit_OD_Q_top_tubes.SetLineColor(kRed);
  TH1F h_digitized_hit_OD_Q_side_tubes("h_digitized_hit_OD_Q_side_tubes","side; hit charge [p.e.]; n of events",100,-1.,300.);
  h_digitized_hit_OD_Q_side_tubes.SetLineWidth(2);
  h_digitized_hit_OD_Q_side_tubes.SetLineColor(kBlack);
  TH1F h_digitized_hit_OD_Q_bottom_tubes("h_digitized_hit_OD_Q_bottom_tubes","bottom; hit charge [p.e.]; n of events",100,-1.,300.);
  h_digitized_hit_OD_Q_bottom_tubes.SetLineWidth(2);
  h_digitized_hit_OD_Q_bottom_tubes.SetLineColor(kBlue);

  TH1F h_charge_vs_z("h_charge_vs_z","h_charge_vs_z; z [cm]; max charge [p.e.]",
		     				       nbins_z,-(z_limit*1.1),(z_limit*1.1));
  h_charge_vs_z.SetLineWidth(2);
  h_charge_vs_z.SetLineColor(kBlack);

  TH1F h_charge("h_charge","h_charge; charge [p.e.]",
		800,-0.5,799.5);
  h_charge.SetLineWidth(2);
  h_charge.SetLineColor(kBlack);

  TH1F h_hit_time("h_hit_time","h_hit_time; hit time [ns]",
		  2000,-0.5,1999.5);
  h_hit_time.SetLineWidth(2);
  h_hit_time.SetLineColor(kBlack);

  TH1F h_total_charge("h_total_charge","h_total_charge; charge [p.e.]",
		      800,1,-1);
  h_total_charge.SetLineWidth(2);
  h_total_charge.SetLineColor(kBlack);

  TH2F h_n_hits_vs_radius("h_n_hits_vs_radius","h_n_hits_vs_radius; r [cm]; n hits",
			  nbins_z,0,1.2*detector_radius,
			  100,0,2000);
  h_n_hits_vs_radius.SetLineWidth(2);
  h_n_hits_vs_radius.SetLineColor(kBlack);

  double max_nhits = 1000.; // dynamic range
  max_nhits = 99.5; // horizontal muon
  double min_nhits = -0.5;
  min_nhits=1;
  max_nhits=-1;

  TH1F h_digitized_nhits_top_tubes("h_digitized_nhits_top_tubes","top; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_top_tubes.SetLineWidth(2);
  h_digitized_nhits_top_tubes.SetLineColor(kRed);
  TH1F h_digitized_nhits_side_tubes("h_digitized_nhits_side_tubes","side; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_side_tubes.SetLineWidth(2);
  h_digitized_nhits_side_tubes.SetLineColor(kBlack);
  TH1F h_digitized_nhits_bottom_tubes("h_digitized_nhits_bottom_tubes","bottom; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_bottom_tubes.SetLineWidth(2);
  h_digitized_nhits_bottom_tubes.SetLineColor(kBlue);
  TH1F h_digitized_nhits_all_tubes("h_digitized_nhits_all_tubes","all; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_all_tubes.SetLineWidth(2);
  h_digitized_nhits_all_tubes.SetLineColor(kBlue);
  TH1F h_digitized_nhits_nonzero_tubes("h_digitized_nhits_nonzero_tubes","nonzero; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_nonzero_tubes.SetLineWidth(2);
  h_digitized_nhits_nonzero_tubes.SetLineColor(kBlue);

  TH2F h_npes_vs_digitized_nhits("h_npes_vs_digitized_nhits","npes vs n hits; n of hits; n of pes",100,min_nhits,max_nhits,800,1,-1);
  TH2F h_nhits_vs_muon_energy("h_nhits_vs_muon_energy","n hits vs muon energy; muon energy [MeV]; n of hits",100,1,-1,100,min_nhits,max_nhits);

  TH2F h_hit_position_top_tubes("h_hit_position_top_tubes","top; x; y", 20,-(r_limit*1.1),(r_limit*1.1), 20,-(r_limit*1.1),(r_limit*1.1));
  TH2F h_hit_position_bottom_tubes("h_hit_position_bottom_tubes","bottom; x; y", 20,-(r_limit*1.1),(r_limit*1.1), 20,-(r_limit*1.1),(r_limit*1.1));
  TH2F h_hit_position_side_tubes("h_hit_position_side_tubes","side; phi; z", 300,-180.,180., 40,-(z_limit*1.1),(z_limit*1.1));
  
  TH1F h_distance_pmt_vertex("h_distance_pmt_vertex","distance PMT vertex; distance PMT vertex [cm]",100,1,-1);
  h_distance_pmt_vertex.SetLineWidth(2);
  h_distance_pmt_vertex.SetLineColor(kBlack);

  TH1F h_nhits_OD_cluster("h_nhits_OD_cluster","nhits OD cluster; nhits OD cluster",100,1,-1);
  h_nhits_OD_cluster.SetLineWidth(2);
  h_nhits_OD_cluster.SetLineColor(kBlack);

  TH1F h_npes_OD_cluster("h_npes_OD_cluster","npes OD cluster; npes OD cluster",100,1,-1);
  h_npes_OD_cluster.SetLineWidth(2);
  h_npes_OD_cluster.SetLineColor(kBlack);

  TH3F *h_muon_start_x_y_z = new TH3F("h_muon_start_x_y_z","h_muon_start_x_y_z",
				       		 nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
				       nbins_y,-(r_limit*1.1),(r_limit*1.1),
				       nbins_z,-(z_limit*1.1),z_limit*1.1);
  h_muon_start_x_y_z->SetFillColor(kBlack);
  h_muon_start_x_y_z->SetMarkerColor(kBlue);
  h_muon_start_x_y_z->SetMarkerStyle(2);

  TH2F *h_muon_start_x_y = new TH2F("h_muon_start_x_y","h_muon_start_x_y",
				     nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
				       nbins_y,-(r_limit*1.1),(r_limit*1.1));
  h_muon_start_x_y->SetFillColor(kBlack);
  h_muon_start_x_y->SetMarkerColor(kBlue);
  h_muon_start_x_y->SetMarkerStyle(1);

  TH1F *h_muon_start_z = new TH1F("h_muon_start_z","h_muon_start_z",
				       nbins_z,-(z_limit*1.1),(z_limit*1.2));
  h_muon_start_z->SetLineColor(kBlack);
  h_muon_start_z->SetLineWidth(2);

  TH3F *h_muon_stop_x_y_z = new TH3F("h_muon_stop_x_y_z","h_muon_stop_x_y_z",
				       		 nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
				       nbins_y,-(r_limit*1.1),(r_limit*1.1),
				       nbins_z,-(z_limit*2),z_limit*1.1);
  h_muon_stop_x_y_z->SetFillColor(kBlack);
  h_muon_stop_x_y_z->SetMarkerColor(kBlue);
  h_muon_stop_x_y_z->SetMarkerStyle(2);

  TH2F *h_muon_stop_x_y = new TH2F("h_muon_stop_x_y","h_muon_stop_x_y",
				     nbins_x,-(r_limit_with_source*1.1),(r_limit*1.1),
				       nbins_y,-(r_limit*1.1),(r_limit*1.1));
  h_muon_stop_x_y->SetFillColor(kBlack);
  h_muon_stop_x_y->SetMarkerColor(kBlue);
  h_muon_stop_x_y->SetMarkerStyle(1);

  TH1F *h_muon_stop_z = new TH1F("h_muon_stop_z","h_muon_stop_z",
				       nbins_z,-(z_limit*2),(z_limit*1.1));
  h_muon_stop_z->SetLineColor(kBlack);
  h_muon_stop_z->SetLineWidth(2);

  TH3F *h_muon_direction_x_y_z = new TH3F("h_muon_direction_x_y_z","h_muon_direction_x_y_z",
					  nbins_x,-1,1,
					  nbins_y,-1,1,
					  nbins_z,-1,1);
  h_muon_direction_x_y_z->SetFillColor(kBlack);
  h_muon_direction_x_y_z->SetMarkerColor(kBlue);
  h_muon_direction_x_y_z->SetMarkerStyle(2);

  TH2F *h_muon_direction_x_y = new TH2F("h_muon_direction_x_y","h_muon_direction_x_y",
					nbins_x,-1,1,
					nbins_y,-1,1);
  h_muon_direction_x_y->SetFillColor(kBlack);
  h_muon_direction_x_y->SetMarkerColor(kBlue);
  h_muon_direction_x_y->SetMarkerStyle(1);

  TH1F *h_muon_direction_z = new TH1F("h_muon_direction_z","h_muon_direction_z; muon cos#theta",
				      nbins_z,-1,1);
  h_muon_direction_z->SetLineColor(kBlack);
  h_muon_direction_z->SetLineWidth(2);

  TH3F *h_muon_momentum_x_y_z = new TH3F("h_muon_momentum_x_y_z","h_muon_momentum_x_y_z",
					  nbins_x,1,-1,
					  nbins_y,1,-1,
					  nbins_z,1,-1);
  h_muon_momentum_x_y_z->SetFillColor(kBlack);
  h_muon_momentum_x_y_z->SetMarkerColor(kBlue);
  h_muon_momentum_x_y_z->SetMarkerStyle(2);

  TH2F *h_muon_momentum_x_y = new TH2F("h_muon_momentum_x_y","h_muon_momentum_x_y",
					nbins_x,1,-1,
					nbins_y,1,-1);
  h_muon_momentum_x_y->SetFillColor(kBlack);
  h_muon_momentum_x_y->SetMarkerColor(kBlue);
  h_muon_momentum_x_y->SetMarkerStyle(1);

  TH1F *h_muon_momentum_z = new TH1F("h_muon_momentum_z","h_muon_momentum_z",
				      nbins_z,1,-1);
  h_muon_momentum_z->SetLineColor(kBlack);
  h_muon_momentum_z->SetLineWidth(2);

  TH1F h_muon_M("h_muon_M","h_muon_M;muon M",100,1,-1);
  h_muon_M.SetLineWidth(2);
  h_muon_M.SetLineColor(kBlack);

  TH1F h_muon_E("h_muon_E","h_muon_E;muon E [MeV]",100,1,-1);
  h_muon_E.SetLineWidth(2);
  h_muon_E.SetLineColor(kBlack);

  TH1F h_muon_time("h_muon_time","h_muon_time;muon time",100,1,-1);
  h_muon_time.SetLineWidth(2);
  h_muon_time.SetLineColor(kBlack);

  TH1F h_muon_phi("h_muon_phi","h_muon_phi;muon #phi",100,1,-1);
  h_muon_phi.SetLineWidth(2);
  h_muon_phi.SetLineColor(kBlack);

  int tube_id;
  double charge;
  int ibin;
  double maxcharge;
  double vtx_radius, vtx_x, vtx_y, vtx_z;
  double hit_time;
  int modulo = 100;
  double cluster_radius=800.;  // cm
  double distance_pmt_vertex;
  double nhits_OD_cluster;
  double npes_OD_cluster;
  double muon_energy;

  for(int ievent=0; ievent<primary_events_tree->GetEntries(); ievent++){

    if( ievent%modulo == 0 )
      std::cout << " event " << ievent << " of " << primary_events_tree->GetEntries() << std::endl;
    // loop on primary events
    primary_events_tree->GetEvent(ievent); 

    for(size_t itrigger=0; itrigger<trigger_ntrack->size(); itrigger++){
      // loop on triggers in the event

      for(int itrack=0; itrack<trigger_ntrack->at(itrigger); itrack++){
	// loop on tracks in the event
	if( track_ipnu->at(itrigger).at(itrack) == 13 && track_M->at(itrigger).at(itrack) > 0 ){ // muon
	  h_muon_start_x_y_z->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	  h_muon_start_x_y->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack));
	  h_muon_start_z->Fill(track_start_z->at(itrigger).at(itrack));

	  h_muon_stop_x_y_z->Fill(track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack));
	  h_muon_stop_x_y->Fill(track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack));
	  h_muon_stop_z->Fill(track_stop_z->at(itrigger).at(itrack));

	  h_muon_direction_x_y_z->Fill(track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack));
	  h_muon_direction_x_y->Fill(track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack));
	  h_muon_direction_z->Fill(track_uz->at(itrigger).at(itrack));

	  h_muon_momentum_x_y_z->Fill(track_px->at(itrigger).at(itrack),track_py->at(itrigger).at(itrack),track_pz->at(itrigger).at(itrack));
	  h_muon_momentum_x_y->Fill(track_px->at(itrigger).at(itrack),track_py->at(itrigger).at(itrack));
	  h_muon_momentum_z->Fill(track_pz->at(itrigger).at(itrack));

	  h_muon_M.Fill(track_M->at(itrigger).at(itrack));
	  muon_energy = track_E->at(itrigger).at(itrack);
	  h_muon_E.Fill(muon_energy);
	  h_muon_time.Fill(track_time->at(itrigger).at(itrack));
	  h_muon_phi.Fill(atan2(track_uy->at(itrigger).at(itrack),track_ux->at(itrigger).at(itrack))*180./acos(-1.));
	}
      }

      for(int ivertex=0; ivertex<trigger_nvertex->at(itrigger); ivertex++){
	// loop on vertices in the event

	vtx_x = trigger_vtx_x->at(itrigger).at(ivertex);
	vtx_y = trigger_vtx_y->at(itrigger).at(ivertex);
	vtx_z = trigger_vtx_z->at(itrigger).at(ivertex);

	h_trigger_vtx_x_y_z->Fill(trigger_vtx_x->at(itrigger).at(ivertex), trigger_vtx_y->at(itrigger).at(ivertex), trigger_vtx_z->at(itrigger).at(ivertex));

	h_trigger_vtx_x_y->Fill(trigger_vtx_x->at(itrigger).at(ivertex), trigger_vtx_y->at(itrigger).at(ivertex));
	h_trigger_vtx_z->Fill(trigger_vtx_z->at(itrigger).at(ivertex));

	vtx_radius = pow(trigger_vtx_x->at(itrigger).at(ivertex),2) + pow(trigger_vtx_y->at(itrigger).at(ivertex),2);
	vtx_radius = sqrt(vtx_radius);
	h_n_hits_vs_radius.Fill(vtx_radius, trigger_number_digitized_hits_OD->at(itrigger));

      }

      int nhits_top = 0;
      int nhits_bottom = 0;
      int nhits_side = 0;
      int nhits_all = 0;
      double total_charge = 0.;
      nhits_OD_cluster = 0;
      npes_OD_cluster = 0;


      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_OD_tube_id->at(itrigger)).size(); idigitizedhit++){
	// loop on digitized hits in the trigger
	
	tube_id = (digitized_hit_OD_tube_id->at(itrigger)).at(idigitizedhit);
	all_pmts_tree_OD->GetEntry(tube_id - 1);

	charge = (digitized_hit_OD_Q->at(itrigger)).at(idigitizedhit);
	total_charge += charge;
	h_charge.Fill(charge);

	hit_time = (digitized_hit_OD_time->at(itrigger)).at(idigitizedhit);
	h_hit_time.Fill(hit_time);

	h_digitized_hit_OD_Q.Fill(charge);
	if( pmt_location_OD == 5 ){
	  h_digitized_hit_OD_Q_top_tubes.Fill(charge);
	  h_hit_position_top_tubes.Fill(pmt_x_OD, pmt_y_OD);
	}
	else if( pmt_location_OD == 4 ){
	  h_digitized_hit_OD_Q_side_tubes.Fill(charge);
	  h_hit_position_side_tubes.Fill(atan2(pmt_x_OD,pmt_y_OD)*180/3.14, pmt_z_OD);
	}
	else if( pmt_location_OD == 3 ){
	  h_digitized_hit_OD_Q_bottom_tubes.Fill(charge);
	  h_hit_position_bottom_tubes.Fill(pmt_x_OD, pmt_y_OD);
	}

	h_hit_x_y_z->Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	h_hit_x_y->Fill(pmt_x_OD, pmt_y_OD);

	if( pmt_location_OD == 5 )
	  nhits_top ++;
	else if( pmt_location_OD == 4 )
	  nhits_side ++;
	else if( pmt_location_OD == 3 )
	  nhits_bottom ++;

	ibin = h_charge_vs_z.FindBin(pmt_z_OD);
	maxcharge = h_charge_vs_z.GetBinContent(ibin);
	if( charge > maxcharge )
	  h_charge_vs_z.SetBinContent(ibin,charge);

	distance_pmt_vertex = sqrt(pow(pmt_y_OD - vtx_y,2) + pow(pmt_z_OD - vtx_z,2));
	h_distance_pmt_vertex.Fill(distance_pmt_vertex);

	if( distance_pmt_vertex <= cluster_radius ){
	  nhits_OD_cluster ++;
	  npes_OD_cluster += charge;
	}

      }

      h_total_charge.Fill(total_charge);
      h_nhits_OD_cluster.Fill(nhits_OD_cluster);
      h_npes_OD_cluster.Fill(npes_OD_cluster);


      h_digitized_nhits_top_tubes.Fill(nhits_top);
      h_digitized_nhits_side_tubes.Fill(nhits_side);
      h_digitized_nhits_bottom_tubes.Fill(nhits_bottom);

      nhits_all = nhits_top + nhits_side + nhits_bottom;

      h_digitized_nhits_all_tubes.Fill(nhits_all);
      h_npes_vs_digitized_nhits.Fill(nhits_all, total_charge);
      h_nhits_vs_muon_energy.Fill(muon_energy,nhits_all);
      if( nhits_all > 0 )
	h_digitized_nhits_nonzero_tubes.Fill(nhits_all);

    }
  }


  f->Close();


  of->cd();
  PMT_x_y_z.Write();
  PMT_OD_x_y_z.Write();
  PMT_OD_x_y.Write();
  PMT_OD_z.Write();
  h_digitized_hit_OD_Q.Write();
  h_digitized_hit_OD_Q_top_tubes.Write();
  h_digitized_hit_OD_Q_side_tubes.Write();
  h_digitized_hit_OD_Q_bottom_tubes.Write();
  h_trigger_vtx_x_y_z->Write();
  h_trigger_vtx_x_y->Write();
  h_trigger_vtx_z->Write();
  h_charge_vs_z.Write();
  h_charge.Write();
  h_hit_time.Write();
  h_total_charge.Write();
  h_n_hits_vs_radius.Write();
  h_distance_pmt_vertex.Write();
  h_nhits_OD_cluster.Write();
  h_npes_OD_cluster.Write();
  h_digitized_nhits_top_tubes.Write();
  h_digitized_nhits_side_tubes.Write();
  h_digitized_nhits_bottom_tubes.Write();
  h_digitized_nhits_all_tubes.Write();
  h_digitized_nhits_nonzero_tubes.Write();
  h_npes_vs_digitized_nhits.Write();
  h_nhits_vs_muon_energy.Write();
  h_hit_position_top_tubes.Write();
  h_hit_position_bottom_tubes.Write();
  h_hit_position_side_tubes.Write();
  h_hit_x_y_z->Write();
  h_hit_x_y->Write();
  h_muon_M.Write();
  h_muon_E.Write();
  h_muon_start_x_y_z->Write();
  h_muon_start_x_y->Write();
  h_muon_start_z->Write();
  h_muon_stop_x_y_z->Write();
  h_muon_stop_x_y->Write();
  h_muon_stop_z->Write();
  h_muon_direction_x_y_z->Write();
  h_muon_direction_x_y->Write();
  h_muon_direction_z->Write();
  h_muon_momentum_x_y_z->Write();
  h_muon_momentum_x_y->Write();
  h_muon_momentum_z->Write();
  h_muon_time.Write();
  h_muon_phi.Write();

  delete trigger_number ; delete  trigger_date ; delete  trigger_mode ; delete  trigger_vtxvol ; delete  trigger_vec_rec_number ; delete  trigger_jmu ; delete  trigger_jp ; delete  trigger_npar ; delete  trigger_ntrack ; delete  trigger_number_raw_hits ; delete  trigger_number_digitized_hits ; delete  trigger_number_raw_hits_OD ; delete  trigger_number_digitized_hits_OD ; delete  trigger_number_times ;

  delete  trigger_vtx_x ; delete  trigger_vtx_y ; delete  trigger_vtx_z ; delete  trigger_sum_q ;

  delete  track_ipnu ; delete    track_parent_type ; delete    track_flag ; delete    track_start_volume ; delete    track_stop_volume ; delete    track_id ;

  delete  track_ux ; delete    track_uy ; delete    track_uz ; delete    track_M ; delete    track_P ; delete    track_E ; delete    track_px ; delete    track_py ; delete    track_pz ; delete    track_stop_x ; delete    track_stop_y ; delete    track_stop_z ; delete    track_start_x ; delete    track_start_y ; delete    track_start_z ; delete    track_time ;

  delete raw_hit_tube_id; delete raw_hit_tube_times_indexes; delete raw_hit_tube_pe; delete raw_hit_times; delete raw_hit_parent_ids;

  delete digitized_hit_tube_id;  delete digitized_hit_Q;  delete digitized_hit_time;  

  delete raw_hit_OD_tube_id; delete raw_hit_OD_tube_times_indexes; delete raw_hit_OD_tube_pe; delete raw_hit_OD_times; delete raw_hit_OD_parent_ids;

  delete digitized_hit_OD_tube_id;  delete digitized_hit_OD_Q;  delete digitized_hit_OD_time;  

  return 1;

}

std::string number_with_digits(int i, int n){

  std::string output="";

  int number_of_non_zero_digits = (int)(log10(i));

  for(int j=0; j<(n - number_of_non_zero_digits - 1); j++)
    output+="0";

  char c[20];
  sprintf(c,"%d",i);
  std::string sval(c);
  output+=sval;


  return output;

}


