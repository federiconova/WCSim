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
bool is_impact_horizontal_plane(double x, double y, double z,
				double ux, double uy, double uz,
				double R, double zplane, double E,
				double * impact_time, double *impact_x, double *impact_y, double *impact_z);
bool is_impact_side(double x, double y, double z,
		    double ux, double uy, double uz,
		    double R, double H, double E,
		    double * impact_time, double *impact_x, double *impact_y, double *impact_z);

bool is_crossing_horizontal_plane(double x0, double y0, double z0,
				  double x1, double y1, double z1,
				  double R, double zplane);

void is_crossing_cylindrical_shell_side(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double Rout, double ztop, double zbottom,
					int * n_intersections_out);

int point_position(double x0, double y0, double z0,
		   double Rout, double zout, double Rin, double zin);

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

  const int n_muon_topologies = 6;
  std::string topology_name[n_muon_topologies]={"OD(2) ID(2)","OD(1) ID(2)","OD(1) ID(1)","OD(2) ID(0)","OD(1) ID(0)","OD(0) ID(1)"};
  int topology_color[n_muon_topologies]={1,2,3,4,6,7};
  int muon_topology_true; 
  int muon_topology_reco; 
  // -2 -> unknown
  // -1 -> muon misses tank

  bool true_muon_enters_OD_top;

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
  // OD PMT side radius: 3300
  // so water side radius: 3400                                                                                                                                      
  // muon start max radius 3500
  // so side muons start out of water
  // OD PMT z max: 3350
  // so water z max: 3550
  // vtx z max: 3700
  // so top muons start out of water
  double ID_radius = 3250.;
  double ID_height = 3300.;
  double OD_radius = 3400.;
  double OD_height = 3550.;
  std::clog << " detector_length " << detector_length << " detector_radius " << detector_radius << " pmt_radius " << pmt_radius << " number_of_pmts " << number_of_pmts << " number_of_pmts_OD " << number_of_pmts_OD << std::endl;

  double r_limit=detector_radius;
  double z_limit=detector_length/2.;
  double z_limit_with_source = 4400.;
  double r_limit_with_source = 5200.;
  double pi = acos(-1.);

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
		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
		 nbins_y,-(r_limit*1.3),(r_limit*1.3),
		 nbins_z,-(z_limit*1.3),z_limit*1.3);
  PMT_x_y_z.SetFillColor(kBlack);
  PMT_x_y_z.SetMarkerColor(kBlack);
  PMT_x_y_z.SetMarkerStyle(1);
  for(int ipmt = 0; ipmt < all_pmts_tree->GetEntries(); ipmt++){
    all_pmts_tree->GetEntry(ipmt);
    PMT_x_y_z.Fill(pmt_x, pmt_y, pmt_z);
  }


  TH3F PMT_OD_x_y_z("PMT_OD_x_y_z","OD pmt; x [cm]; y [cm]; z [cm]",
		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
		 nbins_y,-(r_limit*1.3),(r_limit*1.3),
		 nbins_z,-(z_limit*1.3),z_limit*1.3);
  PMT_OD_x_y_z.SetFillColor(kBlack);
  PMT_OD_x_y_z.SetMarkerColor(kBlack);
  PMT_OD_x_y_z.SetMarkerStyle(1);
  TH2F PMT_OD_x_y("PMT_OD_x_y","OD pmt; x [cm]; y [cm]",
		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
		 nbins_y,-(r_limit*1.3),(r_limit*1.3));
  PMT_OD_x_y.SetFillColor(kBlack);
  PMT_OD_x_y.SetMarkerColor(kBlack);
  PMT_OD_x_y.SetMarkerStyle(1);

  TH1F PMT_OD_z("PMT_OD_z","OD pmt; z [cm]",
				       nbins_z,-(z_limit*1.3),(z_limit*1.3));
  PMT_OD_z.SetLineColor(kBlack);
  PMT_OD_z.SetLineWidth(2);

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    PMT_OD_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
    PMT_OD_x_y.Fill(pmt_x_OD, pmt_y_OD);
    PMT_OD_z.Fill(pmt_z_OD);
  }

  TH3F *h_trigger_vtx_x_y_z = new TH3F("h_trigger_vtx_x_y_z","muon starting point; x [cm]; y [cm]; z [cm]",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_trigger_vtx_x_y_z->SetFillColor(kBlack);
  h_trigger_vtx_x_y_z->SetMarkerColor(kBlue);
  h_trigger_vtx_x_y_z->SetMarkerStyle(2);

  TH2F *h_trigger_vtx_x_y = new TH2F("h_trigger_vtx_x_y","muon starting point; x [cm]; y [cm]",
				     nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3));
  h_trigger_vtx_x_y->SetFillColor(kBlack);
  h_trigger_vtx_x_y->SetMarkerColor(kBlue);
  h_trigger_vtx_x_y->SetMarkerStyle(1);

  TH1F *h_trigger_vtx_z = new TH1F("h_trigger_vtx_z","muon starting point; z [cm]",
				       nbins_z,-(z_limit*1.3),(z_limit*1.3));
  h_trigger_vtx_z->SetLineColor(kBlue);
  h_trigger_vtx_z->SetLineWidth(2);

  TH3F *h_hit_x_y_z = new TH3F("h_hit_x_y_z","h_hit_x_y_z",
			       nbins_x,-(r_limit*1.3),(r_limit*1.3),
			       nbins_y,-(r_limit*1.3),(r_limit*1.3),
			       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_hit_x_y_z->SetFillColor(kBlack);
  h_hit_x_y_z->SetMarkerColor(kRed);
  h_hit_x_y_z->SetMarkerStyle(2);

  TH2F *h_hit_x_y = new TH2F("h_hit_x_y","h_hit_x_y",
			       nbins_x,-(r_limit*1.3),(r_limit*1.3),
			       nbins_y,-(r_limit*1.3),(r_limit*1.3));
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
		     				       nbins_z,-(z_limit*1.3),(z_limit*1.3));
  h_charge_vs_z.SetLineWidth(2);
  h_charge_vs_z.SetLineColor(kBlack);

  TH1F h_charge("h_charge","h_charge; charge [p.e.]",
		800,-0.5,799.5);
  h_charge.SetLineWidth(2);
  h_charge.SetLineColor(kBlack);

  double tmax = 4000.;
  double time_offset = 950.; // ns
  TH1F h_hit_time("h_hit_time","h_hit_time; hit time [ns]",
		  2000,-0.5-time_offset,tmax-time_offset);
  h_hit_time.SetLineWidth(2);
  h_hit_time.SetLineColor(kBlack);

  double speedlight=0.3/1.33; // m/ns
  TH1F h_pathlength("h_pathlength","h_pathlength; path length [m]",
		    2000,(-0.5-time_offset)*speedlight,(tmax-time_offset)*speedlight);
  h_pathlength.SetLineWidth(2);
  h_pathlength.SetLineColor(kBlack);

  TH1F h_pathlength_top("h_pathlength_top","h_pathlength_top; path length [m]",
		    2000,(-0.5-time_offset)*speedlight,(tmax-time_offset)*speedlight);
  h_pathlength_top.SetLineWidth(2);
  h_pathlength_top.SetLineColor(kBlack);

  TH1F h_pathlength_side("h_pathlength_side","h_pathlength_side; path length [m]",
		    2000,(-0.5-time_offset)*speedlight,(tmax-time_offset)*speedlight);
  h_pathlength_side.SetLineWidth(2);
  h_pathlength_side.SetLineColor(kBlack);

  TH1F h_muon_topology_true("h_muon_topology_true","muon topology true; topology (true)",
			    n_muon_topologies+2, -2.5, n_muon_topologies-0.5);
  h_muon_topology_true.SetLineWidth(3);
  h_muon_topology_true.SetLineColor(kBlack);
  h_muon_topology_true.GetXaxis()->SetBinLabel(1,"unknown");
  h_muon_topology_true.GetXaxis()->SetBinLabel(2,"OD(0) ID(0)");
  for(int i=0; i<n_muon_topologies; i++)
    h_muon_topology_true.GetXaxis()->SetBinLabel(i+3,topology_name[i].c_str());


  TH1F h_muon_n_intersections_OD_true("h_muon_n_intersections_OD_true","muon n intersections OD true; n OD intersections",
				      10,-0.5,9.5);
  h_muon_n_intersections_OD_true.SetLineWidth(3);
  h_muon_n_intersections_OD_true.SetLineColor(kBlack);

  TH1F h_muon_n_intersections_ID_true("h_muon_n_intersections_ID_true","muon n intersections ID true; n ID intersections",
				      10,-0.5,9.5);
  h_muon_n_intersections_ID_true.SetLineWidth(3);
  h_muon_n_intersections_ID_true.SetLineColor(kBlack);

  double n_hits_limit = 300.;
  TH1F h_hit_time_few_nhits("h_hit_time_few_nhits",Form("nhits < %.0f; hit time [ns]",n_hits_limit),
		  2000,-0.5-time_offset,tmax-time_offset);
  h_hit_time_few_nhits.SetLineWidth(2);
  h_hit_time_few_nhits.SetLineColor(kBlue);

  TH1F h_hit_time_many_nhits("h_hit_time_many_nhits",Form("nhits > %.0f; hit time [ns]",n_hits_limit),
		  2000,-0.5-time_offset,tmax-time_offset);
  h_hit_time_many_nhits.SetLineWidth(2);
  h_hit_time_many_nhits.SetLineColor(kRed);

  double max_path_length = 1000;

  TH1F h_muon_path_length_impact("h_muon_path_length_impact","h_muon_path_length_impact; path length [cm]",
				    100,0,max_path_length);
  h_muon_path_length_impact.SetLineWidth(2);
  h_muon_path_length_impact.SetLineColor(kBlack);

  TH1F h_muon_path_length_impact_few_nhits("h_muon_path_length_impact_few_nhits",Form("nhits < %.0f; path length [cm]",n_hits_limit),
				    100,0,max_path_length);
  h_muon_path_length_impact_few_nhits.SetLineWidth(2);
  h_muon_path_length_impact_few_nhits.SetLineColor(kBlue);

  TH1F h_muon_path_length_impact_many_nhits("h_muon_path_length_impact_many_nhits",Form("nhits > %.0f; path length [cm]",n_hits_limit),
				    100,0,max_path_length);
  h_muon_path_length_impact_many_nhits.SetLineWidth(2);
  h_muon_path_length_impact_many_nhits.SetLineColor(kRed);

  TH1F h_total_charge("h_total_charge","h_total_charge; charge [p.e.]",
		      800,1,-1);
  h_total_charge.SetLineWidth(2);
  h_total_charge.SetLineColor(kBlack);

  TH2F h_n_hits_vs_radius("h_n_hits_vs_radius","h_n_hits_vs_radius; r [cm]; n hits",
			  nbins_z,0,1.2*detector_radius,
			  100,0,2000);
  h_n_hits_vs_radius.SetLineWidth(2);
  h_n_hits_vs_radius.SetLineColor(kBlack);

  double max_nhits = 1200.5; // dynamic range
  //  max_nhits = 99.5; // horizontal muon
  double min_nhits = -0.5;
  //  min_nhits=1;
  //  max_nhits=-1;

  TH1F h_digitized_photons_id0("h_digitized_photons_id0","top; hit id[0]; n of events",100,1,-1);
  h_digitized_photons_id0.SetLineWidth(2);
  h_digitized_photons_id0.SetLineColor(kRed);
  TH1F h_digitized_n_photons_ids("h_digitized_n_photons_ids","top; n hit ids; n of events",100,1,-1);
  h_digitized_n_photons_ids.SetLineWidth(2);
  h_digitized_n_photons_ids.SetLineColor(kRed);
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

  TH1F * h_digitized_nhits_all_tubes_by_topology[n_muon_topologies];
  TH1F * h_digitized_nhits_zoom_all_tubes_by_topology[n_muon_topologies];
  for(int i=0; i<n_muon_topologies; i++){
    h_digitized_nhits_all_tubes_by_topology[i] = new TH1F(Form("h_digitized_nhits_all_tubes_topology_%d",i),Form("%s; n of OD hits",topology_name[i].c_str()),100,min_nhits,max_nhits);
    h_digitized_nhits_all_tubes_by_topology[i]->SetLineWidth(3);
    h_digitized_nhits_all_tubes_by_topology[i]->SetLineColor(topology_color[i]);
    h_digitized_nhits_zoom_all_tubes_by_topology[i] = new TH1F(Form("h_digitized_nhits_zoom_all_tubes_topology_%d",i),Form("%s; n of OD hits",topology_name[i].c_str()),500,min_nhits,max_nhits);
    h_digitized_nhits_zoom_all_tubes_by_topology[i]->SetLineWidth(3);
    h_digitized_nhits_zoom_all_tubes_by_topology[i]->SetLineColor(topology_color[i]);
  }




  TH1F h_digitized_nhits_top_physics_tubes("h_digitized_nhits_top_physics_tubes","top; n of hits; n of events",100,1,-1);
  h_digitized_nhits_top_physics_tubes.SetLineWidth(2);
  h_digitized_nhits_top_physics_tubes.SetLineColor(kRed);
  TH1F h_digitized_nhits_side_physics_tubes("h_digitized_nhits_side_physics_tubes","side; n of hits; n of events",100,1,-1);
  h_digitized_nhits_side_physics_tubes.SetLineWidth(2);
  h_digitized_nhits_side_physics_tubes.SetLineColor(kBlack);
  TH1F h_digitized_nhits_bottom_physics_tubes("h_digitized_nhits_bottom_physics_tubes","bottom; n of hits; n of events",100,1,-1);
  h_digitized_nhits_bottom_physics_tubes.SetLineWidth(2);
  h_digitized_nhits_bottom_physics_tubes.SetLineColor(kBlue);
  TH1F h_digitized_nhits_physics_tubes("h_digitized_nhits_physics_tubes","nonzero; n of hits; n of events",100,1,-1);
  h_digitized_nhits_physics_tubes.SetLineWidth(2);
  h_digitized_nhits_physics_tubes.SetLineColor(kBlue);

  TH1F h_digitized_nhits_physics_tubes_impact_few_nhits("h_digitized_nhits_physics_tubes_impact_few_nhits",Form("nhits < %.0f; n of hits",n_hits_limit),100,0,max_nhits);
  h_digitized_nhits_physics_tubes_impact_few_nhits.SetLineWidth(2);
  h_digitized_nhits_physics_tubes_impact_few_nhits.SetLineColor(kBlue);

  TH1F h_digitized_nhits_physics_tubes_impact_many_nhits("h_digitized_nhits_physics_tubes_impact_many_nhits",Form("nhits > %.0f; n of hits",n_hits_limit),100,0,max_nhits);
  h_digitized_nhits_physics_tubes_impact_many_nhits.SetLineWidth(2);
  h_digitized_nhits_physics_tubes_impact_many_nhits.SetLineColor(kRed);

  TH2F h_npes_vs_digitized_nhits("h_npes_vs_digitized_nhits","npes vs n hits; n of hits; n of pes",100,min_nhits,max_nhits,800,1,-1);
  TH2F h_nhits_vs_muon_energy("h_nhits_vs_muon_energy","n hits vs muon energy; muon energy [MeV]; n of hits",100,1,-1,100,min_nhits,max_nhits);

  TH2F h_hit_position_top_tubes("h_hit_position_top_tubes","top; x; y", 20,-(r_limit*1.3),(r_limit*1.3), 20,-(r_limit*1.3),(r_limit*1.3));
  TH2F h_hit_position_bottom_tubes("h_hit_position_bottom_tubes","bottom; x; y", 20,-(r_limit*1.3),(r_limit*1.3), 20,-(r_limit*1.3),(r_limit*1.3));
  TH2F h_hit_position_side_tubes("h_hit_position_side_tubes","side; phi; z", 300,-180.,180., 40,-(z_limit*1.3),(z_limit*1.3));

  TH1F h_angle_muon_direction_normal_top_impact_few_nhits("h_angle_muon_direction_normal_top_impact_few_nhits",Form("nhits < %.0f; cos(#mu direction, normal to top)",n_hits_limit),100,-1,1);
  h_angle_muon_direction_normal_top_impact_few_nhits.SetLineWidth(2);
  h_angle_muon_direction_normal_top_impact_few_nhits.SetLineColor(kBlue);
  
  TH1F h_angle_muon_direction_normal_top_impact_many_nhits("h_angle_muon_direction_normal_top_impact_many_nhits",Form("nhits > %.0f; cos(#mu direction, normal to top)",n_hits_limit),100,-1,1);
  h_angle_muon_direction_normal_top_impact_many_nhits.SetLineWidth(2);
  h_angle_muon_direction_normal_top_impact_many_nhits.SetLineColor(kRed);
  
  TH1F h_angle_muon_direction_normal_side_impact_few_nhits("h_angle_muon_direction_normal_side_impact_few_nhits",Form("nhits < %.0f; cos(#mu direction, normal to side)",n_hits_limit),100,-1,1);
  h_angle_muon_direction_normal_side_impact_few_nhits.SetLineWidth(2);
  h_angle_muon_direction_normal_side_impact_few_nhits.SetLineColor(kBlue);
  
  TH1F h_angle_muon_direction_normal_side_impact_many_nhits("h_angle_muon_direction_normal_side_impact_many_nhits",Form("nhits > %.0f; cos(#mu direction, normal to side)",n_hits_limit),100,-1,1);
  h_angle_muon_direction_normal_side_impact_many_nhits.SetLineWidth(2);
  h_angle_muon_direction_normal_side_impact_many_nhits.SetLineColor(kRed);
  
  TH1F h_distance_pmt_impact("h_distance_pmt_impact","distance PMT impact point; distance PMT impact point [cm]",100,1,-1);
  h_distance_pmt_impact.SetLineWidth(2);
  h_distance_pmt_impact.SetLineColor(kBlack);
  TH1F h_distance_pmt_impact_few_nhits("h_distance_pmt_impact_few_nhits",Form("distance when nhits < %.0f; distance PMT impact point [cm]",n_hits_limit),100,1,-1);
  h_distance_pmt_impact_few_nhits.SetLineWidth(2);
  h_distance_pmt_impact_few_nhits.SetLineColor(kBlue);
  TH1F h_distance_pmt_impact_many_nhits("h_distance_pmt_impact_many_nhits",Form("distance when nhits > %.0f; distance PMT impact point [cm]",n_hits_limit),100,1,-1);
  h_distance_pmt_impact_many_nhits.SetLineWidth(2);
  h_distance_pmt_impact_many_nhits.SetLineColor(kRed);

  TH1F h_distance_pmt_center_of_mass("h_distance_pmt_center_of_mass","distance PMT - center of mass; distance PMT - center of mass [deg]",100,1,-1);
  h_distance_pmt_center_of_mass.SetLineWidth(2);
  h_distance_pmt_center_of_mass.SetLineColor(kBlack);

  TH1F * h_distance_pmt_center_of_mass_by_topology[n_muon_topologies];
  TH1D * h_muon_stop_by_tooplogy[n_muon_topologies];
  for(int i=0; i<n_muon_topologies; i++){
    h_distance_pmt_center_of_mass_by_topology[i] = new TH1F(Form("h_distance_pmt_center_of_mass_topology_%d",i),Form("%s; distance PMT - center of mass [deg]",topology_name[i].c_str()),100,1,-1);
    h_distance_pmt_center_of_mass_by_topology[i]->SetLineWidth(3);
    h_distance_pmt_center_of_mass_by_topology[i]->SetLineColor(topology_color[i]);
    h_muon_stop_by_tooplogy[i] = new TH1D(Form("h_muon_stop_by_topology_%d",i),Form("%s",topology_name[i].c_str()),3,-0.5,2.5);
    h_muon_stop_by_tooplogy[i]->SetLineWidth(3);
    h_muon_stop_by_tooplogy[i]->SetLineColor(topology_color[i]);
    h_muon_stop_by_tooplogy[i]->GetXaxis()->SetBinLabel(1,"OD");
    h_muon_stop_by_tooplogy[i]->GetXaxis()->SetBinLabel(2,"ID");
    h_muon_stop_by_tooplogy[i]->GetXaxis()->SetBinLabel(3,"out of tank");
  }



  TH3F *h_muon_start_x_y_z_impact_few_nhits = new TH3F("h_muon_start_x_y_z_impact_few_nhits",Form("when nhits < %.0f",n_hits_limit),
						       nbins_x,-(r_limit*1.3),(r_limit*1.3),
						       nbins_y,-(r_limit*1.3),(r_limit*1.3),
						       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_start_x_y_z_impact_few_nhits->SetFillColor(kBlack);
  h_muon_start_x_y_z_impact_few_nhits->SetMarkerColor(kBlue);
  h_muon_start_x_y_z_impact_few_nhits->SetMarkerStyle(2);

  TH3F *h_muon_start_x_y_z_impact_many_nhits = new TH3F("h_muon_start_x_y_z_impact_many_nhits",Form("when nhits > %.0f",n_hits_limit),
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_start_x_y_z_impact_many_nhits->SetFillColor(kBlack);
  h_muon_start_x_y_z_impact_many_nhits->SetMarkerColor(kRed);
  h_muon_start_x_y_z_impact_many_nhits->SetMarkerStyle(2);

  double cluster_radius_1=800.;  // cm
  double cluster_radius_2=1600.;
  double cluster_radius_cm=45.;  // cm
  TH1F h_nhits_OD_cluster_1("h_nhits_OD_cluster_1",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),100,min_nhits,max_nhits);
  h_nhits_OD_cluster_1.SetLineWidth(3);
  h_nhits_OD_cluster_1.SetLineColor(kBlack);
  h_nhits_OD_cluster_1.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_center_of_mass("h_nhits_OD_cluster_1_center_of_mass",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),100,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_center_of_mass.SetLineWidth(3);
  h_nhits_OD_cluster_1_center_of_mass.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_center_of_mass.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_center_of_mass_few_nhits("h_nhits_OD_cluster_1_center_of_mass_few_nhits",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),100,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_center_of_mass_few_nhits.SetLineWidth(3);
  h_nhits_OD_cluster_1_center_of_mass_few_nhits.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_center_of_mass_few_nhits.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_center_of_mass_many_nhits("h_nhits_OD_cluster_1_center_of_mass_many_nhits",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_cm),100,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_center_of_mass_many_nhits.SetLineWidth(3);
  h_nhits_OD_cluster_1_center_of_mass_many_nhits.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_center_of_mass_many_nhits.SetFillColor(kBlack);

  TH1F * h_nhits_OD_cluster_1_center_of_mass_by_topology[n_muon_topologies];
  TH1F * h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[n_muon_topologies];
  for(int i=0; i<n_muon_topologies; i++){
    h_nhits_OD_cluster_1_center_of_mass_by_topology[i] = new TH1F(Form("h_nhits_OD_cluster_1_center_of_mass_topology_%d",i),Form("cluster (< %.0f) %s; n OD hits in cluster",cluster_radius_cm, topology_name[i].c_str()),max_nhits-min_nhits+1,min_nhits,max_nhits);
    h_nhits_OD_cluster_1_center_of_mass_by_topology[i]->SetLineWidth(3);
    h_nhits_OD_cluster_1_center_of_mass_by_topology[i]->SetLineColor(topology_color[i]);
    //    h_nhits_OD_cluster_1_center_of_mass_by_topology[i]->SetFillColor(topology_color[i]);
    h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[i] = new TH1F(Form("h_nhits_OD_out_of_cluster_1_center_of_mass_topology_%d",i),Form("cluster (< %.0f) %s; n OD hits in cluster",cluster_radius_cm, topology_name[i].c_str()),100,min_nhits,max_nhits);
    h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[i]->SetLineWidth(3);
    h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[i]->SetLineColor(topology_color[i]);
  }

  double min_npes=0;
  double max_npes=4000;

  TH1F h_npes_OD_cluster_1("h_npes_OD_cluster_1",Form("cluster (< %.0f); n OD p.e.'s in cluster",cluster_radius_1),100,min_npes,max_npes);
  h_npes_OD_cluster_1.SetLineWidth(3);
  h_npes_OD_cluster_1.SetLineColor(kBlack);
  h_npes_OD_cluster_1.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_2("h_nhits_OD_cluster_2",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_2),100,min_nhits,max_nhits);
  h_nhits_OD_cluster_2.SetLineWidth(3);
  h_nhits_OD_cluster_2.SetLineColor(kBlue);
  h_nhits_OD_cluster_2.SetFillColor(kBlue);

  TH1F h_npes_OD_cluster_2("h_npes_OD_cluster_2",Form("cluster (< %.0f); n OD p.e.'s in cluster",cluster_radius_2),100,min_npes,max_npes);
  h_npes_OD_cluster_2.SetLineWidth(3);
  h_npes_OD_cluster_2.SetLineColor(kBlue);
  h_npes_OD_cluster_2.SetFillColor(kBlue);

  TH3F *h_muon_start_x_y_z = new TH3F("h_muon_start_x_y_z","h_muon_start_x_y_z",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_start_x_y_z->SetFillColor(kBlack);
  h_muon_start_x_y_z->SetMarkerColor(kBlue);
  h_muon_start_x_y_z->SetMarkerStyle(2);

  TH3F *h_muon_start_x_y_z_impact = new TH3F("h_muon_start_x_y_z_impact","h_muon_start_x_y_z_impact",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_start_x_y_z_impact->SetFillColor(kBlack);
  h_muon_start_x_y_z_impact->SetMarkerColor(kBlue);
  h_muon_start_x_y_z_impact->SetMarkerStyle(2);

  TH1F *h_muon_theta_impact = new TH1F("h_muon_theta_impact","h_muon_theta_impact",nbins_x,0.,180.);
  h_muon_theta_impact->SetLineColor(kBlack);
  h_muon_theta_impact->SetLineWidth(2);

  TH1F *h_muon_phi_impact = new TH1F("h_muon_phi_impact","h_muon_phi_impact",nbins_x,-180.,180.);
  h_muon_phi_impact->SetLineColor(kBlack);
  h_muon_phi_impact->SetLineWidth(2);

  TH1F *h_muon_theta_vertex = new TH1F("h_muon_theta_vertex","MC vertex; #theta [deg]",nbins_x,0.,180.);
  h_muon_theta_vertex->SetLineColor(kBlue);
  h_muon_theta_vertex->SetLineWidth(2);

  TH1F *h_muon_phi_vertex = new TH1F("h_muon_phi_vertex","MC vertex;#phi [deg]",nbins_x,-180.,180.);
  h_muon_phi_vertex->SetLineColor(kBlue);
  h_muon_phi_vertex->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_theta = new TH1F("h_muon_center_of_mass_theta","hits center of mass; #theta",nbins_x,0.,180.);
  h_muon_center_of_mass_theta->SetLineColor(kBlack);
  h_muon_center_of_mass_theta->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_phi = new TH1F("h_muon_center_of_mass_phi","hits center of mass; #phi",nbins_x,-180.,180.);
  h_muon_center_of_mass_phi->SetLineColor(kBlack);
  h_muon_center_of_mass_phi->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_theta_impact_residual = new TH1F("h_muon_center_of_mass_theta_impact_residual","h_muon_center_of_mass_theta_impact_residual",nbins_x,-180.,180.);
  h_muon_center_of_mass_theta_impact_residual->SetLineColor(kBlack);
  h_muon_center_of_mass_theta_impact_residual->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_phi_impact_residual = new TH1F("h_muon_center_of_mass_phi_impact_residual","h_muon_center_of_mass_phi_impact_residual",nbins_x,-180.,180.);
  h_muon_center_of_mass_phi_impact_residual->SetLineColor(kBlack);
  h_muon_center_of_mass_phi_impact_residual->SetLineWidth(2);

  TH2F *h_muon_theta_impact_scatter = new TH2F("h_muon_theta_impact_scatter","h_muon_theta_impact_scatter;MC;center of mass",nbins_x,0.,180.,nbins_x,0.,180.);
  h_muon_theta_impact_scatter->SetLineColor(kBlack);
  h_muon_theta_impact_scatter->SetLineWidth(2);

  TH2F *h_muon_phi_impact_scatter = new TH2F("h_muon_phi_impact_scatter","h_muon_phi_impact_scatter;MC;center of mass",nbins_x,-180.,180.,nbins_x,-180.,180.);
  h_muon_phi_impact_scatter->SetLineColor(kBlack);
  h_muon_phi_impact_scatter->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_vertex_theta_residual = new TH1F("h_muon_center_of_mass_vertex_theta_residual","#theta (cm - vertex)",nbins_x,-180.,180.);
  h_muon_center_of_mass_vertex_theta_residual->SetLineColor(kRed);
  h_muon_center_of_mass_vertex_theta_residual->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_vertex_theta_residual_few_nhits = new TH1F("h_muon_center_of_mass_vertex_theta_residual_few_nhits","#theta (cm - vertex)",nbins_x,-180.,180.);
  h_muon_center_of_mass_vertex_theta_residual_few_nhits->SetLineColor(kRed);
  h_muon_center_of_mass_vertex_theta_residual_few_nhits->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_vertex_theta_residual_many_nhits = new TH1F("h_muon_center_of_mass_vertex_theta_residual_many_nhits","#theta (cm - vertex)",nbins_x,-180.,180.);
  h_muon_center_of_mass_vertex_theta_residual_many_nhits->SetLineColor(kRed);
  h_muon_center_of_mass_vertex_theta_residual_many_nhits->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_vertex_phi_residual = new TH1F("h_muon_center_of_mass_vertex_phi_residual","#phi (cm - vertex)",nbins_x,-180.,180.);
  h_muon_center_of_mass_vertex_phi_residual->SetLineColor(kBlack);
  h_muon_center_of_mass_vertex_phi_residual->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_vertex_phi_residual_few_nhits = new TH1F("h_muon_center_of_mass_vertex_phi_residual_few_nhits","#phi (cm - vertex)",nbins_x,-180.,180.);
  h_muon_center_of_mass_vertex_phi_residual_few_nhits->SetLineColor(kRed);
  h_muon_center_of_mass_vertex_phi_residual_few_nhits->SetLineWidth(2);

  TH1F *h_muon_center_of_mass_vertex_phi_residual_many_nhits = new TH1F("h_muon_center_of_mass_vertex_phi_residual_many_nhits","#phi (cm - vertex)",nbins_x,-180.,180.);
  h_muon_center_of_mass_vertex_phi_residual_many_nhits->SetLineColor(kRed);
  h_muon_center_of_mass_vertex_phi_residual_many_nhits->SetLineWidth(2);

  TH2F *h_muon_center_of_mass_vertex_theta_scatter = new TH2F("h_muon_center_of_mass_vertex_theta_scatter","h_muon_center_of_mass_vertex_theta_scatter;MC;center of mass",nbins_x,0.,180.,nbins_x,0.,180.);
  h_muon_center_of_mass_vertex_theta_scatter->SetLineColor(kBlack);
  h_muon_center_of_mass_vertex_theta_scatter->SetLineWidth(2);

  TH2F *h_muon_center_of_mass_vertex_phi_scatter = new TH2F("h_muon_center_of_mass_vertex_phi_scatter","h_muon_center_of_mass_vertex_phi_scatter;MC;center of mass",nbins_x,-180.,180.,nbins_x,-180.,180.);
  h_muon_center_of_mass_vertex_phi_scatter->SetLineColor(kBlack);
  h_muon_center_of_mass_vertex_phi_scatter->SetLineWidth(2);

  TH2F *h_muon_center_of_mass_vertex_theta_scatter_few_nhits = new TH2F("h_muon_center_of_mass_vertex_theta_scatter_few_nhits","h_muon_center_of_mass_vertex_theta_scatter_few_nhits;MC;center of mass",nbins_x,0.,180.,nbins_x,0.,180.);
  h_muon_center_of_mass_vertex_theta_scatter_few_nhits->SetLineColor(kBlack);
  h_muon_center_of_mass_vertex_theta_scatter_few_nhits->SetLineWidth(2);

  TH2F *h_muon_center_of_mass_vertex_theta_scatter_many_nhits = new TH2F("h_muon_center_of_mass_vertex_theta_scatter_many_nhits","h_muon_center_of_mass_vertex_theta_scatter_many_nhits;MC;center of mass",nbins_x,0.,180.,nbins_x,0.,180.);
  h_muon_center_of_mass_vertex_theta_scatter_many_nhits->SetLineColor(kBlack);
  h_muon_center_of_mass_vertex_theta_scatter_many_nhits->SetLineWidth(2);

  TH3F *h_muon_start_x_y_z_impact_top = new TH3F("h_muon_start_x_y_z_impact_top","h_muon_start_x_y_z_impact_top",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_start_x_y_z_impact_top->SetFillColor(kBlack);
  h_muon_start_x_y_z_impact_top->SetMarkerColor(kBlue);
  h_muon_start_x_y_z_impact_top->SetMarkerStyle(2);

  TH3F *h_muon_start_x_y_z_impact_side = new TH3F("h_muon_start_x_y_z_impact_side","h_muon_start_x_y_z_impact_side",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_start_x_y_z_impact_side->SetFillColor(kBlack);
  h_muon_start_x_y_z_impact_side->SetMarkerColor(kBlue);
  h_muon_start_x_y_z_impact_side->SetMarkerStyle(2);

  TH3F *h_muon_start_x_y_z_impact_bottom = new TH3F("h_muon_start_x_y_z_impact_bottom","h_muon_start_x_y_z_impact_bottom",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_start_x_y_z_impact_bottom->SetFillColor(kBlack);
  h_muon_start_x_y_z_impact_bottom->SetMarkerColor(kBlue);
  h_muon_start_x_y_z_impact_bottom->SetMarkerStyle(2);

  TH2F *h_muon_start_x_y = new TH2F("h_muon_start_x_y","h_muon_start_x_y",
				     nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3));
  h_muon_start_x_y->SetFillColor(kBlack);
  h_muon_start_x_y->SetMarkerColor(kBlue);
  h_muon_start_x_y->SetMarkerStyle(1);

  TH2F *h_muon_start_x_y_impact_few_nhits = new TH2F("h_muon_start_x_y_impact_few_nhits",Form("when nhits < %.0f",n_hits_limit),
						     nbins_x,-(r_limit*1.3),(r_limit*1.3),
						     nbins_y,-(r_limit*1.3),(r_limit*1.3));
  h_muon_start_x_y_impact_few_nhits->SetFillColor(kBlack);
  h_muon_start_x_y_impact_few_nhits->SetMarkerColor(kBlue);
  h_muon_start_x_y_impact_few_nhits->SetMarkerStyle(1);

  TH2F *h_muon_start_x_y_impact_many_nhits = new TH2F("h_muon_start_x_y_impact_many_nhits",Form("when nhits > %.0f",n_hits_limit),
						     nbins_x,-(r_limit*1.3),(r_limit*1.3),
						     nbins_y,-(r_limit*1.3),(r_limit*1.3));
  h_muon_start_x_y_impact_many_nhits->SetFillColor(kBlack);
  h_muon_start_x_y_impact_many_nhits->SetMarkerColor(kRed);
  h_muon_start_x_y_impact_many_nhits->SetMarkerStyle(1);

  TH1F *h_muon_start_z = new TH1F("h_muon_start_z","h_muon_start_z",
				       nbins_z,-(z_limit*1.3),(z_limit*1.3));
  h_muon_start_z->SetLineColor(kBlack);
  h_muon_start_z->SetLineWidth(2);

  TH1F *h_muon_start_z_impact_few_nhits = new TH1F("h_muon_start_z_impact_few_nhits",Form("when nhits < %.0f",n_hits_limit),
						   nbins_z,-(z_limit*1.3),(z_limit*1.3));
  h_muon_start_z_impact_few_nhits->SetLineColor(kBlue);
  h_muon_start_z_impact_few_nhits->SetLineWidth(2);

  TH1F *h_muon_start_z_impact_many_nhits = new TH1F("h_muon_start_z_impact_many_nhits",Form("when nhits > %.0f",n_hits_limit),
						   nbins_z,-(z_limit*1.3),(z_limit*1.3));
  h_muon_start_z_impact_many_nhits->SetLineColor(kRed);
  h_muon_start_z_impact_many_nhits->SetLineWidth(2);

  TH3F *h_muon_stop_x_y_z = new TH3F("h_muon_stop_x_y_z","h_muon_stop_x_y_z",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*2),z_limit*1.3);
  h_muon_stop_x_y_z->SetFillColor(kBlack);
  h_muon_stop_x_y_z->SetMarkerColor(kBlue);
  h_muon_stop_x_y_z->SetMarkerStyle(2);

  TH3F * h_muon_stop_x_y_z_by_topology[n_muon_topologies];
  TH2F * h_muon_stop_x_y_by_topology[n_muon_topologies];
  TH1F * h_muon_stop_z_by_topology[n_muon_topologies];
  for(int i=0; i<n_muon_topologies; i++){
    h_muon_stop_x_y_z_by_topology[i] = new TH3F(Form("h_muon_stop_x_y_z_topology_%d",i),Form("%s",topology_name[i].c_str()),
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
						nbins_y,-(r_limit*1.3),(r_limit*1.3),
						nbins_z,-(z_limit*2),z_limit*1.3);
    h_muon_stop_x_y_z_by_topology[i]->SetMarkerStyle(2);
    h_muon_stop_x_y_z_by_topology[i]->SetMarkerColor(topology_color[i]);
    h_muon_stop_x_y_by_topology[i] = new TH2F(Form("h_muon_stop_x_y_topology_%d",i),Form("%s",topology_name[i].c_str()),
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
						nbins_y,-(r_limit*1.3),(r_limit*1.3));
    h_muon_stop_x_y_by_topology[i]->SetMarkerStyle(6);
    h_muon_stop_x_y_by_topology[i]->SetMarkerColor(topology_color[i]);
    h_muon_stop_z_by_topology[i] = new TH1F(Form("h_muon_stop_z_topology_%d",i),Form("%s",topology_name[i].c_str()), nbins_z,-(z_limit*2),(z_limit*1.3));
    h_muon_stop_z_by_topology[i]->SetLineWidth(3);
    h_muon_stop_z_by_topology[i]->SetLineColor(topology_color[i]);
  }

  TH2F *h_muon_stop_x_y = new TH2F("h_muon_stop_x_y","h_muon_stop_x_y",
				     nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3));
  h_muon_stop_x_y->SetFillColor(kBlack);
  h_muon_stop_x_y->SetMarkerColor(kBlue);
  h_muon_stop_x_y->SetMarkerStyle(1);

  TH1F *h_muon_stop_z = new TH1F("h_muon_stop_z","h_muon_stop_z",
				       nbins_z,-(z_limit*2),(z_limit*1.3));
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

  TH3F *h_muon_impact_x_y_z = new TH3F("h_muon_impact_x_y_z","h_muon_impact_x_y_z",
				       		 nbins_x,-(r_limit*1.3),(r_limit*1.3),
				       nbins_y,-(r_limit*1.3),(r_limit*1.3),
				       nbins_z,-(z_limit*1.3),z_limit*1.3);
  h_muon_impact_x_y_z->SetFillColor(kBlack);
  h_muon_impact_x_y_z->SetMarkerColor(kBlue);
  h_muon_impact_x_y_z->SetMarkerStyle(2);

  TH2F *h_muon_impact_vs_start_x = new TH2F("h_muon_impact_vs_start_x","h_muon_impact_vs_start_x",
					    nbins_x,-(r_limit*1.3),(r_limit*1.3),
					    nbins_x,-(r_limit*1.3),(r_limit*1.3));
  h_muon_impact_vs_start_x->SetFillColor(kBlack);
  h_muon_impact_vs_start_x->SetMarkerColor(kBlue);
  h_muon_impact_vs_start_x->SetMarkerStyle(1);

  TH2F *h_muon_impact_vs_start_y = new TH2F("h_muon_impact_vs_start_y","h_muon_impact_vs_start_y",
					    nbins_x,-(r_limit*1.3),(r_limit*1.3),
					    nbins_x,-(r_limit*1.3),(r_limit*1.3));
  h_muon_impact_vs_start_y->SetFillColor(kBlack);
  h_muon_impact_vs_start_y->SetMarkerColor(kBlue);
  h_muon_impact_vs_start_y->SetMarkerStyle(1);

  TH2F *h_muon_impact_vs_start_z = new TH2F("h_muon_impact_vs_start_z","h_muon_impact_vs_start_z",
					    nbins_x,-(r_limit*1.3),(r_limit*1.3),
					    nbins_x,-(r_limit*1.3),(r_limit*1.3));
  h_muon_impact_vs_start_z->SetFillColor(kBlack);
  h_muon_impact_vs_start_z->SetMarkerColor(kBlue);
  h_muon_impact_vs_start_z->SetMarkerStyle(1);

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
  double vtx_radius, vtx_x, vtx_y, vtx_z, dir_x, dir_y, dir_z;
  double hit_time;
  int modulo = 100;
  double distance_pmt_impact;
  double distance_pmt_cm;
  double nhits_OD_cluster_1;
  double npes_OD_cluster_1;
  double nhits_OD_cluster_2;
  double npes_OD_cluster_2;
  double nhits_OD_cluster_1_center_of_mass;
  double nhits_OD_cluster_1_center_of_mass_by_topology[n_muon_topologies];
  double muon_energy;
  int muon_pdg_id = 13;

  for(int ievent=0; ievent<primary_events_tree->GetEntries(); ievent++){

    if( ievent%modulo == 0 )
      std::cout << " event " << ievent << " of " << primary_events_tree->GetEntries() << std::endl;
    // loop on primary events
    primary_events_tree->GetEvent(ievent); 

    for(size_t itrigger=0; itrigger<trigger_ntrack->size(); itrigger++){
      // loop on triggers in the event

      bool impact_top = false;
      bool impact_bottom = false;
      bool impact_side = false;
      bool impact = false;
      bool OD_impact_top = false;
      bool OD_impact_bottom = false;
      bool OD_impact_side = false;
      bool OD_impact = false;
      double impact_top_time,  impact_top_x,  impact_top_y,  impact_top_z;
      double impact_bottom_time,  impact_bottom_x,  impact_bottom_y,  impact_bottom_z;
      double impact_side_time,  impact_side_x,  impact_side_y,  impact_side_z;
      double impact_time,  impact_x,  impact_y,  impact_z;
      double OD_impact_top_time,  OD_impact_top_x,  OD_impact_top_y,  OD_impact_top_z;
      double OD_impact_bottom_time,  OD_impact_bottom_x,  OD_impact_bottom_y,  OD_impact_bottom_z;
      double OD_impact_side_time,  OD_impact_side_x,  OD_impact_side_y,  OD_impact_side_z;
      double OD_impact_time,  OD_impact_x,  OD_impact_y,  OD_impact_z;

      for(int itrack=0; itrack<trigger_ntrack->at(itrigger); itrack++){
	// loop on tracks in the event
	if( track_ipnu->at(itrigger).at(itrack) == muon_pdg_id && track_M->at(itrigger).at(itrack) > 0 ){ // muon
	  h_muon_start_x_y_z->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	  h_muon_start_x_y->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack));
	  h_muon_start_z->Fill(track_start_z->at(itrigger).at(itrack));

	  h_muon_stop_x_y_z->Fill(track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack));
	  h_muon_stop_x_y->Fill(track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack));
	  h_muon_stop_z->Fill(track_stop_z->at(itrigger).at(itrack));

	  h_muon_direction_x_y_z->Fill(track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack));
	  dir_x = track_ux->at(itrigger).at(itrack);
	  dir_y = track_uy->at(itrigger).at(itrack);
	  dir_z = track_uz->at(itrigger).at(itrack);
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

	  bool is_muon_crossing_OD_top = is_crossing_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
										track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack),
										OD_radius, OD_height
										);

	  bool is_muon_crossing_OD_bottom = is_crossing_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
										track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack),
										OD_radius, -OD_height
										);

	  int n_intersections_out;
	  is_crossing_cylindrical_shell_side(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
										track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack),
					     OD_radius, OD_height, - OD_height,
					     &n_intersections_out);

	  int n_intersections_OD = 0;
	  if( is_muon_crossing_OD_top ) n_intersections_OD++;
	  if( is_muon_crossing_OD_bottom ) n_intersections_OD++;
	  n_intersections_OD += n_intersections_out;



	  bool is_muon_crossing_ID_top = is_crossing_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
										track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack),
										ID_radius, ID_height
										);

	  bool is_muon_crossing_ID_bottom = is_crossing_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
										track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack),
										ID_radius, -ID_height
										);
	  int n_intersections_in;
	  is_crossing_cylindrical_shell_side(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
										track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack),
					     ID_radius, ID_height, - ID_height,
					     &n_intersections_in);

	  int n_intersections_ID = 0;
	  if( is_muon_crossing_ID_top ) n_intersections_ID++;
	  if( is_muon_crossing_ID_bottom ) n_intersections_ID++;
	  n_intersections_ID += n_intersections_in;


	  muon_topology_true = -2;
	  if( n_intersections_OD == 0  && n_intersections_ID == 0 ) muon_topology_true = -1;
	  if( n_intersections_OD == 2  && n_intersections_ID == 2 ) muon_topology_true = 0;
	  if( n_intersections_OD == 1  && n_intersections_ID == 2 ) muon_topology_true = 1;
	  if( n_intersections_OD == 1  && n_intersections_ID == 1 ) muon_topology_true = 2;
	  if( n_intersections_OD == 2  && n_intersections_ID == 0 ) muon_topology_true = 3;
	  if( n_intersections_OD == 1  && n_intersections_ID == 0 ) muon_topology_true = 4;
	  if( n_intersections_OD == 0  && n_intersections_ID == 1 ) muon_topology_true = 5;

	  h_muon_topology_true.Fill(muon_topology_true);
	  h_muon_n_intersections_OD_true.Fill(n_intersections_OD);
	  h_muon_n_intersections_ID_true.Fill(n_intersections_ID);

	  if( muon_topology_true >= 0 ){
	    h_muon_stop_x_y_z_by_topology[muon_topology_true]->Fill(track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack));
	    h_muon_stop_x_y_by_topology[muon_topology_true]->Fill(track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack));
	    h_muon_stop_z_by_topology[muon_topology_true]->Fill(track_stop_z->at(itrigger).at(itrack));
	    h_muon_stop_by_tooplogy[muon_topology_true]->Fill(point_position(
									     track_stop_x->at(itrigger).at(itrack),track_stop_y->at(itrigger).at(itrack),track_stop_z->at(itrigger).at(itrack),
									     OD_radius, OD_height, ID_radius, ID_height
									     ));
	  }

	  impact_top = is_impact_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
						  track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
						  ID_radius, ID_height, muon_energy,
						     & impact_top_time, & impact_top_x, & impact_top_y, & impact_top_z);
	  impact_bottom = is_impact_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
						  track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
						     ID_radius, -ID_height, muon_energy,
						     & impact_bottom_time, & impact_bottom_x, & impact_bottom_y, & impact_bottom_z);
	  if( ! impact_top ) // muons always move down, so if they go through the top the impact point is there (even if later they cross the side)
	    impact_side = is_impact_side(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
					 track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
					 ID_radius, ID_height, muon_energy,
					 & impact_side_time, & impact_side_x, & impact_side_y, & impact_side_z);

	  if( impact_top && impact_bottom ) impact_bottom = false;
	  if( impact_top && impact_side ) impact_side = false;
	  if( impact_side && impact_bottom ){
	    if( impact_bottom_time < impact_side_time ) impact_side = false;
	    else impact_bottom = false;
	  }
	  if( impact_top || impact_side || impact_bottom ) impact = true;


	  OD_impact_top = is_impact_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
						  track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
						  OD_radius, OD_height, muon_energy,
						     & OD_impact_top_time, & OD_impact_top_x, & OD_impact_top_y, & OD_impact_top_z);
	  OD_impact_bottom = is_impact_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
						  track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
						     OD_radius, -OD_height, muon_energy,
						     & OD_impact_bottom_time, & OD_impact_bottom_x, & OD_impact_bottom_y, & OD_impact_bottom_z);
	  if( ! OD_impact_top ) // muons always move down, so if they go through the top the impact point is there (even if later they cross the side)
	    OD_impact_side = is_impact_side(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
					 track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
					 OD_radius, OD_height, muon_energy,
					 & OD_impact_side_time, & OD_impact_side_x, & OD_impact_side_y, & OD_impact_side_z);

	  if( OD_impact_top && OD_impact_bottom ) OD_impact_bottom = false;
	  if( OD_impact_top && OD_impact_side ) OD_impact_side = false;
	  if( OD_impact_side && OD_impact_bottom ){
	    if( OD_impact_bottom_time < OD_impact_side_time ) OD_impact_side = false;
	    else OD_impact_bottom = false;
	  }
	  if( OD_impact_top || OD_impact_side || OD_impact_bottom ) OD_impact = true;

	  if( impact ){
	    h_muon_start_x_y_z_impact->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	    if( impact_top ){
	      h_muon_start_x_y_z_impact_top->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	      impact_time = impact_top_time;
	      impact_x = impact_top_x;
	      impact_y = impact_top_y;
	      impact_z = impact_top_z;
	    }
	    if( impact_side ){
	      h_muon_start_x_y_z_impact_side->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	      impact_time = impact_side_time;
	      impact_x = impact_side_x;
	      impact_y = impact_side_y;
	      impact_z = impact_side_z;
	    }
	    if( impact_bottom ){
	      h_muon_start_x_y_z_impact_bottom->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	      impact_time = impact_bottom_time;
	      impact_x = impact_bottom_x;
	      impact_y = impact_bottom_y;
	      impact_z = impact_bottom_z;
	    }

	    h_muon_impact_x_y_z->Fill(impact_x, impact_y, impact_z);
	    h_muon_impact_vs_start_x->Fill(track_start_x->at(itrigger).at(itrack),impact_x);
	    h_muon_impact_vs_start_y->Fill(track_start_y->at(itrigger).at(itrack),impact_y);
	    h_muon_impact_vs_start_z->Fill(impact_z,impact_z);

	    h_muon_phi_impact->Fill(180./pi*atan2(impact_y,impact_x));
	    h_muon_theta_impact->Fill(180./pi*atan2(sqrt(pow(impact_x,2) + pow(impact_y,2)),impact_z));

	  }
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
      int nhits_top_physics = 0;
      int nhits_bottom_physics = 0;
      int nhits_side_physics = 0;
      int nhits_all_physics = 0;
      double total_charge = 0.;
      nhits_OD_cluster_1 = 0;
      npes_OD_cluster_1 = 0;
      nhits_OD_cluster_2 = 0;
      npes_OD_cluster_2 = 0;
      nhits_OD_cluster_1_center_of_mass = 0;
      for(int i=0; i<n_muon_topologies; i++)
	nhits_OD_cluster_1_center_of_mass_by_topology[i] = 0;

      double cm_x = 0.;
      double cm_y = 0.;
      double cm_z = 0.;


      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_OD_tube_id->at(itrigger)).size(); idigitizedhit++){
	// loop on digitized hits in the trigger
	
	tube_id = (digitized_hit_OD_tube_id->at(itrigger)).at(idigitizedhit);
	all_pmts_tree_OD->GetEntry(tube_id - 1);

	charge = (digitized_hit_OD_Q->at(itrigger)).at(idigitizedhit);
	total_charge += charge;
	h_charge.Fill(charge);

	hit_time = (digitized_hit_OD_time->at(itrigger)).at(idigitizedhit)-time_offset;
	h_hit_time.Fill(hit_time);
	h_pathlength.Fill(hit_time*speedlight);

	cm_x += pmt_x_OD;
	cm_y += pmt_y_OD;
	cm_z += pmt_z_OD;

	if( impact_top )
	  h_pathlength_top.Fill(hit_time*speedlight);
	if( impact_side )
	  h_pathlength_side.Fill(hit_time*speedlight);

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

	if( impact ){
	  distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	  h_distance_pmt_impact.Fill(distance_pmt_impact);

	  if( distance_pmt_impact <= cluster_radius_1 ){
	    nhits_OD_cluster_1 ++;
	    npes_OD_cluster_1 += charge;
	  }
	  if( distance_pmt_impact <= cluster_radius_2 ){
	    nhits_OD_cluster_2 ++;
	    npes_OD_cluster_2 += charge;
	  }
	}

	if( pmt_location_OD == 5 ){
	  nhits_top ++;
	  if( impact )
	    nhits_top_physics ++;
	}
	else if( pmt_location_OD == 4 ){
	  nhits_side ++;
	  if( impact )
	    nhits_side_physics ++;
	}
	else if( pmt_location_OD == 3 ){
	  nhits_bottom ++;
	  if( impact )
	    nhits_bottom_physics ++;
	}
	ibin = h_charge_vs_z.FindBin(pmt_z_OD);
	maxcharge = h_charge_vs_z.GetBinContent(ibin);
	if( charge > maxcharge )
	  h_charge_vs_z.SetBinContent(ibin,charge);


	h_digitized_n_photons_ids.Fill(((digitized_hit_OD_photon_ids->at(itrigger)).at(idigitizedhit)).size());
	if( digitized_hit_OD_photon_ids->size() )
	  h_digitized_photons_id0.Fill(((digitized_hit_OD_photon_ids->at(itrigger)).at(idigitizedhit))[0]);

      }

      h_total_charge.Fill(total_charge);


      h_digitized_nhits_top_tubes.Fill(nhits_top);
      h_digitized_nhits_side_tubes.Fill(nhits_side);
      h_digitized_nhits_bottom_tubes.Fill(nhits_bottom);

      nhits_all = nhits_top + nhits_side + nhits_bottom;
    
      h_digitized_nhits_all_tubes.Fill(nhits_all);
      h_npes_vs_digitized_nhits.Fill(nhits_all, total_charge);
      h_nhits_vs_muon_energy.Fill(muon_energy,nhits_all);

      double phi_vertex = atan2(vtx_y,vtx_x);
      double theta_vertex = atan2(sqrt(pow(vtx_x,2) + pow(vtx_y,2)),vtx_z);
      double phi_cm = atan2(cm_y,cm_x);
      double theta_cm = atan2(sqrt(pow(cm_x,2) + pow(cm_y,2)),cm_z);
      if( OD_impact ){


	h_muon_phi_vertex->Fill(180./pi*phi_vertex);
	h_muon_theta_vertex->Fill(180./pi*theta_vertex);

	h_muon_center_of_mass_phi->Fill(180./pi*phi_cm);
	h_muon_center_of_mass_theta->Fill(180./pi*theta_cm);
      
	h_muon_center_of_mass_vertex_phi_residual->Fill(180./pi*(phi_cm - phi_vertex));
	h_muon_center_of_mass_vertex_theta_residual->Fill(180./pi*(theta_cm - theta_vertex));
      
	h_muon_center_of_mass_vertex_phi_scatter->Fill(180./pi*phi_vertex,180./pi*phi_cm);
	h_muon_center_of_mass_vertex_theta_scatter->Fill(180./pi*theta_vertex,180./pi*theta_cm);

      }

      if( impact ){
	h_digitized_nhits_top_physics_tubes.Fill(nhits_top_physics);
	h_digitized_nhits_side_physics_tubes.Fill(nhits_side_physics);
	h_digitized_nhits_bottom_physics_tubes.Fill(nhits_bottom_physics);

	nhits_all_physics = nhits_top_physics + nhits_side_physics + nhits_bottom_physics;
	h_digitized_nhits_physics_tubes.Fill(nhits_all_physics);
	if( nhits_all_physics < n_hits_limit ){
	  h_digitized_nhits_physics_tubes_impact_few_nhits.Fill(nhits_all_physics);
	}
	else{
	  h_digitized_nhits_physics_tubes_impact_many_nhits.Fill(nhits_all_physics);
	}

	h_muon_center_of_mass_phi_impact_residual->Fill(180./pi*(phi_cm - atan2(impact_y,impact_x)));
	h_muon_center_of_mass_theta_impact_residual->Fill(180./pi*(theta_cm - atan2(sqrt(pow(impact_x,2) + pow(impact_y,2)),impact_z)));

	h_muon_phi_impact_scatter->Fill(180./pi*atan2(impact_y,impact_x),180./pi*phi_cm);
	h_muon_theta_impact_scatter->Fill(180./pi*atan2(sqrt(pow(impact_x,2) + pow(impact_y,2)),impact_z),180./pi*theta_cm);

	if( nhits_all_physics < n_hits_limit ){
	  h_muon_center_of_mass_vertex_theta_scatter_few_nhits->Fill(180./pi*theta_vertex,180./pi*theta_cm);
	  h_muon_center_of_mass_vertex_theta_residual_few_nhits->Fill(180./pi*(theta_cm - theta_vertex));
	  h_muon_center_of_mass_vertex_phi_residual_few_nhits->Fill(180./pi*(phi_cm - phi_vertex));
	}
	else{
	  h_muon_center_of_mass_vertex_theta_scatter_many_nhits->Fill(180./pi*theta_vertex,180./pi*theta_cm);
	  h_muon_center_of_mass_vertex_theta_residual_many_nhits->Fill(180./pi*(theta_cm - theta_vertex));
	  h_muon_center_of_mass_vertex_phi_residual_many_nhits->Fill(180./pi*(phi_cm - phi_vertex));
	}
	if( nhits_OD_cluster_1 == 0 ){
	  std::cout << " event " << ievent << " of " << primary_events_tree->GetEntries() << std::endl;
	  std::cout << " impact " << impact << " impact_top " << impact_top << " impact_side " << impact_side << " impact_bottom " << impact_bottom << std::endl;
	  std::cout << " impact point( " << impact_x << ", " << impact_y << ", " << impact_z << ") impact radius " << sqrt(pow(impact_x,2) + pow(impact_y,2)) << std::endl;
	  std::cout << " vtx point( " << vtx_x << ", " << vtx_y << ", " << vtx_z << ") vtx_radius " << vtx_radius << std::endl;
	  std::cout << " dir point( " << dir_x << ", " << dir_y << ", " << dir_z << ") " << std::endl;
	  std::cout << " nhits_top " << nhits_top << " nhits_side " << nhits_side << " nhits_bottom " << nhits_bottom << " nhits_all_physics " << nhits_all_physics <<  std::endl;
	  std::cout << " path length " << sqrt(pow(vtx_x-impact_x,2) + pow(vtx_y-impact_y,2) + pow(vtx_z-impact_z,2)) << " E " << muon_energy << " possible length " << muon_energy/200.e-2 << std::endl;

	}

	h_nhits_OD_cluster_1.Fill(nhits_OD_cluster_1);
	h_npes_OD_cluster_1.Fill(npes_OD_cluster_1);
	h_nhits_OD_cluster_2.Fill(nhits_OD_cluster_2);
	h_npes_OD_cluster_2.Fill(npes_OD_cluster_2);
      
	for(int itrack=0; itrack<trigger_ntrack->at(itrigger); itrack++){
	  // loop on tracks in the event
	  if( track_ipnu->at(itrigger).at(itrack) == muon_pdg_id && track_M->at(itrigger).at(itrack) > 0 ){ // muon
	    double muon_path_length = sqrt(pow(track_start_x->at(itrigger).at(itrack) - impact_x,2) + pow(track_start_y->at(itrigger).at(itrack) - impact_y,2) + pow(track_start_z->at(itrigger).at(itrack) - impact_z,2));
	    h_muon_path_length_impact.Fill(muon_path_length);
	    if( nhits_all_physics < n_hits_limit ){
	      h_muon_start_x_y_z_impact_few_nhits->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	      h_muon_start_x_y_impact_few_nhits->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack));
	      h_muon_start_z_impact_few_nhits->Fill(track_start_z->at(itrigger).at(itrack));
	      if( impact_top )
		h_angle_muon_direction_normal_top_impact_few_nhits.Fill(track_uz->at(itrigger).at(itrack));
	      if( impact_side )
		h_angle_muon_direction_normal_side_impact_few_nhits.Fill(
									  (track_ux->at(itrigger).at(itrack)*impact_x + track_uy->at(itrigger).at(itrack)*impact_y)/sqrt(pow(impact_x,2) + pow(impact_y,2))
									  );
	      h_muon_path_length_impact_few_nhits.Fill(muon_path_length);
	    }
	    else{
	      h_muon_start_x_y_z_impact_many_nhits->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack));
	      h_muon_start_x_y_impact_many_nhits->Fill(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack));
	      h_muon_start_z_impact_many_nhits->Fill(track_start_z->at(itrigger).at(itrack));
	      if( impact_top )
		h_angle_muon_direction_normal_top_impact_many_nhits.Fill(track_uz->at(itrigger).at(itrack));
	      if( impact_side )
		h_angle_muon_direction_normal_side_impact_many_nhits.Fill(
									  (track_ux->at(itrigger).at(itrack)*impact_x + track_uy->at(itrigger).at(itrack)*impact_y)/sqrt(pow(impact_x,2) + pow(impact_y,2))
									   );
	      h_muon_path_length_impact_many_nhits.Fill(muon_path_length);
	    }
	  }
	}

	for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_OD_tube_id->at(itrigger)).size(); idigitizedhit++){
	  // loop on digitized hits in the trigger
	  
	  tube_id = (digitized_hit_OD_tube_id->at(itrigger)).at(idigitizedhit);
	  all_pmts_tree_OD->GetEntry(tube_id - 1);

	  distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	  hit_time = (digitized_hit_OD_time->at(itrigger)).at(idigitizedhit)-time_offset;
	  
	  if( nhits_all_physics < n_hits_limit ){
	    h_distance_pmt_impact_few_nhits.Fill(distance_pmt_impact);
	    h_hit_time_few_nhits.Fill(hit_time);
	  }
	  else{
	    h_distance_pmt_impact_many_nhits.Fill(distance_pmt_impact);
	    h_hit_time_many_nhits.Fill(hit_time);
	  }
	}
      }
      if( OD_impact ){
	for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_OD_tube_id->at(itrigger)).size(); idigitizedhit++){
	  // loop on digitized hits in the trigger
	  
	  tube_id = (digitized_hit_OD_tube_id->at(itrigger)).at(idigitizedhit);
	  all_pmts_tree_OD->GetEntry(tube_id - 1);

	  double phi_PMT = atan2(pmt_y_OD,pmt_x_OD);
	  double theta_PMT = atan2(sqrt(pow(pmt_x_OD,2) + pow(pmt_y_OD,2)),pmt_z_OD);
	  distance_pmt_cm = sin(theta_PMT)*sin(theta_cm)*cos(phi_PMT-phi_cm) + cos(theta_PMT)*cos(theta_cm);
	  h_distance_pmt_center_of_mass.Fill(180./pi*acos(distance_pmt_cm));
	  if( 180./pi*acos(distance_pmt_cm) <= cluster_radius_cm ){
	    nhits_OD_cluster_1_center_of_mass ++;
	  }
	}
	h_nhits_OD_cluster_1_center_of_mass.Fill(nhits_OD_cluster_1_center_of_mass);
	if( nhits_all_physics < n_hits_limit )
	  h_nhits_OD_cluster_1_center_of_mass_few_nhits.Fill(nhits_OD_cluster_1_center_of_mass);
	else
	  h_nhits_OD_cluster_1_center_of_mass_many_nhits.Fill(nhits_OD_cluster_1_center_of_mass);
      }
      if( muon_topology_true >= 0 ){
	h_digitized_nhits_all_tubes_by_topology[muon_topology_true]->Fill((digitized_hit_OD_tube_id->at(itrigger)).size());
	h_digitized_nhits_zoom_all_tubes_by_topology[muon_topology_true]->Fill((digitized_hit_OD_tube_id->at(itrigger)).size());
	for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_OD_tube_id->at(itrigger)).size(); idigitizedhit++){
	  // loop on digitized hits in the trigger
	  
	  tube_id = (digitized_hit_OD_tube_id->at(itrigger)).at(idigitizedhit);
	  all_pmts_tree_OD->GetEntry(tube_id - 1);

	  double phi_PMT = atan2(pmt_y_OD,pmt_x_OD);
	  double theta_PMT = atan2(sqrt(pow(pmt_x_OD,2) + pow(pmt_y_OD,2)),pmt_z_OD);
	  distance_pmt_cm = sin(theta_PMT)*sin(theta_cm)*cos(phi_PMT-phi_cm) + cos(theta_PMT)*cos(theta_cm);
	  h_distance_pmt_center_of_mass_by_topology[muon_topology_true]->Fill(180./pi*acos(distance_pmt_cm));
	  if( 180./pi*acos(distance_pmt_cm) <= cluster_radius_cm ){
	    nhits_OD_cluster_1_center_of_mass_by_topology[muon_topology_true] ++;
	  }
	}
	h_nhits_OD_cluster_1_center_of_mass_by_topology[muon_topology_true]->Fill(nhits_OD_cluster_1_center_of_mass_by_topology[muon_topology_true]);
	h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[muon_topology_true]->Fill((digitized_hit_OD_tube_id->at(itrigger)).size() - nhits_OD_cluster_1_center_of_mass_by_topology[muon_topology_true]);

      }

    }
  }


  f->Close();



  TH1F h_dynamic_range("h_dynamic_range","h_dynamic_range; dynamic range [p.e.]; fraction of saturated electronics",h_digitized_hit_OD_Q_side_tubes.GetNbinsX(),h_digitized_hit_OD_Q_side_tubes.GetXaxis()->GetXmin(),h_digitized_hit_OD_Q_side_tubes.GetXaxis()->GetXmax());
  double charge_integral = h_digitized_hit_OD_Q_side_tubes.Integral();
  for(int i=1; i<=h_digitized_hit_OD_Q_side_tubes.GetNbinsX(); i++){
    double local_integral = h_digitized_hit_OD_Q_side_tubes.Integral(i,h_digitized_hit_OD_Q_side_tubes.GetNbinsX());
    double charge_ratio = local_integral/charge_integral;
    h_dynamic_range.SetBinContent(i,charge_ratio);
  }

  TH1F h_nhits_cluster_1_efficiency("h_nhits_cluster_1_efficiency",Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1),h_nhits_OD_cluster_1.GetNbinsX(),h_nhits_OD_cluster_1.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1.GetXaxis()->GetXmax());
  h_nhits_cluster_1_efficiency.SetLineColor(h_nhits_OD_cluster_1.GetLineColor());
  h_nhits_cluster_1_efficiency.SetLineWidth(2);
  double nhits_integral = h_nhits_OD_cluster_1.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1.GetNbinsX(); i++){
    double local_integral = h_nhits_OD_cluster_1.Integral(i,h_nhits_OD_cluster_1.GetNbinsX());
    double nhits_ratio = local_integral/nhits_integral;
    h_nhits_cluster_1_efficiency.SetBinContent(i,nhits_ratio);
  }


  TH1F h_nhits_cluster_2_efficiency("h_nhits_cluster_2_efficiency",Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_2),h_nhits_OD_cluster_2.GetNbinsX(),h_nhits_OD_cluster_2.GetXaxis()->GetXmin(),h_nhits_OD_cluster_2.GetXaxis()->GetXmax());
  h_nhits_cluster_2_efficiency.SetLineColor(h_nhits_OD_cluster_2.GetLineColor());
  h_nhits_cluster_2_efficiency.SetLineWidth(2);
  nhits_integral = h_nhits_OD_cluster_2.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_2.GetNbinsX(); i++){
    double local_integral = h_nhits_OD_cluster_2.Integral(i,h_nhits_OD_cluster_2.GetNbinsX());
    double nhits_ratio = local_integral/nhits_integral;
    h_nhits_cluster_2_efficiency.SetBinContent(i,nhits_ratio);
  }


  TH1F h_npes_cluster_1_efficiency("h_npes_cluster_1_efficiency",Form("cluster (< %.0f); threshold [npes]; efficiency",cluster_radius_1),h_npes_OD_cluster_1.GetNbinsX(),h_npes_OD_cluster_1.GetXaxis()->GetXmin(),h_npes_OD_cluster_1.GetXaxis()->GetXmax());
  h_npes_cluster_1_efficiency.SetLineColor(h_npes_OD_cluster_1.GetLineColor());
  h_npes_cluster_1_efficiency.SetLineWidth(2);
  double npes_integral = h_npes_OD_cluster_1.Integral();
  for(int i=1; i<=h_npes_OD_cluster_1.GetNbinsX(); i++){
    double local_integral = h_npes_OD_cluster_1.Integral(i,h_npes_OD_cluster_1.GetNbinsX());
    double npes_ratio = local_integral/npes_integral;
    h_npes_cluster_1_efficiency.SetBinContent(i,npes_ratio);
  }


  TH1F h_npes_cluster_2_efficiency("h_npes_cluster_2_efficiency",Form("cluster (< %.0f); threshold [npes]; efficiency",cluster_radius_2),h_npes_OD_cluster_2.GetNbinsX(),h_npes_OD_cluster_2.GetXaxis()->GetXmin(),h_npes_OD_cluster_2.GetXaxis()->GetXmax());
  h_npes_cluster_2_efficiency.SetLineColor(h_npes_OD_cluster_2.GetLineColor());
  h_npes_cluster_2_efficiency.SetLineWidth(2);
  npes_integral = h_npes_OD_cluster_2.Integral();
  for(int i=1; i<=h_npes_OD_cluster_2.GetNbinsX(); i++){
    double local_integral = h_npes_OD_cluster_2.Integral(i,h_npes_OD_cluster_2.GetNbinsX());
    double npes_ratio = local_integral/npes_integral;
    h_npes_cluster_2_efficiency.SetBinContent(i,npes_ratio);
  }


  TH1F h_nhits_cluster_1_center_of_mass_efficiency("h_nhits_cluster_1_center_of_mass_efficiency",Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1),h_nhits_OD_cluster_1_center_of_mass.GetNbinsX(),h_nhits_OD_cluster_1_center_of_mass.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_center_of_mass.GetXaxis()->GetXmax());
  h_nhits_cluster_1_center_of_mass_efficiency.SetLineColor(h_nhits_OD_cluster_1_center_of_mass.GetLineColor());
  h_nhits_cluster_1_center_of_mass_efficiency.SetLineWidth(2);
  double nhits_integral_center_of_mass = h_nhits_OD_cluster_1_center_of_mass.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_center_of_mass.GetNbinsX(); i++){
    double local_integral = h_nhits_OD_cluster_1_center_of_mass.Integral(i,h_nhits_OD_cluster_1_center_of_mass.GetNbinsX());
    double nhits_ratio = local_integral/nhits_integral_center_of_mass;
    h_nhits_cluster_1_center_of_mass_efficiency.SetBinContent(i,nhits_ratio);
  }


  TH1F h_nhits_cluster_1_center_of_mass_efficiency_few_nhits("h_nhits_cluster_1_center_of_mass_efficiency_few_nhits",Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1),h_nhits_OD_cluster_1_center_of_mass_few_nhits.GetNbinsX(),h_nhits_OD_cluster_1_center_of_mass_few_nhits.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_center_of_mass_few_nhits.GetXaxis()->GetXmax());
  h_nhits_cluster_1_center_of_mass_efficiency_few_nhits.SetLineColor(h_nhits_OD_cluster_1_center_of_mass_few_nhits.GetLineColor());
  h_nhits_cluster_1_center_of_mass_efficiency_few_nhits.SetLineWidth(2);
  double nhits_integral_center_of_mass_few_nhits = h_nhits_OD_cluster_1_center_of_mass_few_nhits.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_center_of_mass_few_nhits.GetNbinsX(); i++){
    double local_integral = h_nhits_OD_cluster_1_center_of_mass_few_nhits.Integral(i,h_nhits_OD_cluster_1_center_of_mass_few_nhits.GetNbinsX());
    double nhits_ratio = local_integral/nhits_integral_center_of_mass_few_nhits;
    h_nhits_cluster_1_center_of_mass_efficiency_few_nhits.SetBinContent(i,nhits_ratio);
  }


  TH1F h_nhits_cluster_1_center_of_mass_efficiency_many_nhits("h_nhits_cluster_1_center_of_mass_efficiency_many_nhits",Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1),h_nhits_OD_cluster_1_center_of_mass_many_nhits.GetNbinsX(),h_nhits_OD_cluster_1_center_of_mass_many_nhits.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_center_of_mass_many_nhits.GetXaxis()->GetXmax());
  h_nhits_cluster_1_center_of_mass_efficiency_many_nhits.SetLineColor(h_nhits_OD_cluster_1_center_of_mass_many_nhits.GetLineColor());
  h_nhits_cluster_1_center_of_mass_efficiency_many_nhits.SetLineWidth(2);
  double nhits_integral_center_of_mass_many_nhits = h_nhits_OD_cluster_1_center_of_mass_many_nhits.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_center_of_mass_many_nhits.GetNbinsX(); i++){
    double local_integral = h_nhits_OD_cluster_1_center_of_mass_many_nhits.Integral(i,h_nhits_OD_cluster_1_center_of_mass_many_nhits.GetNbinsX());
    double nhits_ratio = local_integral/nhits_integral_center_of_mass_many_nhits;
    h_nhits_cluster_1_center_of_mass_efficiency_many_nhits.SetBinContent(i,nhits_ratio);
  }


  TH1F * h_nhits_cluster_1_center_of_mass_by_topology_efficiency[n_muon_topologies];
  for(int j=0; j<n_muon_topologies; j++){
    h_nhits_cluster_1_center_of_mass_by_topology_efficiency[j] = new TH1F(Form("h_nhits_cluster_1_center_of_mass_topology_%d_efficiency",j),Form("cluster (< %.0f) %s; threshold [nhits]; efficiency",cluster_radius_cm,topology_name[j].c_str()),h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->GetNbinsX(),h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->GetXaxis()->GetXmax());
    h_nhits_cluster_1_center_of_mass_by_topology_efficiency[j]->SetLineColor(h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->GetLineColor());
    h_nhits_cluster_1_center_of_mass_by_topology_efficiency[j]->SetLineWidth(3);
    double nhits_integral_center_of_mass = h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->Integral();
    if( nhits_integral_center_of_mass )
      for(int i=1; i<=h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->GetNbinsX(); i++){
	double local_integral = h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->Integral(i,h_nhits_OD_cluster_1_center_of_mass_by_topology[j]->GetNbinsX());
	double nhits_ratio = local_integral/nhits_integral_center_of_mass;
	h_nhits_cluster_1_center_of_mass_by_topology_efficiency[j]->SetBinContent(i,nhits_ratio);
      }
  }


  TH1F * h_nhits_out_of_cluster_1_center_of_mass_by_topology_efficiency[n_muon_topologies];
  for(int j=0; j<n_muon_topologies; j++){
    h_nhits_out_of_cluster_1_center_of_mass_by_topology_efficiency[j] = new TH1F(Form("h_nhits_out_of_cluster_1_center_of_mass_topology_%d_efficiency",j),Form("out of cluster (< %.0f) %s; threshold [nhits]; efficiency",cluster_radius_cm,topology_name[j].c_str()),h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->GetNbinsX(),h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->GetXaxis()->GetXmin(),h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->GetXaxis()->GetXmax());
    h_nhits_out_of_cluster_1_center_of_mass_by_topology_efficiency[j]->SetLineColor(h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->GetLineColor());
    h_nhits_out_of_cluster_1_center_of_mass_by_topology_efficiency[j]->SetLineWidth(2);
    double nhits_integral_center_of_mass = h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->Integral();
    if( nhits_integral_center_of_mass )
      for(int i=1; i<=h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->GetNbinsX(); i++){
	double local_integral = h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->Integral(i,h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[j]->GetNbinsX());
	double nhits_ratio = local_integral/nhits_integral_center_of_mass;
	h_nhits_out_of_cluster_1_center_of_mass_by_topology_efficiency[j]->SetBinContent(i,nhits_ratio);
      }
  }

  TH1F * h_digitized_nhits_all_tubes_by_topology_efficiency[n_muon_topologies];
  for(int j=0; j<n_muon_topologies; j++){
    h_digitized_nhits_all_tubes_by_topology_efficiency[j] = new TH1F(Form("h_digitized_nhits_all_tubes_by_topology_%d_efficiency",j),Form("%s; threshold [nhits]; efficiency",topology_name[j].c_str()),h_digitized_nhits_all_tubes_by_topology[j]->GetNbinsX(),h_digitized_nhits_all_tubes_by_topology[j]->GetXaxis()->GetXmin(),h_digitized_nhits_all_tubes_by_topology[j]->GetXaxis()->GetXmax());
    h_digitized_nhits_all_tubes_by_topology_efficiency[j]->SetLineColor(h_digitized_nhits_all_tubes_by_topology[j]->GetLineColor());
    h_digitized_nhits_all_tubes_by_topology_efficiency[j]->SetLineWidth(3);
    double nhits_integral_center_of_mass = h_digitized_nhits_all_tubes_by_topology[j]->Integral();
    if( nhits_integral_center_of_mass )
      for(int i=1; i<=h_digitized_nhits_all_tubes_by_topology[j]->GetNbinsX(); i++){
	double local_integral = h_digitized_nhits_all_tubes_by_topology[j]->Integral(i,h_digitized_nhits_all_tubes_by_topology[j]->GetNbinsX());
	double nhits_ratio = local_integral/nhits_integral_center_of_mass;
	h_digitized_nhits_all_tubes_by_topology_efficiency[j]->SetBinContent(i,nhits_ratio);
      }
  }


  TH1F * h_digitized_nhits_zoom_all_tubes_by_topology_efficiency[n_muon_topologies];
  for(int j=0; j<n_muon_topologies; j++){
    h_digitized_nhits_zoom_all_tubes_by_topology_efficiency[j] = new TH1F(Form("h_digitized_nhits_zoom_all_tubes_by_topology_%d_efficiency",j),Form("%s; threshold [nhits]; efficiency",topology_name[j].c_str()),h_digitized_nhits_zoom_all_tubes_by_topology[j]->GetNbinsX(),h_digitized_nhits_zoom_all_tubes_by_topology[j]->GetXaxis()->GetXmin(),h_digitized_nhits_zoom_all_tubes_by_topology[j]->GetXaxis()->GetXmax());
    h_digitized_nhits_zoom_all_tubes_by_topology_efficiency[j]->SetLineColor(h_digitized_nhits_zoom_all_tubes_by_topology[j]->GetLineColor());
    h_digitized_nhits_zoom_all_tubes_by_topology_efficiency[j]->SetLineWidth(3);
    double nhits_integral_center_of_mass = h_digitized_nhits_zoom_all_tubes_by_topology[j]->Integral();
    if( nhits_integral_center_of_mass )
      for(int i=1; i<=h_digitized_nhits_zoom_all_tubes_by_topology[j]->GetNbinsX(); i++){
	double local_integral = h_digitized_nhits_zoom_all_tubes_by_topology[j]->Integral(i,h_digitized_nhits_zoom_all_tubes_by_topology[j]->GetNbinsX());
	double nhits_ratio = local_integral/nhits_integral_center_of_mass;
	h_digitized_nhits_zoom_all_tubes_by_topology_efficiency[j]->SetBinContent(i,nhits_ratio);
      }
  }




  of->cd();
  PMT_x_y_z.Write();
  PMT_OD_x_y_z.Write();
  PMT_OD_x_y.Write();
  PMT_OD_z.Write();
  h_muon_topology_true.Write();
  h_muon_n_intersections_OD_true.Write();
  h_muon_n_intersections_ID_true.Write();
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
  h_hit_time_few_nhits.Write();
  h_hit_time_many_nhits.Write();
  h_pathlength.Write();
  h_pathlength_top.Write();
  h_pathlength_side.Write();
  h_total_charge.Write();
  h_n_hits_vs_radius.Write();
  h_distance_pmt_impact.Write();
  h_distance_pmt_impact_few_nhits.Write();
  h_distance_pmt_impact_many_nhits.Write();
  h_distance_pmt_center_of_mass.Write();
  for(int i=0; i<n_muon_topologies; i++)
    h_distance_pmt_center_of_mass_by_topology[i]->Write();
  for(int i=0; i<n_muon_topologies; i++){
    h_muon_stop_x_y_z_by_topology[i]->Write();
    h_muon_stop_x_y_by_topology[i]->Write();
    h_muon_stop_z_by_topology[i]->Write();
    h_muon_stop_by_tooplogy[i]->Write();
  }
  h_nhits_OD_cluster_1.Write();
  h_npes_OD_cluster_1.Write();
  h_nhits_OD_cluster_2.Write();
  h_npes_OD_cluster_2.Write();
  h_nhits_OD_cluster_1_center_of_mass.Write();
  h_nhits_OD_cluster_1_center_of_mass_few_nhits.Write();
  h_nhits_OD_cluster_1_center_of_mass_many_nhits.Write();
  for(int i=0; i<n_muon_topologies; i++){
    h_nhits_OD_cluster_1_center_of_mass_by_topology[i]->Write();
    h_nhits_OD_out_of_cluster_1_center_of_mass_by_topology[i]->Write();
  }
  h_digitized_n_photons_ids.Write();
  h_digitized_photons_id0.Write();
  h_digitized_nhits_top_tubes.Write();
  h_digitized_nhits_side_tubes.Write();
  h_digitized_nhits_bottom_tubes.Write();
  h_digitized_nhits_all_tubes.Write();
  for(int i=0; i<n_muon_topologies; i++){
    h_digitized_nhits_all_tubes_by_topology[i]->Write();
    h_digitized_nhits_zoom_all_tubes_by_topology[i]->Write();
  }
  h_digitized_nhits_top_physics_tubes.Write();
  h_digitized_nhits_side_physics_tubes.Write();
  h_digitized_nhits_bottom_physics_tubes.Write();
  h_digitized_nhits_physics_tubes.Write();
  h_digitized_nhits_physics_tubes_impact_few_nhits.Write();
  h_digitized_nhits_physics_tubes_impact_many_nhits.Write();
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
  h_muon_start_x_y_impact_few_nhits->Write();
  h_muon_start_x_y_impact_many_nhits->Write();
  h_muon_start_z->Write();
  h_muon_start_z_impact_few_nhits->Write();
  h_muon_start_z_impact_many_nhits->Write();
  h_muon_start_x_y_z_impact->Write();
  h_muon_phi_vertex->Write();
  h_muon_theta_vertex->Write();
  h_muon_phi_impact->Write();
  h_muon_theta_impact->Write();
  h_muon_center_of_mass_phi->Write();
  h_muon_center_of_mass_theta->Write();
  h_muon_center_of_mass_phi_impact_residual->Write();
  h_muon_center_of_mass_theta_impact_residual->Write();
  h_muon_phi_impact_scatter->Write();
  h_muon_theta_impact_scatter->Write();
  h_muon_center_of_mass_vertex_phi_residual->Write();
  h_muon_center_of_mass_vertex_theta_residual->Write();
  h_muon_center_of_mass_vertex_theta_residual_few_nhits->Write();
  h_muon_center_of_mass_vertex_theta_residual_many_nhits->Write();
  h_muon_center_of_mass_vertex_phi_residual_few_nhits->Write();
  h_muon_center_of_mass_vertex_phi_residual_many_nhits->Write();
  h_muon_center_of_mass_vertex_phi_scatter->Write();
  h_muon_center_of_mass_vertex_theta_scatter->Write();
  h_muon_center_of_mass_vertex_theta_scatter_few_nhits->Write();
  h_muon_center_of_mass_vertex_theta_scatter_many_nhits->Write();
  h_muon_start_x_y_z_impact_top->Write();
  h_muon_start_x_y_z_impact_side->Write();
  h_muon_start_x_y_z_impact_bottom->Write();
  h_muon_start_x_y_z_impact_few_nhits->Write();
  h_muon_start_x_y_z_impact_many_nhits->Write();
  h_angle_muon_direction_normal_top_impact_few_nhits.Write();
  h_angle_muon_direction_normal_top_impact_many_nhits.Write();
  h_angle_muon_direction_normal_side_impact_few_nhits.Write();
  h_angle_muon_direction_normal_side_impact_many_nhits.Write();
  h_muon_impact_x_y_z->Write();
  h_muon_impact_vs_start_x->Write();
  h_muon_impact_vs_start_y->Write();
  h_muon_impact_vs_start_z->Write();
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
  h_muon_path_length_impact.Write();
  h_muon_path_length_impact_few_nhits.Write();
  h_muon_path_length_impact_many_nhits.Write();
  h_dynamic_range.Write();
  h_nhits_cluster_1_efficiency.Write();
  h_nhits_cluster_2_efficiency.Write();
  h_npes_cluster_1_efficiency.Write();
  h_npes_cluster_2_efficiency.Write();
  h_nhits_cluster_1_center_of_mass_efficiency.Write();
  h_nhits_cluster_1_center_of_mass_efficiency_few_nhits.Write();
  h_nhits_cluster_1_center_of_mass_efficiency_many_nhits.Write();
  for(int i=0; i<n_muon_topologies; i++){
    h_nhits_cluster_1_center_of_mass_by_topology_efficiency[i]->Write();
    h_nhits_out_of_cluster_1_center_of_mass_by_topology_efficiency[i]->Write();
    h_digitized_nhits_all_tubes_by_topology_efficiency[i]->Write();
    h_digitized_nhits_zoom_all_tubes_by_topology_efficiency[i]->Write();
  }

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


bool is_impact_horizontal_plane(double x, double y, double z,
				double ux, double uy, double uz,
				double R, double zplane, double E,
				double * impact_time, double *impact_x, double *impact_y, double *impact_z){

  bool is = false;

  double time = (zplane - z)/uz; // time of crossing horizontal plane
  double x_plane = x + ux*time; // position when crossing horizontal plane
  double y_plane = y + uy*time;
  double z_plane = z + uz*time;
  //  std::cout << " pos(" << x << ", " << y << ", " << z << ") dir(" << ux << ", " << uy << ", " << uz << ") R " << R << " zplane " << zplane << " time " << time << std::endl;
  *impact_time = time;
  *impact_x = x_plane;
  *impact_y = y_plane;
  *impact_z = z_plane;	
  double path_length = sqrt(pow(x - x_plane,2) + pow(y - y_plane,2) + pow(z - z_plane,2));
  double dEdx = 400.e-2;// MeV/cm, on the safe side
  double possible_length = E/dEdx;
  if( path_length > possible_length )
    return is;
  if( time > 0. ){ // crossing happens after muon is born
    double r_plane = sqrt(pow(x_plane,2) + pow(y_plane,2));
    //    std::cout << " x_plane " << x_plane << " y_plane " << y_plane << " r_plane " << r_plane << std::endl;
    if( r_plane < R ) // crossing point is within cylinder radius
      is = true;
  }
  //  std::cout << " is " << is << std::endl;



  return is;

}

bool is_impact_side(double x, double y, double z,
		    double ux, double uy, double uz,
		    double R, double H, double E,
		    double * impact_time, double *impact_x, double *impact_y, double *impact_z){

  bool is = false;

  double r02 = pow(x,2) + pow(y,2); // initial distance from axis
  double r0 = sqrt(r02); // initial distance from axis
  //  std::cout << " r0 " << r0 << " z " << z << std::endl;
  if( r0 < R ){ // muon starts within cylinder radius
    if( z > -H && z < H ) // weird: muon starts inside tank
      return is;
    // muon starts above tank within cylinder radius
    // can't go through side without going first through top
    return is;
  }


  // line - circle intersection
  double x1 = x + ux;
  double y1 = y + uy;
  double D = x*y1 - x1*y;
  double u2 = pow(ux,2) + pow(uy,2);
  double Delta = pow(R,2)*u2 - pow(D,2);

  if( Delta <= 0 ) return is; // no intersection

  Delta = sqrt(Delta);

  double sign = 1.;
  if( uy < 0 ) sign = -1.;
  double x_side_1 = (D*uy + sign*ux*Delta)/u2;
  double x_side_2 = (D*uy - sign*ux*Delta)/u2;
  double t1 = (x_side_1 - x)/ux;


  double y_side_1 = (-D*ux + fabs(uy)*Delta)/u2;
  double y_side_2 = (-D*ux - fabs(uy)*Delta)/u2;
  double t2 = (x_side_2 - x)/ux;

  double r1 = sqrt(pow(x_side_1,2) + pow(y_side_1,2));
  double r2 = sqrt(pow(x_side_2,2) + pow(y_side_2,2));
  //  std::clog << " r1  " << r1 << " t1 " << t1 << " r2 " << r2 << " t2 " << t2 << std::endl;

  double time_intersection;
  if( t1 < 0. && t2 < 0. ) // both intersections before muon is born
    return is;

  bool use_1 = false;
  if( t1 < 0. )
    use_1 = false;
  else if( t2 < 0. )
    use_1 = true;
  else{
    if( t1 < t2 ) use_1 = true;
    else use_1 = false;
  }

  double x_side, y_side, z_side, r_side;

  if( use_1 ){
    time_intersection = t1;
    x_side = x_side_1;
    y_side = y_side_1;
  }else{
    time_intersection = t2;
    x_side = x_side_2;
    y_side = y_side_2;
  }
  z_side = z + uz*time_intersection;
    
  r_side = sqrt(pow(x_side,2) + pow(y_side,2)); // distance from cylinder axis

  //  std::cout << " time_intersection " << time_intersection << " x_side " << x_side << " y_side " << y_side << " z_side " << z_side << " r_side " << r_side << std::endl;
  if( z_side < -H )
    return is;
  if( z_side > H )
    return is;
  
  double path_length = sqrt(pow(x - x_side,2) + pow(y - y_side,2) + pow(z - z_side,2));
  double dEdx = 400.e-2;// MeV/cm, on the safe side
  double possible_length = E/dEdx;
  if( path_length > possible_length )
    return is;

  is = true;
  //  std::cout << " is " << is << std::endl;

  *impact_time = time_intersection;
  *impact_x = x_side;
  *impact_y = y_side;
  *impact_z = z_side;

  // double xv = x*ux + y*uy;
  // double time = - xv/u2; // time of maximum approach to cylinder axis
  // double x_side = x + ux*time; // position at maximum approach to cylinder axis
  // double y_side = y + uy*time;
  // double z_side = z + uz*time;
  // double r_side2 = pow(x_side,2) + pow(y_side,2);
  // double r_side = sqrt(r_side2); // min distance from cylinder axis
  // std::cout << " pos(" << x << ", " << y << ", " << z << ") dir(" << ux << ", " << uy << ", " << uz << ") R " << R << " H " << H << " time " << time << " r_side " << r_side << " z_side " << z_side << std::endl;

  // if( r_side > R ) // muon never gets within cylinder radius
  //   return is;

  // // muon starts outside cylinder radius, but does get inside cylinder radius

  // double time_intersection_1 = - xv - sqrt(pow(xv,2) - u2*(r02 - pow(R,2))); // times at which muon intersects cylinder radius
  // double time_intersection_2 = - xv + sqrt(pow(xv,2) - u2*(r02 - pow(R,2)));

  // double time_intersection;
  // std::cout << " time_intersection_1 " << time_intersection_1 << " time_intersection_2 " << time_intersection_2 << std::endl;
  // if( time + time_intersection_1 < 0. ){
  //   if( time + time_intersection_2 < 0. ) // both intersections before muon is born
  //     return is;
  //   time_intersection = time_intersection_2;
  // }else{
  //   time_intersection =time_intersection_1;
  // }

  // x_side += ux*time_intersection; // position when intersecting radius
  // y_side += uy*time_intersection;
  // z_side += uz*time_intersection;
  // r_side = sqrt(pow(x_side,2) + pow(y_side,2)); // distance from cylinder axis

  // std::cout << " time_intersection " << time_intersection << " x_side " << x_side << " y_side " << y_side << " z_side " << z_side << " r_side " << r_side << std::endl;
  // if( z_side < -H )
  //   return is;
  // if( z_side > H )
  //   return is;
  

  // is = true;
  // //  std::cout << " is " << is << std::endl;

  // *impact_time = time_intersection;
  // *impact_x = x_side;
  // *impact_y = y_side;
  // *impact_z = z_side;

  return is;

}

bool is_crossing_horizontal_plane(double x0, double y0, double z0,
				  double x1, double y1, double z1,
				  double R, double zplane){
  bool upwards_track = true;
  if( z1 > z0 ){
    // weird: the track is going up instead of down; it's not a cosmic muon?
    upwards_track = false;
  }



  // make sure muon starts above zplane and stops below zplane
  if( upwards_track ){
    if( zplane > z0 ) return false;
    if( zplane < z1 ) return false;
  }else{
    if( zplane > z1 ) return false;
    if( zplane < z0 ) return false;
  }

  //  x(z)  =  x0 + [(x1 - x0)/(z1 - z0)]*(z - z0)
  //  y(z)  =  y0 + [(y1 - y0)/(z1 - z0)]*(z - z0)

  double xplane = x0 + (x1 - x0)/(z1 - z0)*(zplane - z0);
  double yplane = y0 + (y1 - y0)/(z1 - z0)*(zplane - z0);

  double rplane = sqrt(pow(xplane,2) + pow(yplane,2));
  if( rplane < R ) return true;

  return false;

}


void is_crossing_cylindrical_shell_side(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double Rout, double ztop, double zbottom,
					int * n_intersections_out){


  * n_intersections_out = 0;

  double zhigh = z0;
  double zlow = z1;

  if( z1 > z0 ){
    // weird: the track is going up instead of down; it's not a cosmic muon?
    zhigh = z1;
    zlow = z0;
  }

  // make sure muon is not all above cylinder
  if( ztop < zlow ) return;
  // make sure muon is not all below cylinder
  if( zbottom > zhigh ) return;

  double dx = x1 - x0;
  double dy = y1 - y0;
  double dr = sqrt(pow(dx,2) + pow(dy, 2));
  double D = x0*y1 - x1*y0;

  double discriminant = pow(Rout*dr,2) - pow(D,2);
  // no intersection or tangency between muon and outer circle
  if( discriminant <= 0 ) return;

  int n_intersections = 0;

  double xc1 = (D*dy + (dy/fabs(dy))*dx*sqrt(discriminant))/pow(dr,2);
  double yc1 = (-D*dx + fabs(dy)*sqrt(discriminant))/pow(dr,2);
  double zc1 = z0 + (z1 - z0)/(x1 - x0)*(xc1 - x0);
  if( 
     zc1 < ztop 
     && zc1 > zbottom 
     && zc1 < zhigh 
     && zc1 > zlow ){
    n_intersections ++;
  }

  double xc2 = (D*dy - (dy/fabs(dy))*dx*sqrt(discriminant))/pow(dr,2);
  double yc2 = (-D*dx - fabs(dy)*sqrt(discriminant))/pow(dr,2);
  double zc2 = z0 + (z1 - z0)/(x1 - x0)*(xc2 - x0);
  if( 
     zc2 < ztop 
     && zc2 > zbottom 
     && zc2 < zhigh 
     && zc2 > zlow ){
    n_intersections ++;
  }

  *n_intersections_out = n_intersections;


  return;
}

int point_position(double x0, double y0, double z0,
		   double Rout, double zout, double Rin, double zin){

  double r= sqrt(pow(x0,2) + pow(y0,2));

  if( r < Rin && fabs(z0) < zin )
    return 1;

  if( r < Rout && fabs(z0) < zout )
    return 0;

  return 2;

}


