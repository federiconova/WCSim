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
#include <TMarker.h>
#include <TCanvas.h>
#include <TStyle.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector<std::vector<Int_t> >+;
#pragma link C++ class std::vector<std::vector<Float_t> >+;
#pragma link C++ class std::vector<std::vector<std::vector<Int_t> > >+;
#pragma link C++ class std::vector<std::vector<std::vector<Float_t> > >+;
#endif

int main(){

  int event_min = 6;
  int event_max = 7;
  double angle_max = 85.;
  double hit_time_max = 1300;
  int muon_ipnu = 13; // 13: muon; 111: pi0; 11: electron

  // open input file
  TFile *f = new TFile("output.root","READ");
  
  // create output file
  TFile *of = new TFile("event_display.root", "RECREATE");

  Int_t n_digitized_hits;
  Int_t n_digitized_hits_OD;
  Float_t vertex_x;
  Float_t vertex_y;
  Float_t vertex_z;
  Float_t vertex_phi;
  Float_t vertex_ux;
  Float_t vertex_uy;
  Float_t vertex_uz;
  std::vector<Float_t> v_light_ux;
  std::vector<Float_t> v_light_uy;
  std::vector<Float_t> v_light_uz;
  double light_ux;
  double light_uy;
  double light_uz;
  std::vector<Float_t> v_tube_x;
  std::vector<Float_t> v_tube_y;
  std::vector<Float_t> v_tube_z;
  std::vector<Float_t> v_tube_phi;
  double tube_x;
  double tube_y;
  double tube_z;
  double tube_phi;
  std::vector<Float_t> v_tube_ux;
  std::vector<Float_t> v_tube_uy;
  std::vector<Float_t> v_tube_uz;
  double tube_ux;
  double tube_uy;
  double tube_uz;
  std::vector<Float_t> v_cosangle, v_angle;
  double cosangle, angle;

  std::vector<Float_t> v_light_OD_ux;
  std::vector<Float_t> v_light_OD_uy;
  std::vector<Float_t> v_light_OD_uz;
  std::vector<Float_t> v_tube_OD_x;
  std::vector<Float_t> v_tube_OD_y;
  std::vector<Float_t> v_tube_OD_z;
  std::vector<Float_t> v_tube_OD_phi;
  std::vector<Float_t> v_tube_OD_ux;
  std::vector<Float_t> v_tube_OD_uy;
  std::vector<Float_t> v_tube_OD_uz;
  std::vector<Float_t> v_cosangle_OD, v_angle_OD;


  // geometry tree
  TTree * geom_tree = (TTree*)f->Get("geom_tree");
  Float_t detector_length, detector_radius, pmt_radius;
  Int_t number_of_pmts;
  Int_t number_of_pmts_OD;
  geom_tree->SetBranchAddress("detector_length",&detector_length);
  geom_tree->SetBranchAddress("detector_radius",&detector_radius);
  geom_tree->SetBranchAddress("pmt_radius",&pmt_radius);
  geom_tree->SetBranchAddress("number_of_pmts",&number_of_pmts);
  geom_tree->SetBranchAddress("number_of_pmts_OD",&number_of_pmts_OD);
  geom_tree->GetEntry(0);

  double radius = detector_radius;
  double length = detector_length;

  std::clog << " detector radius " << radius << " length " << length << std::endl;


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

  // all pmts OD tree
  TTree * all_pmts_OD_tree = (TTree*)f->Get("all_pmts_OD_tree");
  Int_t pmt_OD_number, pmt_OD_location;
  Float_t pmt_OD_ux, pmt_OD_uy, pmt_OD_uz, pmt_OD_x, pmt_OD_y, pmt_OD_z;
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_number",&pmt_OD_number);
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_location",&pmt_OD_location);
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_ux",&pmt_OD_ux);
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_uy",&pmt_OD_uy);
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_uz",&pmt_OD_uz);
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_x",&pmt_OD_x);
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_y",&pmt_OD_y);
  all_pmts_OD_tree->SetBranchAddress("pmt_OD_z",&pmt_OD_z);
  for(int ipmt = 0; ipmt < number_of_pmts_OD; ipmt++){
    all_pmts_OD_tree->GetEntry(ipmt);
  }


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

  std::vector<std::vector<Int_t> > * digitized_hit_tube_id = 0;
  std::vector<std::vector<Float_t> > * digitized_hit_Q = 0; std::vector<std::vector<Float_t> > * digitized_hit_time = 0;
  std::vector<std::vector<std::vector<Int_t> > > * digitized_hit_photon_ids = 0;

  std::vector<std::vector<Int_t> > * raw_hit_OD_tube_id = 0; std::vector<std::vector<Int_t> > * raw_hit_OD_tube_times_indexes = 0; std::vector<std::vector<Int_t> > * raw_hit_OD_tube_pe = 0;
  std::vector<std::vector<std::vector<Float_t> > > * raw_hit_OD_times = 0;
  std::vector<std::vector<std::vector<Int_t> > > * raw_hit_OD_parent_ids = 0;

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
  primary_events_tree->SetBranchAddress("trigger_number_raw_hits_OD",&trigger_number_raw_hits_OD);
  primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits);
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

  primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  primary_events_tree->SetBranchAddress("digitized_hit_Q",&digitized_hit_Q);
  primary_events_tree->SetBranchAddress("digitized_hit_time",&digitized_hit_time);
  primary_events_tree->SetBranchAddress("digitized_hit_photon_ids",&digitized_hit_photon_ids);

  primary_events_tree->SetBranchAddress("raw_hit_OD_tube_id",&raw_hit_OD_tube_id);
  primary_events_tree->SetBranchAddress("raw_hit_OD_tube_times_indexes",&raw_hit_OD_tube_times_indexes);
  primary_events_tree->SetBranchAddress("raw_hit_OD_tube_pe",&raw_hit_OD_tube_pe);
  primary_events_tree->SetBranchAddress("raw_hit_OD_times",&raw_hit_OD_times);
  primary_events_tree->SetBranchAddress("raw_hit_OD_parent_ids",&raw_hit_OD_parent_ids);

  primary_events_tree->SetBranchAddress("digitized_hit_OD_tube_id",&digitized_hit_OD_tube_id);
  primary_events_tree->SetBranchAddress("digitized_hit_OD_Q",&digitized_hit_OD_Q);
  primary_events_tree->SetBranchAddress("digitized_hit_OD_time",&digitized_hit_OD_time);
  primary_events_tree->SetBranchAddress("digitized_hit_OD_photon_ids",&digitized_hit_OD_photon_ids);


  double pi = acos(-1.);


  double light_r;

  int ipnu;
  int tube_id, raw_tube_id;

  bool found_raw_hit_id;
  bool found_vertex;


  int nbins = 100;
  double zero = 1.e-3;

  for(int this_event = event_min; this_event < event_max; this_event++){

    TH2F event_display_phi_z("event_display_phi_z","event_display_phi_z; #phi; z", nbins, -pi, pi, nbins, -length/2., length/2.);
    
    TH2F event_display_phi_z_OD("event_display_phi_z_OD","event_display_phi_z_OD; #phi; z", nbins, -pi, pi, nbins, -length/2., length/2.);
    
    TH2F event_display_top("event_display_top","event_display_top; y; -x", nbins, -radius, radius, nbins, -radius, radius);
    
    TH2F event_display_top_OD("event_display_top_OD","event_display_top_OD; y; -x", nbins, -radius, radius, nbins, -radius, radius);
    
    TH2F event_display_bottom("event_display_bottom","event_display_bottom; y; x", nbins, -radius, radius, nbins, -radius, radius);
    
    TH2F event_display_bottom_OD("event_display_bottom_OD","event_display_bottom_OD; y; x", nbins, -radius, radius, nbins, -radius, radius);
    
    for(int i=1; i<=nbins; i++){
      for(int j=1; j<=nbins; j++){
	event_display_phi_z.SetBinContent(i,j,zero);
	event_display_phi_z_OD.SetBinContent(i,j,zero);
	if( std::pow(event_display_top.GetXaxis()->GetBinCenter(i),2) + std::pow(event_display_top.GetYaxis()->GetBinCenter(j),2) <= std::pow(radius,2) ){
	  event_display_top.SetBinContent(i,j,zero);
	  event_display_top_OD.SetBinContent(i,j,zero);
	  event_display_bottom.SetBinContent(i,j,zero);
	  event_display_bottom_OD.SetBinContent(i,j,zero);
	}
      }
    }

  double hit_time;
  int hit_parent_id, primary_track;
  // loop on primary events
  primary_events_tree->GetEvent(this_event); 
  v_light_ux.clear();
  v_light_uy.clear();
  v_light_uz.clear();
  v_tube_x.clear();
  v_tube_y.clear();
  v_tube_z.clear();
  v_tube_phi.clear();
  v_tube_ux.clear();
  v_tube_uy.clear();
  v_tube_uz.clear();
  v_cosangle.clear();
  v_angle.clear();
  
  for(size_t itrigger=0; itrigger<trigger_ntrack->size(); itrigger++){
    // loop on triggers in the event
    
    n_digitized_hits = trigger_number_digitized_hits->at(itrigger);
    n_digitized_hits_OD = trigger_number_digitized_hits_OD->at(itrigger);
    //      h_hits.Fill(n_digitized_hits);
    
    std::clog << " event " << this_event << " nhits " << n_digitized_hits << " nhits_OD " << n_digitized_hits_OD << std::endl;
    
    
    //  reconstruct vertex 
    found_vertex = false;
    for(int itrack=0; itrack<trigger_ntrack->at(itrigger); itrack++){
      // loop on tracks in the trigger
      if( itrack == 0 ) continue;
      ipnu = (track_ipnu->at(itrigger)).at(itrack);
      if( ipnu == muon_ipnu ){
	primary_track = itrack;
	vertex_x = (track_start_x->at(itrigger)).at(itrack);
	vertex_y = (track_start_y->at(itrigger)).at(itrack);
	vertex_z = (track_start_z->at(itrigger)).at(itrack);
	vertex_phi = atan2(vertex_y, vertex_x);
	vertex_ux = (track_ux->at(itrigger)).at(itrack);
	vertex_uy = (track_uy->at(itrigger)).at(itrack);
	vertex_uz = (track_uz->at(itrigger)).at(itrack);
	found_vertex = true;

	std::clog << " itrack " << itrack << " ipnu " << ipnu << " mipnu " << muon_ipnu << " vertex_x " << vertex_x << " vertex_y " << vertex_y << " vertex_z " << vertex_z << " foundv " << found_vertex << " track_parent_type " << (track_parent_type->at(itrigger)).at(itrack) << " track_id " << (track_id->at(itrigger)).at(itrack) << " track_E " << (track_E->at(itrigger)).at(itrack) <<  std::endl;

	//	break;
      }
    }
    if( !found_vertex ) continue;
  
    for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_tube_id->at(itrigger)).size(); idigitizedhit++){
      // loop on digitized hits in the trigger
      
      // get pmt id of digitized hit
      tube_id = (digitized_hit_tube_id->at(itrigger)).at(idigitizedhit);

      // find raw hit associated to digitized hit
      found_raw_hit_id = false;
      for(int irawhit=0; irawhit<trigger_number_raw_hits->at(itrigger); irawhit++){

	// loop on raw hits in the trigger
	raw_tube_id = (raw_hit_tube_id->at(itrigger)).at(irawhit);
	if( raw_tube_id == tube_id ){
	  found_raw_hit_id = true;
	  break;
	}
      }
      if( !found_raw_hit_id ){
	std::clog << " problem: digi tube " << tube_id << " has no raw hit " << std::endl;
	continue;
      }

      
      // load pmt information
      all_pmts_tree->GetEntry(tube_id - 1);
      tube_x = pmt_x;
      tube_y = pmt_y;
      tube_z = pmt_z;
      tube_phi = atan2(pmt_y, pmt_x);
      v_tube_x.push_back(tube_x);
      v_tube_y.push_back(tube_y);
      v_tube_z.push_back(tube_z);
      v_tube_phi.push_back(tube_phi);
      
      // get normal direction into the tube
      tube_ux = - pmt_ux;
      tube_uy = - pmt_uy;
      tube_uz = - pmt_uz;
      // h_tube_ux.Fill(tube_ux);
      // h_tube_uy.Fill(tube_uy);
      // h_tube_uz.Fill(tube_uz);
      v_tube_ux.push_back(tube_ux);
      v_tube_uy.push_back(tube_uy);
      v_tube_uz.push_back(tube_uz);
      
      // retrieve information about light
      light_ux = pmt_x - vertex_x;
      light_uy = pmt_y - vertex_y;
      light_uz = pmt_z - vertex_z;
      light_r = sqrt(pow(light_ux,2) + pow(light_uy,2) + pow(light_uz,2));
      light_ux /= light_r;
      light_uy /= light_r;
      light_uz /= light_r;
      // h_light_ux.Fill(light_ux);
      // h_light_uy.Fill(light_uy);
      // h_light_uz.Fill(light_uz);
      v_light_ux.push_back(light_ux);
      v_light_uy.push_back(light_uy);
      v_light_uz.push_back(light_uz);
      
      cosangle = light_ux*tube_ux + light_uy*tube_uy + light_uz*tube_uz;
      angle = acos(cosangle)*180/pi;
      //	h_cosangle.Fill(cosangle);
      //	h_angle.Fill(angle);
      v_cosangle.push_back(cosangle);
      v_angle.push_back(angle);

      if( fabs(tube_uz) < 0.1 ){
	event_display_phi_z.Fill(tube_phi, tube_z);
      }

      if( tube_uz > 0.1 ){
	event_display_top.Fill(tube_y, -tube_x);
      }

      if( tube_uz < - 0.1 ){
	event_display_bottom.Fill(tube_y, tube_x);
      }

    }


    //    qqq;
    for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_OD_tube_id->at(itrigger)).size(); idigitizedhit++){
      // loop on digitized hits in the trigger
      
      // get pmt id of digitized hit
      tube_id = (digitized_hit_OD_tube_id->at(itrigger)).at(idigitizedhit);

      // find raw hit associated to digitized hit
      found_raw_hit_id = false;
      for(int irawhit=0; irawhit<trigger_number_raw_hits_OD->at(itrigger); irawhit++){

	// loop on raw hits in the trigger
	raw_tube_id = (raw_hit_OD_tube_id->at(itrigger)).at(irawhit);
	if( raw_tube_id == tube_id ){
	  found_raw_hit_id = true;
	  break;
	}
      }
      if( !found_raw_hit_id ){
	std::clog << " problem: digi tube " << tube_id << " has no raw hit " << std::endl;
	continue;
      }

      
      // load pmt information
      all_pmts_OD_tree->GetEntry(tube_id - 1);
      tube_x = pmt_OD_x;
      tube_y = pmt_OD_y;
      tube_z = pmt_OD_z;
      tube_phi = atan2(pmt_OD_y, pmt_OD_x);
      v_tube_OD_x.push_back(tube_x);
      v_tube_OD_y.push_back(tube_y);
      v_tube_OD_z.push_back(tube_z);
      v_tube_OD_phi.push_back(tube_phi);
      
      // get normal direction into the tube
      tube_ux = - pmt_OD_ux;
      tube_uy = - pmt_OD_uy;
      tube_uz = - pmt_OD_uz;
      // h_tube_ux.Fill(tube_ux);
      // h_tube_uy.Fill(tube_uy);
      // h_tube_uz.Fill(tube_uz);
      v_tube_OD_ux.push_back(tube_ux);
      v_tube_OD_uy.push_back(tube_uy);
      v_tube_OD_uz.push_back(tube_uz);
      
      // retrieve information about light
      light_ux = pmt_OD_x - vertex_x;
      light_uy = pmt_OD_y - vertex_y;
      light_uz = pmt_OD_z - vertex_z;
      light_r = sqrt(pow(light_ux,2) + pow(light_uy,2) + pow(light_uz,2));
      light_ux /= light_r;
      light_uy /= light_r;
      light_uz /= light_r;
      // h_light_ux.Fill(light_ux);
      // h_light_uy.Fill(light_uy);
      // h_light_uz.Fill(light_uz);
      v_light_OD_ux.push_back(light_ux);
      v_light_OD_uy.push_back(light_uy);
      v_light_OD_uz.push_back(light_uz);
      
      cosangle = light_ux*tube_ux + light_uy*tube_uy + light_uz*tube_uz;
      angle = acos(cosangle)*180/pi;
      //	h_cosangle.Fill(cosangle);
      //	h_angle.Fill(angle);
      v_cosangle_OD.push_back(cosangle);
      v_angle_OD.push_back(angle);


      if( fabs(tube_uz) < 0.1 ){
	event_display_phi_z_OD.Fill(tube_phi, tube_z);
      }

      if( tube_uz > 0.1 ){
	event_display_top_OD.Fill(tube_y, -tube_x);
      }

      if( tube_uz < - 0.1 ){
	event_display_bottom_OD.Fill(tube_y, tube_x);
      }

    }
    //    qqq;



  }

  std::clog << " vertex_x " << vertex_x << " vertex_y " << vertex_y << " vertex_z " << vertex_z << std::endl;

  TMarker star_phi_z(vertex_phi, vertex_z, 1);
  star_phi_z.SetMarkerStyle(20);
  star_phi_z.SetMarkerColor(kRed);
  star_phi_z.SetMarkerSize(2);

  TMarker star_top(vertex_y, -vertex_x, 1);
  star_top.SetMarkerStyle(20);
  star_top.SetMarkerColor(kRed);
  star_top.SetMarkerSize(2);

  TMarker star_bottom(vertex_y, vertex_x, 1);
  star_bottom.SetMarkerStyle(20);
  star_bottom.SetMarkerColor(kRed);
  star_bottom.SetMarkerSize(2);

  
  of->cd();


  of->Write();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  double BottomMargin = .17;   // fraction of pad height                        
  double LeftMargin = .13;
  TCanvas canvas("canvas","",700,525);
  canvas.SetBottomMargin(BottomMargin);
  canvas.SetLeftMargin(LeftMargin);

  double z_max =  std::max(event_display_phi_z.GetMaximum(),event_display_top.GetMaximum());
  z_max =  std::max(z_max,event_display_bottom.GetMaximum());

  double z_max_OD =  std::max(event_display_phi_z_OD.GetMaximum(),event_display_top_OD.GetMaximum());
  z_max_OD =  std::max(z_max_OD,event_display_bottom_OD.GetMaximum());

  event_display_phi_z.Draw("colz");
  event_display_phi_z.GetZaxis()->SetRangeUser(0.,z_max);
  star_phi_z.Draw("same");
  //  canvas.SetLogz();
  canvas.Print(Form("Event_Display_%d_center_normal.png",this_event));

  event_display_phi_z_OD.Draw("colz");
  event_display_phi_z_OD.GetZaxis()->SetRangeUser(0.,z_max_OD);
  star_phi_z.Draw("same");
  //canvas.SetLogz();
  canvas.Print(Form("Event_Display_%d_center_OD.png",this_event));

  event_display_top.Draw("colz");
  event_display_top.GetZaxis()->SetRangeUser(0.,z_max);
  star_top.Draw("same");
  //  canvas.SetLogz();
  canvas.Print(Form("Event_Display_%d_top_normal.png",this_event));

  event_display_top_OD.Draw("colz");
  event_display_top_OD.GetZaxis()->SetRangeUser(0.,z_max_OD);
  star_top.Draw("same");
  //canvas.SetLogz();
  canvas.Print(Form("Event_Display_%d_top_OD.png",this_event));

  event_display_bottom.Draw("colz");
  event_display_bottom.GetZaxis()->SetRangeUser(0.,z_max);
  star_bottom.Draw("same");
  //  canvas.SetLogz();
  canvas.Print(Form("Event_Display_%d_bottom_normal.png",this_event));

  event_display_bottom_OD.Draw("colz");
  event_display_bottom_OD.GetZaxis()->SetRangeUser(0.,z_max_OD);
  star_bottom.Draw("same");
  //canvas.SetLogz();
  canvas.Print(Form("Event_Display_%d_bottom_OD.png",this_event));

  }

  delete trigger_number ; delete  trigger_date ; delete  trigger_mode ; delete  trigger_vtxvol ; delete  trigger_vec_rec_number ; delete  trigger_jmu ; delete  trigger_jp ; delete  trigger_npar ; delete  trigger_ntrack ; delete  trigger_number_raw_hits ; delete  trigger_number_digitized_hits ; delete  trigger_number_times ;

  delete  trigger_vtx_x ; delete  trigger_vtx_y ; delete  trigger_vtx_z ; delete  trigger_sum_q ;

  delete  track_ipnu ; delete    track_parent_type ; delete    track_flag ; delete    track_start_volume ; delete    track_stop_volume ; delete    track_id ;

  delete  track_ux ; delete    track_uy ; delete    track_uz ; delete    track_M ; delete    track_P ; delete    track_E ; delete    track_px ; delete    track_py ; delete    track_pz ; delete    track_stop_x ; delete    track_stop_y ; delete    track_stop_z ; delete    track_start_x ; delete    track_start_y ; delete    track_start_z ; delete    track_time ;

  delete raw_hit_tube_id; delete raw_hit_tube_times_indexes; delete raw_hit_tube_pe; delete raw_hit_times; delete raw_hit_parent_ids;

  delete digitized_hit_tube_id;  delete digitized_hit_Q;  delete digitized_hit_time;   delete digitized_hit_photon_ids;

  return 1;

}


