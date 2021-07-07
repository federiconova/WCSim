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
#include <TRandom.h>

std::string number_with_digits(int i, int n);
bool is_impact_horizontal_plane(double x, double y, double z,
				double ux, double uy, double uz,
				double R, double zplane, double E,
				double * impact_time, double *impact_x, double *impact_y, double *impact_z);
bool is_impact_side(double x, double y, double z,
		    double ux, double uy, double uz,
		    double R, double H, double E,
		    double * impact_time, double *impact_x, double *impact_y, double *impact_z);

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
  double OD_radius = 3300.;
  double OD_height = 3400.;
  std::clog << " detector_length " << detector_length << " detector_radius " << detector_radius << " pmt_radius " << pmt_radius << " number_of_pmts " << number_of_pmts << " number_of_pmts_OD " << number_of_pmts_OD << std::endl;

  //list of failed pmts IDs
  double failure_rate = 0.05;
  int number_of_failed_pmts = number_of_pmts_OD*failure_rate;
  int* failed_pmts_IDs = (int*)malloc(number_of_failed_pmts*sizeof(int));
  for(int i=0;i<number_of_failed_pmts;i++){
    failed_pmts_IDs[i] = gRandom->Uniform(number_of_pmts_OD);
  }

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


  TH3F PMT_OD_x_y_z("PMT_OD_x_y_z","Topcap; x [cm]; y [cm]; z [cm]",
		    nbins_x,-(r_limit*1.3),(r_limit*1.3),
		    nbins_y,-(r_limit*1.3),(r_limit*1.3),
		    nbins_z,-(z_limit*1.3),z_limit*1.3);
  PMT_OD_x_y_z.SetFillColor(kBlack);
  PMT_OD_x_y_z.SetMarkerColor(kBlack);
  PMT_OD_x_y_z.SetMarkerStyle(1);
  PMT_OD_x_y_z.SetStats(0);

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

  TH3F PMT_OD_off_x_y_z("PMT_OD_off_x_y_z","OD pmt off; x [cm]; y [cm]; z [cm]",
		        nbins_x,-(r_limit*1.3),(r_limit*1.3),
			nbins_y,-(r_limit*1.3),(r_limit*1.3),
			nbins_z,-(z_limit*1.3),z_limit*1.3);
  PMT_OD_off_x_y_z.SetFillColor(kRed);
  PMT_OD_off_x_y_z.SetMarkerColor(kRed);
  PMT_OD_off_x_y_z.SetMarkerStyle(7);
  PMT_OD_off_x_y_z.SetStats(0);

  TH2F PMT_OD_off_x_y("PMT_OD_off_x_y","OD pmt off; x [cm]; y [cm]",
		      nbins_x,-(r_limit*1.3),(r_limit*1.3),
		      nbins_y,-(r_limit*1.3),(r_limit*1.3));
  PMT_OD_off_x_y.SetFillColor(kBlack);
  PMT_OD_off_x_y.SetMarkerColor(kBlack);
  PMT_OD_off_x_y.SetMarkerStyle(1);

  TH1F PMT_OD_off_z("PMT_OD_off_z","OD pmt off; z [cm]",
		    nbins_z,-(z_limit*1.3),(z_limit*1.3));
  PMT_OD_off_z.SetLineColor(kBlack);
  PMT_OD_off_z.SetLineWidth(2);

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    PMT_OD_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
    PMT_OD_x_y.Fill(pmt_x_OD, pmt_y_OD);
    PMT_OD_z.Fill(pmt_z_OD);
    for(int ipmtoff = 0; ipmtoff < number_of_failed_pmts; ipmtoff++){
      if(pmt_number_OD == failed_pmts_IDs[ipmtoff]){
	/*PMT_OD_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	PMT_OD_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	PMT_OD_off_z.Fill(pmt_z_OD);*/
      }
    }
  }
  TH3F PMT_OD_box_off_x_y_z("PMT_OD_box_off_x_y_z","24 Column (8x8 + 8x8(offset 4H)) ; x [cm]; y [cm]; z [cm]",
			    nbins_x,-(r_limit*1.3),(r_limit*1.3),
			    nbins_y,-(r_limit*1.3),(r_limit*1.3),
			    nbins_z,-(z_limit*1.3),z_limit*1.3);
  PMT_OD_box_off_x_y_z.SetFillColor(kRed);
  PMT_OD_box_off_x_y_z.SetMarkerColor(kRed);
  PMT_OD_box_off_x_y_z.SetMarkerStyle(7);
  PMT_OD_box_off_x_y_z.SetStats(0);

  TH2F PMT_OD_box_off_x_y("PMT_OD_box_off_x_y","OD pmt box off; x [cm]; y [cm]",
			  nbins_x,-(r_limit*1.3),(r_limit*1.3),
			  nbins_y,-(r_limit*1.3),(r_limit*1.3));
  PMT_OD_box_off_x_y.SetFillColor(kRed);
  PMT_OD_box_off_x_y.SetMarkerColor(kRed);
  PMT_OD_box_off_x_y.SetMarkerStyle(7);

  TH1F PMT_OD_box_off_z("PMT_OD_box_off_z","OD pmt box off; z [cm]",
			nbins_z,-(z_limit*1.3),(z_limit*1.3));
  PMT_OD_box_off_z.SetLineColor(kRed);
  PMT_OD_box_off_z.SetLineWidth(7);
  
  // 8x8 hole in topcap
  int topcap_failure_block_vertical = 8;
  int topcap_failure_block_horizontal = 8;
  int topcap_failure_block_size = topcap_failure_block_vertical*topcap_failure_block_horizontal;
  int* topcap_pmt_off_id = (int*)malloc(topcap_failure_block_size*sizeof(int));
  int* interleaved_topcap_pmt_off_id = (int*)malloc(0.5*topcap_failure_block_size*sizeof(int));

  int n_topcap_pmts = 1527;
  int topcap_pmt_id[n_topcap_pmts]; //number of topcap pmts
  for(int i = 5617; i< 8672; i++){ //range of topcap pmt IDs
    if(i % 2 !=0) topcap_pmt_id[(i-5617)/2] = i;
    else continue;
  }

  int random_topcap_pmt_id = 6255;
  float pmt_topcap_x_OD_off;
  float pmt_topcap_y_OD_off;

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_uz_OD == 1 && pmt_number_OD == random_topcap_pmt_id){
      pmt_topcap_x_OD_off = pmt_x_OD;
      pmt_topcap_y_OD_off = pmt_y_OD;
    }
  }

  int j = 0;
  int l = 0;
  int r = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_uz_OD ==1){
      if((pmt_x_OD < (pmt_topcap_x_OD_off - 450) || pmt_x_OD > (pmt_topcap_x_OD_off + 600) || pmt_y_OD < (pmt_topcap_y_OD_off - 450) ||  pmt_y_OD > (pmt_topcap_y_OD_off + 600))) continue;
      else{
        topcap_pmt_off_id[j] = pmt_number_OD;
	if(j % 8 ==0 && j!=0 && j!=16 && j!=32 && j!=48) r = 0; //add interleaving to topcap
        else if(j % 8 ==0 && (j==0 || j==16 || j==32 || j==48)) r = 1;
        if((r == 0 && j % 2 ==0) || (r == 1 && j % 2 !=0)){
	  interleaved_topcap_pmt_off_id[l] = topcap_pmt_off_id[j];
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	  PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	  PMT_OD_box_off_z.Fill(pmt_z_OD);
	  */l++;
	}
        j++;
      }
    }
  }

  // 4x36 hole in side
  int failure_block_vertical = 4;
  int failure_block_horizontal = 36;
  int failure_block_size = 0.5*failure_block_vertical*failure_block_horizontal;
  int* low_side_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));
  
  int random_low_pmt_id = 50; //50 for lower region, 2643 for higher region
  float pmt_x_OD_off;
  float pmt_y_OD_off;
  float pmt_z_OD_off;
  
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == random_low_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }  

  int k = 0;
  double z_separation = 145;
  double x_separation = 102.14;//105
  double y_separation = 102.14;//105
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else if(pmt_number_OD < 9000 && pmt_number_OD % 2 != 0) continue; //use for low region
      else if(pmt_number_OD > 9000 && pmt_number_OD % 2 == 0) continue; // use for low region
      else{
	low_side_pmt_off_id[k] = pmt_number_OD;
	/*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	PMT_OD_box_off_z.Fill(pmt_z_OD);
	*/k++;
      }
    }
  }

  //high side 4x36 
  k=0;
  int random_high_side_pmt_id = 2643;
  int* high_side_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));  
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == random_high_side_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else if(pmt_number_OD > 2268 && pmt_number_OD % 2 != 0) continue; //use for high region
      else if(pmt_number_OD <= 2268 && pmt_number_OD % 2 == 0) continue; // use for high region
      else{
        high_side_pmt_off_id[k] = pmt_number_OD;
	/*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	PMT_OD_box_off_z.Fill(pmt_z_OD);
	*/k++;
      }
    }
  }

  // 8x18 hole in side, vertical cabling
  failure_block_vertical = 18;//18
  failure_block_horizontal = 8;//8
  //failure_block_size = 1620;
  int* side_vertical_8x18_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == random_low_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }
  
  int m = 1;//0
  int q = 0;
  int p=0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){ 
      if(pmt_y_OD < pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD > (pmt_z_OD_off + failure_block_vertical*z_separation) || pmt_z_OD < pmt_z_OD_off) continue;
      //else if((m>15 && m<40 && m % 2 ==0)|| (m>39 && m<64 && m % 2 !=0)||(m>63 && m<88 && m % 2 ==0)||(m>87 && m<112 && m % 2 !=0)||(m>111 && m<136 && m % 2 ==0)||(m>135 && m % 2 !=0)|| (m<16 && m!=1 && m!=2 && m!=5 && m!=6 && m!=9 && m!=10 && m!=13 && m!=14)) continue;
      else{ 
	if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){ //switches off odd or even PMTs depending on the row. Will switch every 4th row. 
	  side_vertical_8x18_pmt_off_id[q] = pmt_number_OD;
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	  PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	  PMT_OD_box_off_z.Fill(pmt_z_OD);
	  */q++;
	}
	if((p==0 && m % (2*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){ //factor before failure_block_horizontal may need to change in case p==0. Depends on row of final PMT, extending down the OD (PMTs numbered three rows at a time), will either be 1,2 or 3. 
	  m=0;
	  p++;	  
	}
	m++;
      }
    }
  }  

  //36x4 vertical cabling 
  failure_block_vertical = 36;//18
  failure_block_horizontal = 4;//8
  int* side_vertical_36x4_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));
  int ver_random_pmt_id = 5242;  
  
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }

  p=0;
  m = 1;
  q = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
	if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){
	  side_vertical_36x4_pmt_off_id[q] = pmt_number_OD;
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	  PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	  PMT_OD_box_off_z.Fill(pmt_z_OD);
	  */q++;
	}
	if((p==0 && m % (2*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
	  m=0;
	  p++;
	}
	m++;
      } 
    }
  }

  //overlap of 36x4, 18x8
  int h=0;
  failure_block_vertical = 18;
  failure_block_horizontal = 8;
  int* side_vertical_36x4_overlap_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }
  p = 0;
  m = 1;
  q = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
	if((p % 2 ==0 && pmt_number_OD % 2 != 0) || (p % 2 != 0 && pmt_number_OD % 2 == 0)){
	  side_vertical_36x4_overlap_pmt_off_id[h] = pmt_number_OD;
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	  PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	  PMT_OD_box_off_z.Fill(pmt_z_OD);
	  */h++;
	}
	if((p==0 && m % (2*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
	  m=0;
	  p++;
	}
	m++;
      }
    }
  }
 
  //16x4 vertical cabling, 32 column
  failure_block_vertical = 16;
  failure_block_horizontal = 4;
  z_separation = 143;
  // failure_block_size = (0.5*failure_block_vertical*failure_block_horizontal);
  int* side_vertical_16x4_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }

  p=0;
  m = 1;
  q = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off)continue;
      else{
        if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){
	  side_vertical_16x4_pmt_off_id[q] = pmt_number_OD;
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
          PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
          PMT_OD_box_off_z.Fill(pmt_z_OD);
	  */q++;
        }
        if(m % (3*failure_block_horizontal) == 0) {
          m=0;
          p++;
        }
        m++;
      }
    }
  }
    
  //4x16 vertical cabling, 32 column
  z_separation = 143;
  y_separation = 120;
  failure_block_vertical = 4;
  failure_block_horizontal = 16;
  int* side_vertical_4x16_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }

  p=0;
  m = 1;
  q = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
        if((p % 2 !=0 && pmt_number_OD % 2 == 0) || (p % 2 == 0 && pmt_number_OD % 2 != 0)){
          side_vertical_4x16_pmt_off_id[q] = pmt_number_OD;
          /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
          PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
          PMT_OD_box_off_z.Fill(pmt_z_OD);
      */q++;
        }
	if( m % (3*failure_block_horizontal) == 0){
	  m=0;
	  p++;
	}
	m++;
      }
    }
  }
  
  //4x16 vertical cabling, 32 column cross
  z_separation = 143;  
  y_separation = 120;
  int ver_random_4x16_pmt_id=4392;
  failure_block_vertical = 4;
  failure_block_horizontal = 16;
  int* side_vertical_4x16_cross_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_4x16_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }

  p=0;
  m = 1;
  q = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
	if((p % 2 !=0 && pmt_number_OD % 2 == 0) || (p % 2 == 0 && pmt_number_OD % 2 != 0)){
	  side_vertical_4x16_cross_pmt_off_id[q] = pmt_number_OD;
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
          PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
          PMT_OD_box_off_z.Fill(pmt_z_OD);
          */q++;
        }
        if( m % (3*failure_block_horizontal) == 0){
          m=0;
          p++;
        }
        m++;
      }
    }
  }

  //8x8 vertical cabling, 16 column
  y_separation = 107.40;
  failure_block_vertical = 8;
  failure_block_horizontal = 8;
  int* side_vertical_8x8_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }

  p=0;
  m = 1;
  q = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);  
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
	if((p % 2 !=0 && pmt_number_OD % 2 == 0) || (p % 2 == 0 && pmt_number_OD % 2 != 0)){
	  side_vertical_8x8_pmt_off_id[q] = pmt_number_OD;
	  PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	  PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	  PMT_OD_box_off_z.Fill(pmt_z_OD);
	  q++;
	}
	if((p==0 && m % (failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
	  m=0;
	  p++;
	}
	m++;
      }
    }
  }
  std::cout<<"overlap "<<q<<std::endl;
  
//8x8 vertical cabling, 24 column, overlap
 y_separation = 110;
 failure_block_vertical = 8;
 failure_block_horizontal = 8;
 int ver_random_pmt_8x8_id = 5230;
 int* side_vertical_8x8_offset_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_number_OD == ver_random_pmt_8x8_id){
     pmt_x_OD_off = pmt_x_OD;
     pmt_y_OD_off = pmt_y_OD;
     pmt_z_OD_off = pmt_z_OD;
   }
 }
 p=0;
 m = 1;
 q = 0;
 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
     if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
     else{
       if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){
	 side_vertical_8x8_offset_pmt_off_id[q] = pmt_number_OD;
	 PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	 PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	 PMT_OD_box_off_z.Fill(pmt_z_OD);
	 q++;
       }
       if((p==0 && m % (failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
	 m=0;
	 p++;
       }
       m++;
     }
   }
 }
 std::cout<<"overlap "<<q<<std::endl;


 //4x32 vertical cabling, 32 column 64 readout
 z_separation = 140;
 y_separation = 102;
 failure_block_vertical = 4;
 failure_block_horizontal = 32;
 failure_block_size = 128;
 ver_random_pmt_id = 5278;
 int* side_vertical_4x32_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_number_OD == ver_random_pmt_id){
     pmt_x_OD_off = pmt_x_OD;
     pmt_y_OD_off = pmt_y_OD;
     pmt_z_OD_off = pmt_z_OD;
   }
 }

 p=0;
 m = 1;
 q = 0;
 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
     if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
     else{
       if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){
	 side_vertical_4x32_pmt_off_id[q] = pmt_number_OD;
	 /*	 PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
	 PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
	 PMT_OD_box_off_z.Fill(pmt_z_OD);
	 */q++;
       }
       if(m % (3*failure_block_horizontal) == 0){
	 m=0;
	 p++;
       }
       m++;
     }
   }
 }
 
 //32x4 vertical cabling, 32 column 64 readout
 z_separation = 143;
 y_separation = 40;
 failure_block_vertical = 32;
 failure_block_horizontal = 4;
 int* side_vertical_32x4_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_number_OD == ver_random_pmt_id){
     pmt_x_OD_off = pmt_x_OD;
     pmt_y_OD_off = pmt_y_OD;
     pmt_z_OD_off = pmt_z_OD;
   }
 }

 p=0;
 m = 1;
 q = 0;
 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
     if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
     else{
       if((p % 2 !=0 && pmt_number_OD % 2 == 0) || (p % 2 == 0 && pmt_number_OD % 2 != 0)){
         side_vertical_32x4_pmt_off_id[q] = pmt_number_OD;
         /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
         PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
         PMT_OD_box_off_z.Fill(pmt_z_OD);
         */q++;
       }
       if((p==0 && m % (failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0)){
         m=0;
         p++;
       }
       m++;
     }
   }
 }

 //16x8 vertical cabling, 16 column
 failure_block_vertical = 16;
 failure_block_horizontal = 8;
 y_separation = 50;
 z_separation = 143;
 int* side_vertical_16x8_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_number_OD == ver_random_pmt_id){
     pmt_x_OD_off = pmt_x_OD;
     pmt_y_OD_off = pmt_y_OD;
     pmt_z_OD_off = pmt_z_OD;
   }
 }

 p=0;
 m = 1;
 q = 0;
 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
     if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
     else{
       if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){
         side_vertical_16x8_pmt_off_id[q] = pmt_number_OD;
         /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
         PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
         PMT_OD_box_off_z.Fill(pmt_z_OD);
         */q++;
       }
       if((m % (3*failure_block_horizontal) == 0) ){
         m=0;
         p++;
       }
       m++;
     }
   }
 }
 
 //8x16 vertical cabling, 16 column
 failure_block_vertical = 8;
 failure_block_horizontal = 16;
 y_separation = 70;
 int* side_vertical_8x16_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_number_OD == ver_random_pmt_id){
     pmt_x_OD_off = pmt_x_OD;
     pmt_y_OD_off = pmt_y_OD;
     pmt_z_OD_off = pmt_z_OD;
   }
 }

 p=0;
 m = 1;
 q = 0;
 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
     if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
     else{
       if((p % 2 !=0 && pmt_number_OD % 2 == 0) || (p % 2 == 0 && pmt_number_OD % 2 != 0)){
         side_vertical_8x16_pmt_off_id[q] = pmt_number_OD;
         /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
         PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
         PMT_OD_box_off_z.Fill(pmt_z_OD);
         */q++;
       }
       if((p==0 && m % (failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
         m=0;
         p++;
       }
       m++;
     }
   }
 }

 //16x8 vertical cabling, offset 24 column 64 readout
 failure_block_vertical = 16;
 failure_block_horizontal = 8;
 z_separation = 143;
 y_separation = 70;
 int ver_random_16x8_pmt_id = 5266; //5278 original
 int* side_vertical_16x8_offset_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_number_OD == ver_random_16x8_pmt_id){
     pmt_x_OD_off = pmt_x_OD;
     pmt_y_OD_off = pmt_y_OD;
     pmt_z_OD_off = pmt_z_OD;
   }
 }

 p=0;
 m = 1;
 q = 0;
 for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
   all_pmts_tree_OD->GetEntry(ipmt);
   if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
     if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
     else{
       if((p % 2 !=0 && pmt_number_OD % 2 == 0) || (p % 2 == 0 && pmt_number_OD % 2 != 0)){
         side_vertical_16x8_offset_pmt_off_id[q] = pmt_number_OD;
         /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
         PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
         PMT_OD_box_off_z.Fill(pmt_z_OD);
         */q++;
       }
       if((p==0 && m % (3*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
         m=0;
         p++;
       }
       m++;
     }
   }
 }

  z_separation = 145;
  y_separation = 107.40;
  //12x12 vertical cabling
  ver_random_pmt_id = 5242;
  failure_block_size = 72;
  failure_block_vertical = 12;
  failure_block_horizontal = 12;
  int* side_vertical_12x12_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }
 
  m = 1;
  q = 0;
  p = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
        if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){
          side_vertical_12x12_pmt_off_id[q] = pmt_number_OD;
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
          PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
          PMT_OD_box_off_z.Fill(pmt_z_OD);
  */q++;
        }
        if((p==0 && m % (2*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
	  m=0;
	  p++;
	}
	m++;
      }
    }
  }

  //12x12 overlap vertical cabling, 24x6

  failure_block_vertical = 24;
  failure_block_horizontal = 6;
  int* side_vertical_12x12_overlap_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }

  m = 1;
  q = 0;
  p = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
	if((p % 2 ==0 && pmt_number_OD % 2 != 0) || (p % 2 != 0 && pmt_number_OD % 2 == 0)){
	  side_vertical_12x12_overlap_pmt_off_id[q] = pmt_number_OD;
	  /*	  PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
          PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
          PMT_OD_box_off_z.Fill(pmt_z_OD);
	  */q++;
	}
	if((p==0 && m % (2*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
	  m=0;
	  p++;
	}
        m++;
      }
    }
  }
  
  //12x12 overlap vertical cabling, 36x4
  y_separation = 102.14;
  failure_block_vertical = 36;
  failure_block_horizontal = 4;
  int* side_vertical_12x12_overlap_2_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_random_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }
  
  m = 1;
  q = 0;
  p = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
      	if((p % 2 ==0 && pmt_number_OD % 2 != 0) || (p % 2 != 0 && pmt_number_OD % 2 == 0)){
	  side_vertical_12x12_overlap_2_pmt_off_id[q] = pmt_number_OD;
	  /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
          PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
          PMT_OD_box_off_z.Fill(pmt_z_OD);
	  */q++;
	}	
	if((p==0 && m % (2*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
	  m=0;
	  p++;
	}
	m++;
      }
    }
  }

  //12x12 vertical cabling cross
  failure_block_vertical = 12;
  failure_block_horizontal = 12;
  int ver_cross_pmt_id = 3526;
  int* side_vertical_12x12_cross_pmt_off_id = (int*)malloc(failure_block_size*sizeof(int));

  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_number_OD == ver_cross_pmt_id){
      pmt_x_OD_off = pmt_x_OD;
      pmt_y_OD_off = pmt_y_OD;
      pmt_z_OD_off = pmt_z_OD;
    }
  }

  m = 1;
  q = 0;
  p = 0;
  for(int ipmt = 0; ipmt < all_pmts_tree_OD->GetEntries(); ipmt++){
    all_pmts_tree_OD->GetEntry(ipmt);
    if(pmt_ux_OD > 0 && pmt_uy_OD > 0){
      if(pmt_y_OD > pmt_y_OD_off || pmt_x_OD < (pmt_x_OD_off - failure_block_horizontal*x_separation) || pmt_y_OD < (pmt_y_OD_off - failure_block_horizontal*y_separation) || pmt_z_OD < (pmt_z_OD_off - failure_block_vertical*z_separation) || pmt_z_OD > pmt_z_OD_off) continue;
      else{
        if((p % 2 ==0 && pmt_number_OD % 2 == 0) || (p % 2 != 0 && pmt_number_OD % 2 != 0)){
          side_vertical_12x12_cross_pmt_off_id[q] = pmt_number_OD;
          /*PMT_OD_box_off_x_y_z.Fill(pmt_x_OD, pmt_y_OD, pmt_z_OD);
          PMT_OD_box_off_x_y.Fill(pmt_x_OD, pmt_y_OD);
          PMT_OD_box_off_z.Fill(pmt_z_OD);
          */q++;
        }
        if((p==0 && m % (2*failure_block_horizontal) == 0) || (p>0 && m % (3*failure_block_horizontal) == 0) ){
          m=0;
          p++;
        }
        m++;
      }
    }
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
  TH1F h_hit_time("h_hit_time","h_hit_time; hit time [ns]",
		  2000,-0.5,tmax);
  h_hit_time.SetLineWidth(2);
  h_hit_time.SetLineColor(kBlack);

  double n_hits_limit = 300.;
  TH1F h_hit_time_few_nhits("h_hit_time_few_nhits",Form("nhits < %.0f; hit time [ns]",n_hits_limit),
		  2000,-0.5,tmax);
  h_hit_time_few_nhits.SetLineWidth(2);
  h_hit_time_few_nhits.SetLineColor(kBlue);

  TH1F h_hit_time_many_nhits("h_hit_time_many_nhits",Form("nhits > %.0f; hit time [ns]",n_hits_limit),
		  2000,-0.5,tmax);
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

  double max_nhits = 1200.; // dynamic range
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
  TH1F h_digitized_nhits_top_tubes_masked("h_digitized_nhits_top_tubes_masked","top; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_top_tubes_masked.SetLineWidth(2);
  h_digitized_nhits_top_tubes_masked.SetLineColor(kRed);
  TH1F h_digitized_nhits_side_tubes_masked("h_digitized_nhits_side_tubes_masked","side; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_side_tubes_masked.SetLineWidth(2);
  h_digitized_nhits_side_tubes_masked.SetLineColor(kBlack);
  TH1F h_digitized_nhits_bottom_tubes_masked("h_digitized_nhits_bottom_tubes_masked","bottom; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_bottom_tubes_masked.SetLineWidth(2);
  h_digitized_nhits_bottom_tubes_masked.SetLineColor(kBlue);
  TH1F h_digitized_nhits_all_tubes("h_digitized_nhits_all_tubes","all; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_all_tubes.SetLineWidth(2);
  h_digitized_nhits_all_tubes.SetLineColor(kBlue);
  TH1F h_digitized_nhits_all_tubes_masked("h_digitized_all_tubes_masked","all; n of hits; n of events",100,min_nhits,max_nhits);
  h_digitized_nhits_all_tubes_masked.SetLineWidth(2);
  h_digitized_nhits_all_tubes_masked.SetLineColor(kRed);
  TH1F h_digitized_nhits_all_tubes_physics_masked("h_digitized_all_tubes_physics_masked","all; n of hits; n of events",100,0,1600);//min_nhits,max_nhits);
  h_digitized_nhits_all_tubes_physics_masked.SetLineWidth(2);
  h_digitized_nhits_all_tubes_physics_masked.SetLineColor(kRed);
  TH1F h_digitized_nhits_top_physics_tubes("h_digitized_nhits_top_physics_tubes","top; n of hits; n of events",100,1,-1);
  h_digitized_nhits_top_physics_tubes.SetLineWidth(2);
  h_digitized_nhits_top_physics_tubes.SetLineColor(kRed);
  TH1F h_digitized_nhits_side_physics_tubes("h_digitized_nhits_side_physics_tubes","side; n of hits; n of events",100,1,-1);
  h_digitized_nhits_side_physics_tubes.SetLineWidth(2);
  h_digitized_nhits_side_physics_tubes.SetLineColor(kBlack);
  TH1F h_digitized_nhits_bottom_physics_tubes("h_digitized_nhits_bottom_physics_tubes","bottom; n of hits; n of events",100,1,-1);
  h_digitized_nhits_bottom_physics_tubes.SetLineWidth(2);
  h_digitized_nhits_bottom_physics_tubes.SetLineColor(kBlue);
  TH1F h_digitized_nhits_physics_tubes("h_digitized_nhits_physics_tubes","nonzero; n of hits; n of events",100,0,1600);//1,-1);
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
  TH1F h_nhits_OD_cluster_1("h_nhits_OD_cluster_1",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1.SetLineWidth(2);
  h_nhits_OD_cluster_1.SetLineColor(kBlack);
  h_nhits_OD_cluster_1.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked("h_nhits_OD_cluster_1_masked",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_side_ver("h_nhits_OD_cluster_1_masked_side_ver",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_side_ver.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_side_ver.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_side_ver.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_side_low("h_nhits_OD_cluster_1_masked_side_low",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_side_low.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_side_low.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_side_low.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_side_high("h_nhits_OD_cluster_1_masked_side_high",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_side_high.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_side_high.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_side_high.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_top("h_nhits_OD_cluster_1_masked_top",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_top.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_top.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_top.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_top_int("h_nhits_OD_cluster_1_masked_top_int",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_top_int.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_top_int.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_top_int.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_36x4("h_nhits_OD_cluster_1_masked_ver_36x4",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_36x4.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_36x4.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_36x4.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_12x12("h_nhits_OD_cluster_1_masked_ver_12x12",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_12x12.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_12x12.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_12x12.SetFillColor(kBlack);
  
  TH1F h_nhits_OD_cluster_1_masked_ver_12x12_2("h_nhits_OD_cluster_1_masked_ver_12x12_2",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_12x12_2.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_12x12_2.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_12x12_2.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_12x12_2_cross("h_nhits_OD_cluster_1_masked_ver_12x12_2_cross",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.SetFillColor(kBlack);
 
  TH1F h_nhits_OD_cluster_1_masked_ver_16_column("h_nhits_OD_cluster_1_masked_ver_16_column",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_16_column.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_16_column.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_16_column.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_32_column("h_nhits_OD_cluster_1_masked_ver_32_column",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_32_column.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_32_column.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_32_column.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_32_column_cross("h_nhits_OD_cluster_1_masked_ver_32_column_cross",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_32_column_cross.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_32_column_cross.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_32_column_cross.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_24_column("h_nhits_OD_cluster_1_masked_ver_24_column",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_24_column.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_24_column.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_24_column.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_16_column_64readout("h_nhits_OD_cluster_1_masked_ver_16_column_64readout",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_16_column_64readout.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_16_column_64readout.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_16_column_64readout.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_32_column_64readout("h_nhits_OD_cluster_1_masked_ver_32_column_64readout",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_32_column_64readout.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_32_column_64readout.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_32_column_64readout.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_1_masked_ver_24_column_64readout("h_nhits_OD_cluster_1_masked_ver_24_column_64readout",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_1),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_1_masked_ver_24_column_64readout.SetLineWidth(2);
  h_nhits_OD_cluster_1_masked_ver_24_column_64readout.SetLineColor(kBlack);
  h_nhits_OD_cluster_1_masked_ver_24_column_64readout.SetFillColor(kBlack);

  double min_npes=0;
  double max_npes=4000;

  TH1F h_npes_OD_cluster_1("h_npes_OD_cluster_1",Form("cluster (< %.0f); n OD p.e.'s in cluster",cluster_radius_1),500,min_npes,max_npes);
  h_npes_OD_cluster_1.SetLineWidth(2);
  h_npes_OD_cluster_1.SetLineColor(kBlack);
  h_npes_OD_cluster_1.SetFillColor(kBlack);

  TH1F h_nhits_OD_cluster_2("h_nhits_OD_cluster_2",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_2),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_2.SetLineWidth(2);
  h_nhits_OD_cluster_2.SetLineColor(kBlue);
  h_nhits_OD_cluster_2.SetFillColor(kBlue);

  TH1F h_nhits_OD_cluster_2_masked("h_nhits_OD_cluster_2_masked",Form("cluster (< %.0f); n OD hits in cluster",cluster_radius_2),500,min_nhits,max_nhits);
  h_nhits_OD_cluster_2_masked.SetLineWidth(2);
  h_nhits_OD_cluster_2_masked.SetLineColor(kBlue);
  h_nhits_OD_cluster_2_masked.SetFillColor(kBlue);

  TH1F h_npes_OD_cluster_2("h_npes_OD_cluster_2",Form("cluster (< %.0f); n OD p.e.'s in cluster",cluster_radius_2),500,min_npes,max_npes);
  h_npes_OD_cluster_2.SetLineWidth(2);
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
  double nhits_OD_cluster_1;
  double npes_OD_cluster_1;
  double nhits_OD_cluster_2;
  double npes_OD_cluster_2;
  double nhits_OD_cluster_1_masked;
  double nhits_OD_cluster_1_masked_side_ver;
  double nhits_OD_cluster_1_masked_side_low;
  double nhits_OD_cluster_1_masked_side_high;
  double nhits_OD_cluster_1_masked_top;
  double nhits_OD_cluster_1_masked_top_int;
  double nhits_OD_cluster_2_masked;
  double nhits_OD_cluster_1_masked_ver_36x4;
  double nhits_OD_cluster_1_masked_ver_12x12;
  double nhits_OD_cluster_1_masked_ver_12x12_2;
  double nhits_OD_cluster_1_masked_ver_12x12_2_cross;
  double nhits_OD_cluster_1_masked_ver_16_column;
  double nhits_OD_cluster_1_masked_ver_32_column;
  double nhits_OD_cluster_1_masked_ver_32_column_cross;
  double nhits_OD_cluster_1_masked_ver_24_column;
  double nhits_OD_cluster_1_masked_ver_16_column_64readout;
  double nhits_OD_cluster_1_masked_ver_32_column_64readout;
  double nhits_OD_cluster_1_masked_ver_24_column_64readout;
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
      double impact_top_time,  impact_top_x,  impact_top_y,  impact_top_z;
      double impact_bottom_time,  impact_bottom_x,  impact_bottom_y,  impact_bottom_z;
      double impact_side_time,  impact_side_x,  impact_side_y,  impact_side_z;
      double impact_time,  impact_x,  impact_y,  impact_z;

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


	  impact_top = is_impact_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
						  track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
						  OD_radius, OD_height, muon_energy,
						     & impact_top_time, & impact_top_x, & impact_top_y, & impact_top_z);
	  impact_bottom = is_impact_horizontal_plane(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
						  track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
						     OD_radius, -OD_height, muon_energy,
						     & impact_bottom_time, & impact_bottom_x, & impact_bottom_y, & impact_bottom_z);
	  if( ! impact_top ) // muons always move down, so if they go through the top the impact point is there (even if later they cross the side)
	    impact_side = is_impact_side(track_start_x->at(itrigger).at(itrack),track_start_y->at(itrigger).at(itrack),track_start_z->at(itrigger).at(itrack),
					 track_ux->at(itrigger).at(itrack),track_uy->at(itrigger).at(itrack),track_uz->at(itrigger).at(itrack),
					 OD_radius, OD_height, muon_energy,
					 & impact_side_time, & impact_side_x, & impact_side_y, & impact_side_z);

	  if( impact_top && impact_bottom ) impact_bottom = false;
	  if( impact_top && impact_side ) impact_side = false;
	  if( impact_side && impact_bottom ){
	    if( impact_bottom_time < impact_side_time ) impact_side = false;
	    else impact_bottom = false;
	  }
	  if( impact_top || impact_side || impact_bottom ) impact = true;

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
	    h_muon_impact_vs_start_z->Fill(track_start_z->at(itrigger).at(itrack),impact_z);


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
      int nhits_bottom_physics_off = 0;
      int nhits_side_physics_off = 0;
      int nhits_all_physics_off = 0;
      int nhits_top_physics_off = 0;
      
      double total_charge = 0.;
      nhits_OD_cluster_1 = 0;
      npes_OD_cluster_1 = 0;
      nhits_OD_cluster_2 = 0;
      npes_OD_cluster_2 = 0;
      nhits_OD_cluster_1_masked = 0;
      nhits_OD_cluster_1_masked_side_ver = 0;
      nhits_OD_cluster_1_masked_side_low = 0;
      nhits_OD_cluster_1_masked_side_high = 0;
      nhits_OD_cluster_1_masked_top = 0;
      nhits_OD_cluster_1_masked_top_int = 0;
      nhits_OD_cluster_2_masked = 0;
      nhits_OD_cluster_1_masked_ver_36x4 = 0;
      nhits_OD_cluster_1_masked_ver_12x12 = 0;
      nhits_OD_cluster_1_masked_ver_12x12_2 = 0;
      nhits_OD_cluster_1_masked_ver_12x12_2_cross = 0;
      nhits_OD_cluster_1_masked_ver_16_column = 0;
      nhits_OD_cluster_1_masked_ver_32_column = 0;
      nhits_OD_cluster_1_masked_ver_32_column_cross = 0;
      nhits_OD_cluster_1_masked_ver_24_column = 0;
      nhits_OD_cluster_1_masked_ver_16_column_64readout = 0;
      nhits_OD_cluster_1_masked_ver_32_column_64readout = 0;
      nhits_OD_cluster_1_masked_ver_24_column_64readout = 0;
      int nhits_top_off = 0;
      int nhits_bottom_off = 0;
      int nhits_side_off = 0;
      int nhits_all_off = 0;
      

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

	int pmt_off_side_ver = 0;
        for(int i = 0; i< failure_block_size; i++ ){ //number_of_failed_pmts; i++ ){
	  if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_8x18_pmt_off_id[i] ) // failed_pmts_IDs[i] 
	    pmt_off_side_ver++;
	}
	int pmt_off_side_high = 0;
        for(int i = 0; i< failure_block_size; i++ ){ 
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == high_side_pmt_off_id[i] ) 
            pmt_off_side_high++;
        }
	int pmt_off_side_low = 0;
        for(int i = 0; i< failure_block_size; i++ ){ 
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == low_side_pmt_off_id[i] ) 
            pmt_off_side_low++;
        }
	int pmt_off_top =  0;
        for(int i = 0; i< topcap_failure_block_size; i++ ){ 
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == topcap_pmt_off_id[i] ) 
            pmt_off_top++;
        }
	int pmt_off_top_int = 0;
        for(int i = 0; i< 0.5*topcap_failure_block_size; i++ ){ 
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == interleaved_topcap_pmt_off_id[i] ) 
            pmt_off_top_int++;
        }	
	int pmt_off_ver_36x4 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
	  if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_36x4_pmt_off_id[i] )
	    pmt_off_ver_36x4++;
	}
	int pmt_off_ver_36x4_overlap = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_36x4_overlap_pmt_off_id[i] )
            pmt_off_ver_36x4_overlap++;
        }
	int pmt_off_ver_12x12 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_12x12_pmt_off_id[i] )
            pmt_off_ver_12x12++;
        }
	int pmt_off_ver_12x12_overlap = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_12x12_overlap_pmt_off_id[i] )
            pmt_off_ver_12x12_overlap++;
        }
	int pmt_off_ver_12x12_2_overlap = 0;
        for(int i = 0; i< failure_block_size; i++ ){
	  if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_12x12_overlap_2_pmt_off_id[i] )
            pmt_off_ver_12x12_2_overlap++;
	}
	int pmt_off_ver_12x12_2_overlap_cross = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_12x12_cross_pmt_off_id[i] )
            pmt_off_ver_12x12_2_overlap_cross++;
        }
	int pmt_off_ver_16x4 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_16x4_pmt_off_id[i] )
            pmt_off_ver_16x4++;
        }
	int pmt_off_ver_8x8 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_8x8_pmt_off_id[i] )
            pmt_off_ver_8x8++;
        }
	int pmt_off_ver_4x16 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_4x16_pmt_off_id[i] )
            pmt_off_ver_4x16++;
        }
	int pmt_off_ver_4x16_cross = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_4x16_cross_pmt_off_id[i] )
            pmt_off_ver_4x16_cross++;
        }
	int pmt_off_ver_8x8_overlap = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_8x8_offset_pmt_off_id[i] )
            pmt_off_ver_8x8_overlap++;
        }
	int pmt_off_ver_32x4 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_32x4_pmt_off_id[i] )
            pmt_off_ver_32x4++;
        }
        int pmt_off_ver_4x32 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_4x32_pmt_off_id[i] )
            pmt_off_ver_4x32++;
        }
        int pmt_off_ver_8x16 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_8x16_pmt_off_id[i] )
            pmt_off_ver_8x16++;
        }
        int pmt_off_ver_16x8 = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_16x8_pmt_off_id[i] )
            pmt_off_ver_16x8++;
        }
	int pmt_off_ver_16x8_offset = 0;
        for(int i = 0; i< failure_block_size; i++ ){
          if (digitized_hit_OD_tube_id->at(itrigger).at(idigitizedhit) == side_vertical_16x8_offset_pmt_off_id[i] )
            pmt_off_ver_16x8_offset++;
        }
	
	//if(pmt_off) continue;
	if( pmt_location_OD == 5 ){
	  nhits_top_off ++;
	  if( impact )
	    nhits_top_physics_off ++;
	}
        else if( pmt_location_OD == 4 ){
	  nhits_side_off ++;
	  if( impact )
	    nhits_side_physics_off ++;
        }
        else if( pmt_location_OD == 3 ){
          nhits_bottom_off ++;
	  if( impact )
            nhits_bottom_physics_off ++;
	}
	
	if( impact ){
          distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
          h_distance_pmt_impact.Fill(distance_pmt_impact);

          if( distance_pmt_impact <= cluster_radius_1 ){
            nhits_OD_cluster_1_masked ++;
          }
          if( distance_pmt_impact <= cluster_radius_2 ){
            nhits_OD_cluster_2_masked ++;
	  }
	}
	if(pmt_off_side_ver==0){
	  if( impact ){
	    distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	    if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_side_ver ++;
	  }
	}
	if(pmt_off_side_high==0){
	  if( impact ){
	    distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	    if( distance_pmt_impact <= cluster_radius_1 ) nhits_OD_cluster_1_masked_side_high ++;
	  }
	}
	if(pmt_off_side_low==0){
	  if( impact ){
	    distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	    if( distance_pmt_impact <= cluster_radius_1 ) nhits_OD_cluster_1_masked_side_low ++;
	  }
	}
	if(pmt_off_top==0){
	  if( impact ){
	    distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	    if( distance_pmt_impact <= cluster_radius_1 ) nhits_OD_cluster_1_masked_top ++;
	  }
	}
	if(pmt_off_top_int==0){
	  if( impact ){
	    distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	    if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_top_int ++;
	  }
	}
	if(pmt_off_ver_36x4==0 && pmt_off_ver_36x4_overlap==0){	
	  if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_36x4 ++;
	  }
        }
	if(pmt_off_ver_12x12==0 && pmt_off_ver_12x12_overlap==0){
	  if( impact ){
	    distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	    if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_12x12 ++;
	  }
	}
	if(pmt_off_ver_12x12==0 && pmt_off_ver_12x12_2_overlap==0){
	  if( impact ){
	    distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	    if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_12x12_2 ++;
	  }
	}
	if(pmt_off_ver_12x12_2_overlap_cross==0 && pmt_off_ver_12x12_2_overlap==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_12x12_2_cross ++;
          }
        }	
	if(pmt_off_ver_16x4==0 && pmt_off_ver_8x8==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_16_column ++;
          }
        }
	if(pmt_off_ver_16x4==0 && pmt_off_ver_4x16==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_32_column ++;
          }
        }
	if(pmt_off_ver_4x16_cross==0 && pmt_off_ver_16x4==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_32_column_cross ++;
          }
        }
	if(pmt_off_ver_8x8==0 && pmt_off_ver_8x8_overlap==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_24_column ++;
          }
        }
	if(pmt_off_ver_16x8==0 && pmt_off_ver_8x16==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_16_column_64readout ++;
          }
        }
	if(pmt_off_ver_32x4==0 && pmt_off_ver_4x32==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_32_column_64readout ++;
          }
        }
	if(pmt_off_ver_16x8==0 && pmt_off_ver_16x8_offset==0){
          if( impact ){
            distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
            if( distance_pmt_impact <= cluster_radius_1 )nhits_OD_cluster_1_masked_ver_24_column_64readout ++;
          }
        }
        
      }
      h_total_charge.Fill(total_charge);
      
      
      h_digitized_nhits_top_tubes.Fill(nhits_top);
      h_digitized_nhits_side_tubes.Fill(nhits_side);
      h_digitized_nhits_bottom_tubes.Fill(nhits_bottom);
      
      nhits_all = nhits_top + nhits_side + nhits_bottom;
      nhits_all_off = nhits_top_off + nhits_side_off + nhits_bottom_off;
      h_digitized_nhits_all_tubes_masked.Fill(nhits_all_off);    
      
      h_digitized_nhits_all_tubes.Fill(nhits_all);
      h_npes_vs_digitized_nhits.Fill(nhits_all, total_charge);
      h_nhits_vs_muon_energy.Fill(muon_energy,nhits_all);
      
      
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

	nhits_all_physics_off = nhits_top_physics_off + nhits_side_physics_off + nhits_bottom_physics_off;
	h_digitized_nhits_all_tubes_physics_masked.Fill(nhits_all_physics_off);

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
      
	h_nhits_OD_cluster_1_masked.Fill(nhits_OD_cluster_1_masked);
	h_nhits_OD_cluster_2_masked.Fill(nhits_OD_cluster_2_masked);
	h_nhits_OD_cluster_1_masked_top_int.Fill(nhits_OD_cluster_1_masked_top_int);
	h_nhits_OD_cluster_1_masked_top.Fill(nhits_OD_cluster_1_masked_top);
	h_nhits_OD_cluster_1_masked_side_high.Fill(nhits_OD_cluster_1_masked_side_high);
	h_nhits_OD_cluster_1_masked_side_low.Fill(nhits_OD_cluster_1_masked_side_low);
	h_nhits_OD_cluster_1_masked_side_ver.Fill(nhits_OD_cluster_1_masked_side_ver);
	h_nhits_OD_cluster_1_masked_ver_36x4.Fill(nhits_OD_cluster_1_masked_ver_36x4);
	h_nhits_OD_cluster_1_masked_ver_12x12.Fill(nhits_OD_cluster_1_masked_ver_12x12);
	h_nhits_OD_cluster_1_masked_ver_12x12_2.Fill(nhits_OD_cluster_1_masked_ver_12x12_2);
	h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.Fill(nhits_OD_cluster_1_masked_ver_12x12_2_cross);
	h_nhits_OD_cluster_1_masked_ver_16_column.Fill(nhits_OD_cluster_1_masked_ver_16_column);
	h_nhits_OD_cluster_1_masked_ver_32_column.Fill(nhits_OD_cluster_1_masked_ver_32_column);
	h_nhits_OD_cluster_1_masked_ver_32_column_cross.Fill(nhits_OD_cluster_1_masked_ver_32_column_cross);
	h_nhits_OD_cluster_1_masked_ver_24_column.Fill(nhits_OD_cluster_1_masked_ver_24_column);
	h_nhits_OD_cluster_1_masked_ver_16_column_64readout.Fill(nhits_OD_cluster_1_masked_ver_16_column_64readout);
        h_nhits_OD_cluster_1_masked_ver_32_column_64readout.Fill(nhits_OD_cluster_1_masked_ver_32_column_64readout);
        h_nhits_OD_cluster_1_masked_ver_24_column_64readout.Fill(nhits_OD_cluster_1_masked_ver_24_column_64readout);


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

	// qqq
	for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_OD_tube_id->at(itrigger)).size(); idigitizedhit++){
	  // loop on digitized hits in the trigger
	  
	  tube_id = (digitized_hit_OD_tube_id->at(itrigger)).at(idigitizedhit);
	  all_pmts_tree_OD->GetEntry(tube_id - 1);

	  distance_pmt_impact = sqrt(pow(pmt_x_OD - impact_x,2) + pow(pmt_y_OD - impact_y,2) + pow(pmt_z_OD - impact_z,2));
	  hit_time = (digitized_hit_OD_time->at(itrigger)).at(idigitizedhit);
	  
	  if( nhits_all_physics < n_hits_limit ){
	    h_distance_pmt_impact_few_nhits.Fill(distance_pmt_impact);
	    h_hit_time_few_nhits.Fill(hit_time);
	  }
	  else{
	    h_distance_pmt_impact_many_nhits.Fill(distance_pmt_impact);
	    h_hit_time_many_nhits.Fill(hit_time);
	  }
	}
      // qqq
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

  TH1F h_nhits_cluster_1_efficiency("h_nhits_cluster_1_efficiency","All PMTs; threshold [nhits]; efficiency",h_nhits_OD_cluster_1.GetNbinsX(),h_nhits_OD_cluster_1.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency.SetLineColor(h_nhits_OD_cluster_1.GetLineColor());
  h_nhits_cluster_1_efficiency.SetLineWidth(2);
  h_nhits_cluster_1_efficiency.SetStats(0);
  double nhits_integral = h_nhits_OD_cluster_1.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1.GetNbinsX(); i++){
    double local_integral = h_nhits_OD_cluster_1.Integral(i,h_nhits_OD_cluster_1.GetNbinsX());
    double nhits_ratio = local_integral/nhits_integral;
    h_nhits_cluster_1_efficiency.SetBinContent(i,nhits_ratio);
    if(nhits_ratio!=0)std::cout<<"All cluster 1 nhits "<<nhits_ratio<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_2_efficiency("h_nhits_cluster_2_efficiency","All PMTs Cluster 2; threshold [nhits]; efficiency",h_nhits_OD_cluster_2.GetNbinsX(),h_nhits_OD_cluster_2.GetXaxis()->GetXmin(),h_nhits_OD_cluster_2.GetXaxis()->GetXmax());
  h_nhits_cluster_2_efficiency.SetLineColor(h_nhits_OD_cluster_2.GetLineColor());
  h_nhits_cluster_2_efficiency.SetLineWidth(2);
  nhits_integral = h_nhits_OD_cluster_2.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_2.GetNbinsX(); i++){
    double local_integral = h_nhits_OD_cluster_2.Integral(i,h_nhits_OD_cluster_2.GetNbinsX());
    double nhits_ratio = local_integral/nhits_integral;
    h_nhits_cluster_2_efficiency.SetBinContent(i,nhits_ratio);
    // if(nhits_ratio!=0)std::cout<<"All cluster 2 nhits "<<nhits_ratio<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_2_efficiency.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked("h_nhits_cluster_1_efficiency_masked","Dead PMTs-Side Vertical; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked.GetNbinsX(),h_nhits_OD_cluster_1_masked.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked.SetLineColor(kRed);
  h_nhits_cluster_1_efficiency_masked.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked.SetStats(0);
  double nhits_integral_masked = h_nhits_OD_cluster_1_masked.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked.GetNbinsX(); i++){
    double local_integral_masked = h_nhits_OD_cluster_1_masked.Integral(i,h_nhits_OD_cluster_1_masked.GetNbinsX());
    double nhits_ratio_masked = local_integral_masked/nhits_integral_masked;
    h_nhits_cluster_1_efficiency_masked.SetBinContent(i,nhits_ratio_masked);
    if(nhits_ratio_masked!=0)std::cout<<"Masked cluster 1 nhits "<<nhits_ratio_masked<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_side_ver("h_nhits_cluster_1_efficiency_masked_side_ver","Dead PMTs - Vertical; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_side_ver.GetNbinsX(),h_nhits_OD_cluster_1_masked_side_ver.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_side_ver.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_side_ver.SetLineColor(kRed);
  h_nhits_cluster_1_efficiency_masked_side_ver.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_side_ver.SetStats(0);
  double nhits_integral_masked_side_ver = h_nhits_OD_cluster_1_masked_side_ver.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_side_ver.GetNbinsX(); i++){
    double local_integral_masked_side_ver = h_nhits_OD_cluster_1_masked_side_ver.Integral(i,h_nhits_OD_cluster_1_masked_side_ver.GetNbinsX());
    double nhits_ratio_masked_side_ver = local_integral_masked_side_ver/nhits_integral_masked_side_ver;
    h_nhits_cluster_1_efficiency_masked_side_ver.SetBinContent(i,nhits_ratio_masked_side_ver);
    if(nhits_ratio_masked_side_ver!=0)std::cout<<"Side Ver Masked cluster 1 nhits "<<nhits_ratio_masked_side_ver<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_side_ver.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_side_high("h_nhits_cluster_1_efficiency_masked_side_high","Dead PMTs - Horizontal High; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_side_high.GetNbinsX(),h_nhits_OD_cluster_1_masked_side_high.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_side_high.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_side_high.SetLineColor(kCyan);
  h_nhits_cluster_1_efficiency_masked_side_high.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_side_high.SetStats(0);
  double nhits_integral_masked_side_high = h_nhits_OD_cluster_1_masked_side_high.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_side_high.GetNbinsX(); i++){
    double local_integral_masked_side_high = h_nhits_OD_cluster_1_masked_side_high.Integral(i,h_nhits_OD_cluster_1_masked_side_high.GetNbinsX());
    double nhits_ratio_masked_side_high = local_integral_masked_side_high/nhits_integral_masked_side_high;
    h_nhits_cluster_1_efficiency_masked_side_high.SetBinContent(i,nhits_ratio_masked_side_high);
    if(nhits_ratio_masked_side_high!=0)std::cout<<"Side High Masked cluster 1 nhits "<<nhits_ratio_masked_side_high<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_side_high.GetBinCenter(i)<<std::endl;
  }
  TH1F h_nhits_cluster_1_efficiency_masked_side_low("h_nhits_cluster_1_efficiency_masked_side_low","Dead PMTs - Horizontal Low; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_side_low.GetNbinsX(),h_nhits_OD_cluster_1_masked_side_low.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_side_low.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_side_low.SetLineColor(kGray);
  h_nhits_cluster_1_efficiency_masked_side_low.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_side_low.SetStats(0);
  double nhits_integral_masked_side_low = h_nhits_OD_cluster_1_masked_side_low.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_side_low.GetNbinsX(); i++){
    double local_integral_masked_side_low = h_nhits_OD_cluster_1_masked_side_low.Integral(i,h_nhits_OD_cluster_1_masked_side_low.GetNbinsX());
    double nhits_ratio_masked_side_low = local_integral_masked_side_low/nhits_integral_masked_side_low;
    h_nhits_cluster_1_efficiency_masked_side_low.SetBinContent(i,nhits_ratio_masked_side_low);
    if(nhits_ratio_masked_side_low!=0)std::cout<<"Side low Masked cluster 1 nhits "<<nhits_ratio_masked_side_low<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_side_low.GetBinCenter(i)<<std::endl;
  }
  TH1F h_nhits_cluster_1_efficiency_masked_top("h_nhits_cluster_1_efficiency_masked_top","Dead PMTs - Topcap (8x8 - Not Interleaved); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_top.GetNbinsX(),h_nhits_OD_cluster_1_masked_top.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_top.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_top.SetLineColor(kGreen);
  h_nhits_cluster_1_efficiency_masked_top.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_top.SetStats(0);
  double nhits_integral_masked_top = h_nhits_OD_cluster_1_masked_top.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_top.GetNbinsX(); i++){
    double local_integral_masked_top = h_nhits_OD_cluster_1_masked_top.Integral(i,h_nhits_OD_cluster_1_masked_top.GetNbinsX());
    double nhits_ratio_masked_top = local_integral_masked_top/nhits_integral_masked_top;
    h_nhits_cluster_1_efficiency_masked_top.SetBinContent(i,nhits_ratio_masked_top);
    if(nhits_ratio_masked_top!=0)std::cout<<"Top Masked cluster 1 nhits "<<nhits_ratio_masked_top<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_top.GetBinCenter(i)<<std::endl;
  }
  TH1F h_nhits_cluster_1_efficiency_masked_top_int("h_nhits_cluster_1_efficiency_masked_top_int","Dead PMTs- Topcap (8x8 - Interleaved); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_top_int.GetNbinsX(),h_nhits_OD_cluster_1_masked_top_int.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_top_int.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_top_int.SetLineColor(kBlue);
  h_nhits_cluster_1_efficiency_masked_top_int.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_top_int.SetStats(0);
  double nhits_integral_masked_top_int = h_nhits_OD_cluster_1_masked_top_int.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_top_int.GetNbinsX(); i++){
    double local_integral_masked_top_int = h_nhits_OD_cluster_1_masked_top_int.Integral(i,h_nhits_OD_cluster_1_masked_top_int.GetNbinsX());
    double nhits_ratio_masked_top_int = local_integral_masked_top_int/nhits_integral_masked_top_int;
    h_nhits_cluster_1_efficiency_masked_top_int.SetBinContent(i,nhits_ratio_masked_top_int);
    if(nhits_ratio_masked_top_int!=0)std::cout<<"Top int Masked cluster 1 nhits "<<nhits_ratio_masked_top_int<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_top_int.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_36x4("h_nhits_cluster_1_efficiency_masked_ver_36x4","16 Column, 72 Channel (36x4) + (18x8); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_36x4.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_36x4.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_36x4.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_36x4.SetLineColor(kCyan);
  h_nhits_cluster_1_efficiency_masked_ver_36x4.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_36x4.SetStats(0);
  double nhits_integral_masked_ver_36x4 = h_nhits_OD_cluster_1_masked_ver_36x4.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_36x4.GetNbinsX(); i++){
    double local_integral_masked_ver_36x4 = h_nhits_OD_cluster_1_masked_ver_36x4.Integral(i,h_nhits_OD_cluster_1_masked_ver_36x4.GetNbinsX());
    double nhits_ratio_masked_ver_36x4 = local_integral_masked_ver_36x4/nhits_integral_masked_ver_36x4;
    h_nhits_cluster_1_efficiency_masked_ver_36x4.SetBinContent(i,nhits_ratio_masked_ver_36x4);
    if(nhits_ratio_masked_ver_36x4!=0)std::cout<<"Vertical Cabling 36x4 Masked cluster 1 nhits "<<nhits_ratio_masked_ver_36x4<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_36x4.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_12x12("h_nhits_cluster_1_efficiency_masked_ver_12x12"," (12x12) + (24x6); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_12x12.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_12x12.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_12x12.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_12x12.SetLineColor(kMagenta);
  h_nhits_cluster_1_efficiency_masked_ver_12x12.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_12x12.SetStats(0);
  double nhits_integral_masked_ver_12x12 = h_nhits_OD_cluster_1_masked_ver_12x12.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_12x12.GetNbinsX(); i++){
    double local_integral_masked_ver_12x12 = h_nhits_OD_cluster_1_masked_ver_12x12.Integral(i,h_nhits_OD_cluster_1_masked_ver_12x12.GetNbinsX());
    double nhits_ratio_masked_ver_12x12 = local_integral_masked_ver_12x12/nhits_integral_masked_ver_12x12;
    h_nhits_cluster_1_efficiency_masked_ver_12x12.SetBinContent(i,nhits_ratio_masked_ver_12x12);
    if(nhits_ratio_masked_ver_12x12!=0)std::cout<<"Vertical Cabling 12x12 Masked cluster 1 nhits "<<nhits_ratio_masked_ver_12x12<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_12x12.GetBinCenter(i)<<std::endl;
  }
  
  TH1F h_nhits_cluster_1_efficiency_masked_ver_12x12_2("h_nhits_cluster_1_efficiency_masked_ver_12x12_2","24 Column, 72 Channel (12x12) + (36x4) L-Shape; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_12x12_2.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_12x12_2.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_12x12_2.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2.SetLineColor(kBlue);
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2.SetStats(0);
  double nhits_integral_masked_ver_12x12_2 = h_nhits_OD_cluster_1_masked_ver_12x12_2.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_12x12_2.GetNbinsX(); i++){
    double local_integral_masked_ver_12x12_2 = h_nhits_OD_cluster_1_masked_ver_12x12_2.Integral(i,h_nhits_OD_cluster_1_masked_ver_12x12_2.GetNbinsX());
    double nhits_ratio_masked_ver_12x12_2 = local_integral_masked_ver_12x12_2/nhits_integral_masked_ver_12x12_2;
    h_nhits_cluster_1_efficiency_masked_ver_12x12_2.SetBinContent(i,nhits_ratio_masked_ver_12x12_2);
    if(nhits_ratio_masked_ver_12x12_2!=0)std::cout<<"Vertical Cabling 12x12 Masked cluster 1 nhits "<<nhits_ratio_masked_ver_12x12_2<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_12x12_2.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross("h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross","24 Column, 72 Channel (12x12) + (36x4) - Cross; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross.SetLineColor(kGreen+3);
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross.SetStats(0);
  double nhits_integral_masked_ver_12x12_2_cross = h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.GetNbinsX(); i++){
    double local_integral_masked_ver_12x12_2_cross = h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.Integral(i,h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.GetNbinsX());
    double nhits_ratio_masked_ver_12x12_2_cross = local_integral_masked_ver_12x12_2_cross/nhits_integral_masked_ver_12x12_2_cross;
    h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross.SetBinContent(i,nhits_ratio_masked_ver_12x12_2_cross);
    if(nhits_ratio_masked_ver_12x12_2_cross!=0)std::cout<<"Vertical Cabling 12x12 Masked cluster 1 nhits "<<nhits_ratio_masked_ver_12x12_2_cross<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_16_column("h_nhits_cluster_1_efficiency_masked_ver_16_column","16 Column, 32 Channel (16x4) + (8x8); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_16_column.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_16_column.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_16_column.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_16_column.SetLineColor(kCyan);
  h_nhits_cluster_1_efficiency_masked_ver_16_column.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_16_column.SetStats(0);
  double nhits_integral_masked_ver_16_column = h_nhits_OD_cluster_1_masked_ver_16_column.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_16_column.GetNbinsX(); i++){
    double local_integral_masked_ver_16_column = h_nhits_OD_cluster_1_masked_ver_16_column.Integral(i,h_nhits_OD_cluster_1_masked_ver_16_column.GetNbinsX());
    double nhits_ratio_masked_ver_16_column = local_integral_masked_ver_16_column/nhits_integral_masked_ver_16_column;
    h_nhits_cluster_1_efficiency_masked_ver_16_column.SetBinContent(i,nhits_ratio_masked_ver_16_column);
    if(nhits_ratio_masked_ver_16_column!=0)std::cout<<"Vertical Cabling 16 column Masked cluster 1 nhits "<<nhits_ratio_masked_ver_16_column<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_16_column.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_32_column("h_nhits_cluster_1_efficiency_masked_ver_32_column","32 Column, 32 Channel (16x4) + (4x16) L-Shape; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_32_column.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_32_column.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_32_column.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_32_column.SetLineColor(kRed);
  h_nhits_cluster_1_efficiency_masked_ver_32_column.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_32_column.SetStats(0);
  double nhits_integral_masked_ver_32_column = h_nhits_OD_cluster_1_masked_ver_32_column.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_32_column.GetNbinsX(); i++){
    double local_integral_masked_ver_32_column = h_nhits_OD_cluster_1_masked_ver_32_column.Integral(i,h_nhits_OD_cluster_1_masked_ver_32_column.GetNbinsX());
    double nhits_ratio_masked_ver_32_column = local_integral_masked_ver_32_column/nhits_integral_masked_ver_32_column;
    h_nhits_cluster_1_efficiency_masked_ver_32_column.SetBinContent(i,nhits_ratio_masked_ver_32_column);
    if(nhits_ratio_masked_ver_32_column!=0)std::cout<<"Vertical Cabling 32 column Masked cluster 1 nhits "<<nhits_ratio_masked_ver_32_column<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_32_column.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_32_column_cross("h_nhits_cluster_1_efficiency_masked_ver_32_column_cross","32 Column, 32 Channel (16x4) + (4x16) Cross; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_32_column_cross.SetLineColor(kGreen);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_cross.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_cross.SetStats(0);
  double nhits_integral_masked_ver_32_column_cross = h_nhits_OD_cluster_1_masked_ver_32_column_cross.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX(); i++){
    double local_integral_masked_ver_32_column_cross = h_nhits_OD_cluster_1_masked_ver_32_column_cross.Integral(i,h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX());
    double nhits_ratio_masked_ver_32_column_cross = local_integral_masked_ver_32_column_cross/nhits_integral_masked_ver_32_column_cross;
  h_nhits_cluster_1_efficiency_masked_ver_32_column_cross.SetBinContent(i,nhits_ratio_masked_ver_32_column_cross);
  if(nhits_ratio_masked_ver_32_column_cross!=0)std::cout<<"Vertical Cabling 32 column cross Masked cluster 1 nhits "<<nhits_ratio_masked_ver_32_column_cross<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_32_column_cross.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_32_column_average("h_nhits_cluster_1_efficiency_masked_ver_32_column_average","32 Column, 32 Channel (16x4) + (4x16) Average; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average.SetLineColor(kMagenta);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average.SetStats(0);
  double nhits_integral_1 = h_nhits_OD_cluster_1_masked_ver_32_column.Integral();
  double nhits_integral_2 = h_nhits_OD_cluster_1_masked_ver_32_column_cross.Integral();
  double nhits_integral_average = (((2*nhits_integral_1)+nhits_integral_2)/3);
  std::cout<<"nhits_integral_average"<<nhits_integral_average<<std::endl;
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX(); i++){
    double local_integral_masked_ver_32_column_crossa = h_nhits_OD_cluster_1_masked_ver_32_column_cross.Integral(i,h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX());
    double local_integral_masked_ver_32_columna = h_nhits_OD_cluster_1_masked_ver_32_column.Integral(i,h_nhits_OD_cluster_1_masked_ver_32_column.GetNbinsX()); 
    double local_integral_average = (((2*local_integral_masked_ver_32_columna)+local_integral_masked_ver_32_column_crossa)/3);
    std::cout<<"local_integral_average"<<local_integral_average<<" column_crossa "<<local_integral_masked_ver_32_column_crossa<<std::endl; 
    double nhits_ratio_masked_ver_average = local_integral_average/nhits_integral_average;
    h_nhits_cluster_1_efficiency_masked_ver_32_column_average.SetBinContent(i,nhits_ratio_masked_ver_average);
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_32_column_average_72readout("h_nhits_cluster_1_efficiency_masked_ver_32_column_average_72readout","24 Column, 72 Channel (13x23 + 36x4) Average; threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetXaxis()-> GetXmin(),h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average_72readout.SetLineColor(kMagenta);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average_72readout.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average_72readout.SetStats(0);
  double nhits_integral_12 = h_nhits_OD_cluster_1_masked_ver_12x12_2.Integral();
  double nhits_integral_22 = h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.Integral();
  double nhits_integral_average2 = (((2*nhits_integral_12)+nhits_integral_22)/3);
  std::cout<<"nhits_integral_average"<<nhits_integral_average2<<std::endl;
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_32_column_cross.GetNbinsX(); i++){
    double local_integral_masked_ver_32_column_crossa2 = h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.Integral(i,h_nhits_OD_cluster_1_masked_ver_12x12_2_cross.GetNbinsX());
    double local_integral_masked_ver_32_columna2 = h_nhits_OD_cluster_1_masked_ver_12x12_2.Integral(i,h_nhits_OD_cluster_1_masked_ver_12x12_2.GetNbinsX());
    double local_integral_average2 = (((2*local_integral_masked_ver_32_columna2)+local_integral_masked_ver_32_column_crossa2)/3);
    std::cout<<"local_integral_average"<<local_integral_average2<<" column_crossa "<<local_integral_masked_ver_32_column_crossa2<<std::endl;
    double nhits_ratio_masked_ver_average2 = local_integral_average2/nhits_integral_average2;
    h_nhits_cluster_1_efficiency_masked_ver_32_column_average_72readout.SetBinContent(i,nhits_ratio_masked_ver_average2);
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_24_column("h_nhits_cluster_1_efficiency_masked_ver_24_column","24 Column, 32 Channel (8x8) + (8x8); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_24_column.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_24_column.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_24_column.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_24_column.SetLineColor(kBlue);
  h_nhits_cluster_1_efficiency_masked_ver_24_column.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_24_column.SetStats(0);
  double nhits_integral_masked_ver_24_column = h_nhits_OD_cluster_1_masked_ver_24_column.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_24_column.GetNbinsX(); i++){
    double local_integral_masked_ver_24_column = h_nhits_OD_cluster_1_masked_ver_24_column.Integral(i,h_nhits_OD_cluster_1_masked_ver_24_column.GetNbinsX());
    double nhits_ratio_masked_ver_24_column = local_integral_masked_ver_24_column/nhits_integral_masked_ver_24_column;
    h_nhits_cluster_1_efficiency_masked_ver_24_column.SetBinContent(i,nhits_ratio_masked_ver_24_column);
    if(nhits_ratio_masked_ver_24_column!=0)std::cout<<"Vertical Cabling 24 column Masked cluster 1 nhits "<<nhits_ratio_masked_ver_24_column<<" bin "<<i<<"nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_24_column.GetBinCenter(i)<<std::endl;
  }

  TH1F h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout("h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout","32 Column, 64 Channel (4x32) + (32x4); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_32_column_64readout.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_32_column_64readout.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_32_column_64readout.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
  h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout.SetLineColor(kRed);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout.SetLineWidth(2);
  h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout.SetStats(0);
  double nhits_integral_masked_ver_32_column_64readout = h_nhits_OD_cluster_1_masked_ver_32_column_64readout.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_32_column_64readout.GetNbinsX(); i++){
  double local_integral_masked_ver_32_column_64readout = h_nhits_OD_cluster_1_masked_ver_32_column_64readout.Integral(i,h_nhits_OD_cluster_1_masked_ver_32_column_64readout.GetNbinsX());
  double nhits_ratio_masked_ver_32_column_64readout = local_integral_masked_ver_32_column_64readout/nhits_integral_masked_ver_32_column_64readout;
  h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout.SetBinContent(i,nhits_ratio_masked_ver_32_column_64readout);
  if(nhits_ratio_masked_ver_32_column_64readout!=0)std::cout<<"Vertical Cabling 32 column  64 readout Masked cluster 1 nhits "<<nhits_ratio_masked_ver_32_column_64readout<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout.GetBinCenter(i)<<std::endl;
 }

 TH1F h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout("h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout","24 Column, 64 Channel (16x8) + (16x8, offset); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_24_column_64readout.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_24_column_64readout.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_24_column_64readout.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
 h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout.SetLineColor(kBlue);
h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout.SetLineWidth(2);
h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout.SetStats(0);
double nhits_integral_masked_ver_24_column_64readout = h_nhits_OD_cluster_1_masked_ver_24_column_64readout.Integral();
for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_24_column_64readout.GetNbinsX(); i++){
  double local_integral_masked_ver_24_column_64readout = h_nhits_OD_cluster_1_masked_ver_24_column_64readout.Integral(i,h_nhits_OD_cluster_1_masked_ver_24_column_64readout.GetNbinsX());
  double nhits_ratio_masked_ver_24_column_64readout = local_integral_masked_ver_24_column_64readout/nhits_integral_masked_ver_24_column_64readout;
  h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout.SetBinContent(i,nhits_ratio_masked_ver_24_column_64readout);
  if(nhits_ratio_masked_ver_24_column_64readout!=0)std::cout<<"Vertical Cabling 24 column 64 readout Masked cluster 1 nhits "<<nhits_ratio_masked_ver_24_column_64readout<<" bin "<<i<<"nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout.GetBinCenter(i)<<std::endl;
 }

 TH1F h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout("h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout","16 Column, 64 channel (16x8) + (8x16); threshold [nhits]; efficiency",h_nhits_OD_cluster_1_masked_ver_16_column_64readout.GetNbinsX(),h_nhits_OD_cluster_1_masked_ver_16_column_64readout.GetXaxis()->GetXmin(),h_nhits_OD_cluster_1_masked_ver_16_column_64readout.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency,cluster_radius_1)
 h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout.SetLineColor(kCyan);
h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout.SetLineWidth(2);
h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout.SetStats(0);
double nhits_integral_masked_ver_16_column_64readout = h_nhits_OD_cluster_1_masked_ver_16_column_64readout.Integral();
for(int i=1; i<=h_nhits_OD_cluster_1_masked_ver_16_column_64readout.GetNbinsX(); i++){
  double local_integral_masked_ver_16_column_64readout = h_nhits_OD_cluster_1_masked_ver_16_column_64readout.Integral(i,h_nhits_OD_cluster_1_masked_ver_16_column_64readout.GetNbinsX());
  double nhits_ratio_masked_ver_16_column_64readout = local_integral_masked_ver_16_column_64readout/nhits_integral_masked_ver_16_column_64readout;
  h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout.SetBinContent(i,nhits_ratio_masked_ver_16_column_64readout);
  if(nhits_ratio_masked_ver_16_column_64readout!=0)std::cout<<"Vertical Cabling 16 column 64 readout Masked cluster 1 nhits "<<nhits_ratio_masked_ver_16_column_64readout<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout.GetBinCenter(i)<<std::endl;
 }

  TH1F h_nhits_cluster_2_efficiency_masked("h_nhits_cluster_2_efficiency_masked","Dead PMTs Cluster 2; threshold [nhits]; efficiency",h_nhits_OD_cluster_2_masked.GetNbinsX(),h_nhits_OD_cluster_2_masked.GetXaxis()->GetXmin(),h_nhits_OD_cluster_2.GetXaxis()->GetXmax()); //Form("cluster (< %.0f); threshold [nhits]; efficiency",cluster_radius_2)
  h_nhits_cluster_2_efficiency_masked.SetLineColor(kGreen);
  h_nhits_cluster_2_efficiency_masked.SetLineWidth(2);
  nhits_integral_masked = h_nhits_OD_cluster_2_masked.Integral();
  for(int i=1; i<=h_nhits_OD_cluster_2_masked.GetNbinsX(); i++){
    double local_integral_masked = h_nhits_OD_cluster_2_masked.Integral(i,h_nhits_OD_cluster_2_masked.GetNbinsX());
    double nhits_ratio_masked = local_integral_masked/nhits_integral_masked;
    h_nhits_cluster_2_efficiency_masked.SetBinContent(i,nhits_ratio_masked);
    //if(nhits_ratio_masked!=0)std::cout<<"Masked cluster 2 nhits "<<nhits_ratio_masked<<" bin "<<i<<" nhits threshold "<<h_nhits_cluster_2_efficiency_masked.GetBinCenter(i)<<std::endl;
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

  of->cd();
  PMT_x_y_z.Write();
  PMT_OD_x_y_z.Write();
  PMT_OD_off_x_y_z.Write();
  PMT_OD_x_y.Write();
  PMT_OD_off_x_y.Write();
  PMT_OD_z.Write();
  PMT_OD_off_z.Write();
  PMT_OD_box_off_z.Write();
  PMT_OD_box_off_x_y.Write();
  PMT_OD_box_off_x_y_z.Write();
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
  h_total_charge.Write();
  h_n_hits_vs_radius.Write();
  h_distance_pmt_impact.Write();
  h_distance_pmt_impact_few_nhits.Write();
  h_distance_pmt_impact_many_nhits.Write();
  h_nhits_OD_cluster_1.Write();
  h_npes_OD_cluster_1.Write();
  h_nhits_OD_cluster_1_masked.Write();
  h_nhits_OD_cluster_1_masked_ver_36x4.Write();
  h_nhits_OD_cluster_1_masked_side_ver.Write();
  h_nhits_cluster_1_efficiency_masked_top_int.Write();
  h_nhits_cluster_1_efficiency_masked_side_low.Write();
  h_nhits_cluster_1_efficiency_masked_side_high.Write();
  h_nhits_cluster_1_efficiency_masked_side_ver.Write();
  h_nhits_cluster_1_efficiency_masked_top.Write();
  h_nhits_cluster_1_efficiency_masked_ver_36x4.Write();
  h_nhits_cluster_1_efficiency_masked_ver_12x12.Write();
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2.Write();
  h_nhits_cluster_1_efficiency_masked_ver_12x12_2_cross.Write();
  h_nhits_cluster_1_efficiency_masked_ver_16_column.Write();
  h_nhits_cluster_1_efficiency_masked_ver_32_column.Write();
  h_nhits_cluster_1_efficiency_masked_ver_32_column_cross.Write();
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average.Write();  
  h_nhits_cluster_1_efficiency_masked_ver_32_column_average_72readout.Write();
  h_nhits_cluster_1_efficiency_masked_ver_24_column.Write();
  h_nhits_cluster_1_efficiency_masked_ver_16_column_64readout.Write();
  h_nhits_cluster_1_efficiency_masked_ver_32_column_64readout.Write();
  h_nhits_cluster_1_efficiency_masked_ver_24_column_64readout.Write();
  h_nhits_OD_cluster_2.Write();
  h_npes_OD_cluster_2.Write();
  h_nhits_OD_cluster_2_masked.Write();
  h_digitized_n_photons_ids.Write();
  h_digitized_photons_id0.Write();
  h_digitized_nhits_top_tubes.Write();
  h_digitized_nhits_side_tubes.Write();
  h_digitized_nhits_bottom_tubes.Write();
  h_digitized_nhits_top_tubes_masked.Write();
  h_digitized_nhits_side_tubes_masked.Write();
  h_digitized_nhits_bottom_tubes_masked.Write();
  h_digitized_nhits_all_tubes_masked.Write();
  h_digitized_nhits_all_tubes_physics_masked.Write();
  h_digitized_nhits_all_tubes.Write();
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
  h_nhits_cluster_1_efficiency_masked.Write();
  h_nhits_cluster_2_efficiency_masked.Write();

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
