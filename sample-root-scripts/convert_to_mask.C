#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TFile.h>
#include <TTree.h>

int main(){

  // open input file
  TFile *f = new TFile("output.root","READ");
  
  // create output file
  TFile *of = new TFile("mask.root", "RECREATE");

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

  Int_t i_ncapturecount, i_neutroncount;
  const Int_t MAXNCAPTURES = 25;
  Double_t i_capt_x[MAXNCAPTURES], i_capt_y[MAXNCAPTURES], i_capt_z[MAXNCAPTURES], i_capt_t0[MAXNCAPTURES], i_capt_E[MAXNCAPTURES];
  Int_t i_capt_num[MAXNCAPTURES], i_capt_pid[MAXNCAPTURES], i_capt_nucleus[MAXNCAPTURES], i_capt_nphot[MAXNCAPTURES], i_capt_ngamma[MAXNCAPTURES];
  primary_events_tree->SetBranchAddress("ncapturecount",&i_ncapturecount);
  primary_events_tree->SetBranchAddress("neutroncount",&i_neutroncount);
  primary_events_tree->SetBranchAddress("capt_x",i_capt_x);
  primary_events_tree->SetBranchAddress("capt_y",i_capt_y);
  primary_events_tree->SetBranchAddress("capt_z",i_capt_z);
  primary_events_tree->SetBranchAddress("capt_t0",i_capt_t0);
  primary_events_tree->SetBranchAddress("capt_E",i_capt_E);
  primary_events_tree->SetBranchAddress("capt_num",i_capt_num);
  primary_events_tree->SetBranchAddress("capt_pid",i_capt_pid);
  primary_events_tree->SetBranchAddress("capt_nucleus",i_capt_nucleus);
  primary_events_tree->SetBranchAddress("capt_nphot",i_capt_nphot);
  primary_events_tree->SetBranchAddress("capt_ngamma",i_capt_ngamma);
  Int_t i_nneutrons;
  primary_events_tree->SetBranchAddress("nneutrons",&i_nneutrons);

  double cmTOmm = 10.;

  TTree output_hits_tree("HitsTree","HitsTree");
  Int_t evt;
  output_hits_tree.Branch("evt",&evt,"evt/I");
  const Int_t MAXNHITS = 100000;
  Int_t nhits;
  output_hits_tree.Branch("nhits",&nhits,"nhits/I");
  Double_t hit_time[MAXNHITS], hit_x[MAXNHITS], hit_y[MAXNHITS], hit_z[MAXNHITS];
  Int_t hit_PMTid[MAXNHITS];
  output_hits_tree.Branch("hit_time",hit_time,"hit_time[nhits]/D");
  output_hits_tree.Branch("hit_x",hit_x,"hit_x[nhits]/D");
  output_hits_tree.Branch("hit_y",hit_y,"hit_y[nhits]/D");
  output_hits_tree.Branch("hit_z",hit_z,"hit_z[nhits]/D");
  output_hits_tree.Branch("hit_PMTid",hit_PMTid,"hit_PMTid[nhits]/I");
  const Int_t MAXNPARTS = 15000;
  Int_t npart;
  Double_t part_xStart[MAXNPARTS], part_yStart[MAXNPARTS], part_zStart[MAXNPARTS], part_tStart[MAXNPARTS], part_xEnd[MAXNPARTS], part_yEnd[MAXNPARTS], part_zEnd[MAXNPARTS], part_tEnd[MAXNPARTS], part_pxStart[MAXNPARTS], part_pyStart[MAXNPARTS], part_pzStart[MAXNPARTS], part_pxEnd[MAXNPARTS], part_pyEnd[MAXNPARTS], part_pzEnd[MAXNPARTS], part_KEstart[MAXNPARTS], part_KEend[MAXNPARTS];
  Int_t part_processStart[MAXNPARTS], part_processEnd[MAXNPARTS], part_parentid[MAXNPARTS], part_trackid[MAXNPARTS], part_pid[MAXNPARTS];
  output_hits_tree.Branch("npart",&npart,"npart/I");
  output_hits_tree.Branch("part_xStart",part_xStart,"part_xStart[npart]/D");
  output_hits_tree.Branch("part_yStart",part_yStart,"part_yStart[npart]/D");
  output_hits_tree.Branch("part_zStart",part_zStart,"part_zStart[npart]/D");
  output_hits_tree.Branch("part_tStart",part_tStart,"part_tStart[npart]/D");
  output_hits_tree.Branch("part_xEnd",part_xEnd,"part_xEnd[npart]/D");
  output_hits_tree.Branch("part_yEnd",part_yEnd,"part_yEnd[npart]/D");
  output_hits_tree.Branch("part_zEnd",part_zEnd,"part_zEnd[npart]/D");
  output_hits_tree.Branch("part_tEnd",part_tEnd,"part_tEnd[npart]/D");
  output_hits_tree.Branch("part_pxStart",part_pxStart,"part_pxStart[npart]/D");
  output_hits_tree.Branch("part_pyStart",part_pyStart,"part_pyStart[npart]/D");
  output_hits_tree.Branch("part_pzStart",part_pzStart,"part_pzStart[npart]/D");
  output_hits_tree.Branch("part_pxEnd",part_pxEnd,"part_pxEnd[npart]/D");
  output_hits_tree.Branch("part_pyEnd",part_pyEnd,"part_pyEnd[npart]/D");
  output_hits_tree.Branch("part_pzEnd",part_pzEnd,"part_pzEnd[npart]/D");
  output_hits_tree.Branch("part_KEstart",part_KEstart,"part_KEstart[npart]/D");
  output_hits_tree.Branch("part_KEend",part_KEend,"part_KEend[npart]/D");
  output_hits_tree.Branch("part_processStart",part_processStart,"part_processStart[npart]/I");
  output_hits_tree.Branch("part_processEnd",part_processEnd,"part_processEnd[npart]/I");
  output_hits_tree.Branch("part_parentid",part_parentid,"part_parentid[npart]/I");
  output_hits_tree.Branch("part_trackid",part_trackid,"part_trackid[npart]/I");
  output_hits_tree.Branch("part_pid",part_pid,"part_pid[npart]/I");
  Int_t ncapturecount, neutroncount;
  Double_t capt_x[MAXNCAPTURES], capt_y[MAXNCAPTURES], capt_z[MAXNCAPTURES], capt_t0[MAXNCAPTURES], capt_E[MAXNCAPTURES];
  Int_t capt_num[MAXNCAPTURES], capt_pid[MAXNCAPTURES], capt_nucleus[MAXNCAPTURES], capt_nphot[MAXNCAPTURES], capt_ngamma[MAXNCAPTURES];
  output_hits_tree.Branch("ncapturecount",&ncapturecount,"ncapturecount/I");
  output_hits_tree.Branch("neutroncount",&neutroncount,"neutroncount/I");
  output_hits_tree.Branch("capt_x",capt_x,"capt_x[ncapturecount]/D");
  output_hits_tree.Branch("capt_y",capt_y,"capt_y[ncapturecount]/D");
  output_hits_tree.Branch("capt_z",capt_z,"capt_z[ncapturecount]/D");
  output_hits_tree.Branch("capt_t0",capt_t0,"capt_t0[ncapturecount]/D");
  output_hits_tree.Branch("capt_E",capt_E,"capt_E[ncapturecount]/D");
  output_hits_tree.Branch("capt_num",capt_num,"capt_num[ncapturecount]/I");
  output_hits_tree.Branch("capt_pid",capt_pid,"capt_pid[ncapturecount]/I");
  output_hits_tree.Branch("capt_nucleus",capt_nucleus,"capt_nucleus[ncapturecount]/I");
  output_hits_tree.Branch("capt_nphot",capt_nphot,"capt_nphot[ncapturecount]/I");
  output_hits_tree.Branch("capt_ngamma",capt_ngamma,"capt_ngamma[ncapturecount]/I");
  Int_t mode, neutrino_id, ntrks, nneutrons;
  Double_t neutrino_E, neutrino_px, neutrino_py, neutrino_pz, vtxx, vtxy, vtxz;
  output_hits_tree.Branch("mode",&mode,"mode/I");
  output_hits_tree.Branch("neutrino_E",&neutrino_E,"neutrino_E/D");
  output_hits_tree.Branch("neutrino_id",&neutrino_id,"neutrino_id/I");
  output_hits_tree.Branch("neutrino_px",&neutrino_px,"neutrino_px/D");
  output_hits_tree.Branch("neutrino_py",&neutrino_py,"neutrino_py/D");
  output_hits_tree.Branch("neutrino_pz",&neutrino_pz,"neutrino_pz/D");
  output_hits_tree.Branch("ntrks",&ntrks,"ntrks/I");
  output_hits_tree.Branch("nneutrons",&nneutrons,"nneutrons/I");
  output_hits_tree.Branch("vtxx",&vtxx,"vtxx/D");
  output_hits_tree.Branch("vtxy",&vtxy,"vtxy/D");
  output_hits_tree.Branch("vtxz",&vtxz,"vtxz/D");
  const Int_t MAXNTRACKS = 40000;
  Int_t mpid[MAXNTRACKS];
  Double_t px[MAXNTRACKS], py[MAXNTRACKS], pz[MAXNTRACKS], KE[MAXNTRACKS];
  output_hits_tree.Branch("mpid",mpid,"mpid[ntrks]/I");
  output_hits_tree.Branch("px",px,"px[ntrks]/D");
  output_hits_tree.Branch("py",py,"py[ntrks]/D");
  output_hits_tree.Branch("pz",pz,"pz[ntrks]/D");
  output_hits_tree.Branch("KE",KE,"KE[ntrks]/D");

  int modulo = 100;

  int local_tube_id;

  bool offset_times_to_zero = true;
  double min_time;

  int ntrks_counter;

  if( offset_times_to_zero ){
    double local_time;
    min_time = 99999999999999999999999.;
    for(int ievent=0; ievent<primary_events_tree->GetEntries(); ievent++){
      // loop on primary events
      primary_events_tree->GetEntry(ievent); 

      for(size_t itrigger=0; itrigger<trigger_ntrack->size(); itrigger++){
	// loop on triggers in the event

	for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_tube_id->at(itrigger)).size(); idigitizedhit++){
	  // loop on digitized hits in the trigger
	  
	  local_time= (digitized_hit_time->at(itrigger)).at(idigitizedhit);
	  if( local_time < min_time ) min_time = local_time;
	}
      }
    }
    std::clog << " subtracting " << min_time << " from measured times " << std::endl;
  }

  for(int ievent=0; ievent<primary_events_tree->GetEntries(); ievent++){

    if( ievent % modulo == 0 ){
      std::clog << " converting event " << ievent << " of " << primary_events_tree->GetEntries() << std::endl;
    }

    evt = ievent;

    // loop on primary events
    primary_events_tree->GetEntry(ievent); 

    ntrks_counter = 0;
    for(size_t itrigger=0; itrigger<trigger_ntrack->size(); itrigger++){
      // loop on triggers in the event

      nhits = (digitized_hit_tube_id->at(itrigger)).size();
      vtxx = trigger_vtx_x->at(itrigger).at(0);
      vtxy = trigger_vtx_y->at(itrigger).at(0);
      vtxz = trigger_vtx_z->at(itrigger).at(0);

      for(int itrack=0; itrack<trigger_ntrack->at(itrigger); itrack++){
	// loop on tracks in the trigger

	if( itrack >= MAXNTRACKS || itrack >= MAXNPARTS) {
	  std::clog << " problem: itrack " << itrack << " max " << MAXNTRACKS << std::endl;
	  continue;
	}

	if( (track_ipnu->at(itrigger)).at(itrack) == 14 ){
	  neutrino_id = (track_ipnu->at(itrigger)).at(itrack);
	  neutrino_px = (track_ux->at(itrigger)).at(itrack);
	  neutrino_py = (track_uy->at(itrigger)).at(itrack);
	  neutrino_pz = (track_uz->at(itrigger)).at(itrack);
	  neutrino_E = (track_E->at(itrigger)).at(itrack);
	}else{

	  part_xStart[ntrks_counter] = (track_start_x->at(itrigger)).at(itrack)*cmTOmm;
	  part_yStart[ntrks_counter] = (track_start_y->at(itrigger)).at(itrack)*cmTOmm;
	  part_zStart[ntrks_counter] = (track_start_z->at(itrigger)).at(itrack)*cmTOmm;
	  part_tStart[ntrks_counter] = (track_time->at(itrigger)).at(itrack);
	  part_xEnd[ntrks_counter] = (track_stop_x->at(itrigger)).at(itrack)*cmTOmm;
	  part_yEnd[ntrks_counter] = (track_stop_y->at(itrigger)).at(itrack)*cmTOmm;
	  part_zEnd[ntrks_counter] = (track_stop_z->at(itrigger)).at(itrack)*cmTOmm;
	  part_pxStart[ntrks_counter] = (track_ux->at(itrigger)).at(itrack);
	  part_pyStart[ntrks_counter] = (track_uy->at(itrigger)).at(itrack);
	  part_pzStart[ntrks_counter] = (track_uz->at(itrigger)).at(itrack);
	  part_pxEnd[ntrks_counter] = 0.;
	  part_pyEnd[ntrks_counter] = 0.;
	  part_pzEnd[ntrks_counter] = 0.;
	  part_KEstart[ntrks_counter] = (track_E->at(itrigger)).at(itrack) - (track_M->at(itrigger)).at(itrack);
	  part_KEend[ntrks_counter] = 0.;
	  part_parentid[ntrks_counter] = (track_parent_type->at(itrigger)).at(itrack);
	  part_trackid[ntrks_counter] = itrack;
	  part_pid[ntrks_counter] = (track_ipnu->at(itrigger)).at(itrack);

	  mpid[ntrks_counter] = (track_ipnu->at(itrigger)).at(itrack);
	  px[ntrks_counter] = (track_ux->at(itrigger)).at(itrack);
	  py[ntrks_counter] = (track_uy->at(itrigger)).at(itrack);
	  pz[ntrks_counter] = (track_uz->at(itrigger)).at(itrack);
	  KE[ntrks_counter] = (track_E->at(itrigger)).at(itrack) - (track_M->at(itrigger)).at(itrack);
	  ntrks_counter++;
	}


      }

      ncapturecount = i_ncapturecount;
      neutroncount = i_neutroncount;
      for(int icap=0; icap<i_ncapturecount; icap++){
	capt_x[icap] = i_capt_x[icap];
	capt_y[icap] = i_capt_y[icap];
	capt_z[icap] = i_capt_z[icap];
	capt_t0[icap] = i_capt_t0[icap];
	capt_E[icap] = i_capt_E[icap];
	capt_num[icap] = i_capt_num[icap];
	capt_pid[icap] = i_capt_pid[icap];
	capt_nucleus[icap] = i_capt_nucleus[icap];
	capt_nphot[icap] = i_capt_nphot[icap];
	capt_ngamma[icap] = i_capt_ngamma[icap];
      }

      mode = 0;

      nneutrons = i_nneutrons;

    
      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_tube_id->at(itrigger)).size(); idigitizedhit++){
	// loop on digitized hits in the trigger

	if( idigitizedhit >= MAXNHITS ) {
	  std::clog << " problem: idigitizedhit " << idigitizedhit << " max " << MAXNHITS << std::endl;
	  continue;
	}

	hit_time[idigitizedhit] = (digitized_hit_time->at(itrigger)).at(idigitizedhit);
	if( offset_times_to_zero ){
	  hit_time[idigitizedhit] -= min_time;
	}

	local_tube_id = (digitized_hit_tube_id->at(itrigger)).at(idigitizedhit);
	all_pmts_tree->GetEntry(local_tube_id);
	hit_x[idigitizedhit] = pmt_x*cmTOmm;
	hit_y[idigitizedhit] = pmt_y*cmTOmm;
	hit_z[idigitizedhit] = pmt_z*cmTOmm;
	hit_PMTid[idigitizedhit] = local_tube_id;


      }

    }

    ntrks = ntrks_counter;
    npart = ntrks;


    output_hits_tree.Fill();
  }


  // output pmts tree
  TTree output_pmts_tree("PMTsTree","PMTsTree");
  Int_t pmt_id;
  Double_t pmt_pos_x, pmt_pos_y, pmt_pos_z, pmt_size, pmt_qe, pmt_time_res;
  Int_t is_lappd;
  output_pmts_tree.Branch("pmt_id",&pmt_id,"pmt_id/I");
  output_pmts_tree.Branch("pmt_pos_x",&pmt_pos_x,"pmt_pos_x/D");
  output_pmts_tree.Branch("pmt_pos_y",&pmt_pos_y,"pmt_pos_y/D");
  output_pmts_tree.Branch("pmt_pos_z",&pmt_pos_z,"pmt_pos_z/D");
  output_pmts_tree.Branch("pmt_size",&pmt_size,"pmt_size/D");
  output_pmts_tree.Branch("pmt_qe",&pmt_qe,"pmt_qe/D");
  output_pmts_tree.Branch("pmt_time_res",&pmt_time_res,"pmt_time_res/D");
  output_pmts_tree.Branch("is_lappd",&is_lappd,"is_lappd/I");

  for(int ipmt=0; ipmt<all_pmts_tree->GetEntries(); ipmt++){

    all_pmts_tree->GetEntry(ipmt); 
    pmt_id = pmt_number;
    pmt_pos_x = pmt_x*cmTOmm;
    pmt_pos_y = pmt_y*cmTOmm;
    pmt_pos_z = pmt_z*cmTOmm;
    pmt_size = 304.8;
    pmt_qe = 22.;
    pmt_time_res = 2.;
    is_lappd = -1.718e9;

    output_pmts_tree.Fill();

  }





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


