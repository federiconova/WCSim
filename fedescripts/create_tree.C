//  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WCSIM_BASE_DIR
// C++ includes
#include  <iostream>
#include  <stdlib.h>
#include <iomanip>
#include <cmath>
#include <limits>

//ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "math.h"
#include "TSpectrum2.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH2.h"
#include "TF2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"


//WCSim Includes
#include "../../include/WCSimRootEvent.hh"
#include "../../include/WCSimRootGeom.hh"
#include "../../include/WCSimEnumerations.hh"

#ifdef __MAKECINT__
#pragma link C++ class std::vector<std::vector<Int_t> >+;
#pragma link C++ class std::vector<std::vector<Float_t> >+;
#pragma link C++ class std::vector<std::vector<std::vector<Int_t> > >+;
#pragma link C++ class std::vector<std::vector<std::vector<Float_t> > >+;
#endif

void load_library();
TFile * get_input_file(char *filename);
TString create_filename(const char * prefix, TString& filename_string);

//Definitions
#define PI 3.141592654
#define GridX 20350 // This is the half-width of the grid which will divide the unfolded cylinder 
#define GridY 10000 // This is the half-length of the grid which will divide the unfolded cylinder
#define GridXBin 40 // The number of bins in the horizontal
#define GridYBin 20 // The number of bins in the vertical
/*
*How to run:
*
* enter in the terminal root -l llib.C 'ODAnalysis.C("WCSim.root","outputFile.root",false)' to run the code
* where you replace WCSim.root with your file name and outputfile with the name you wish to save it under
*/

// a structure to hold a pmt id along with its charge
typedef struct {
  double charge;
  double id;
} pmt;
// a structure to hold the details of each grid point
typedef struct {
  double charge;
  int xbin;
  int ybin;
  double x;
  double y;
  double z;
  double totCharge;
} section;

// a function to be used to sort pmt charges using pmt structures
bool CompPMT(pmt a, pmt b)
{
  return a.charge > b.charge;
}
// a function to be used to sort pmt charges using section structures
bool CompSection(section a, section b)
{
  return a.charge > b.charge;
}
// Degrees to Radians conversions
double RadToDeg(double x){
  return x*180/PI;
}

//Radians to Degrees conversions
double DegToRad(double x){
  return x*PI/180;
}
// Print the details of a pmt to the console
void PrintPMT(int pmt_ID1, WCSimRootGeom *geo){

  if (pmt_ID1 <0) exit(1);
    WCSimRootPMT pmt1 = geo->GetPMT(pmt_ID1);

    double pmt1_x = pmt1.GetPosition(0);
    double pmt1_y = pmt1.GetPosition(1);
    double pmt1_z = pmt1.GetPosition(2);

    std::cout << "PMT Position:" << std::endl;
    std::cout << "x: " << pmt1_x << std::endl;
    std::cout << "y: " << pmt1_y << std::endl;
    std::cout << "z: " << pmt1_z << std::endl;

}
// Get the distance between two PMTs
double GetDistance(int pmt_ID1, int pmt_ID2, WCSimRootGeom *geo){

  double distance;

  if ( pmt_ID2 == -1){
    distance = -1;
  }

  else{
    WCSimRootPMT pmt1 = geo->GetPMT(pmt_ID1);
    WCSimRootPMT pmt2 = geo->GetPMT(pmt_ID2);

    double pmt1_x = pmt1.GetPosition(0);
    double pmt1_y = pmt1.GetPosition(1);
    double pmt1_z = pmt1.GetPosition(2);

    double pmt2_x = pmt2.GetPosition(0);
    double pmt2_y = pmt2.GetPosition(1);
    double pmt2_z = pmt2.GetPosition(2);

    distance = sqrt( (pow((pmt2_x - pmt1_x), 2)) + (pow((pmt2_y - pmt1_y), 2)) + (pow((pmt2_z - pmt1_z), 2)) );
    distance = distance / 100. ;
  }

  return distance;
}
// Check if a PMT is within a range of PMT ID numbers, if so, return true.
bool checkPMT(int pmt, int low, int high) {

	if (pmt >= low && pmt <= high) return true;
	else return false;

}
// finds angles using dot product
double Angles(double a[2], double b[2]){

  double angle = 0;
  double dot = a[0]*b[0] + a[1]*b[1];
  double aval = sqrt( pow(a[0],2) + pow(a[1],2) );
  double bval = sqrt( pow(b[0],2) + pow(b[1],2) );

  angle = acos( dot/(aval*bval) );

  return angle;

}

void Reorder(double charges[3], double ids[3]){
  std::vector<pmt> vec;
  pmt a;
  for (int i = 0; i < 3; i++){
    a.charge = charges[i];
    a.id = ids[i];
    vec.push_back(a);
  }

  std::sort(vec.begin(), vec.end(), CompPMT); // sorts highest to lowest

  for (int i = 0; i < 3; i++){
    charges[i] = vec[i].charge;
    ids[i] = vec[i].id;
  }


}

// Convert cylinder coordinates to those on an unfolded plane
void CylinderToSquare(double square[2], int tubeID, double charge, WCSimRootGeom *geo, double radius, double height){

  // Find out where the PMT is in the tank
  double tube[3];
  int cylLoc = geo->GetPMT(tubeID).GetCylLoc();
  tube[0] = geo->GetPMT(tubeID).GetPosition(0);
  tube[1] = geo->GetPMT(tubeID).GetPosition(1);
  tube[2] = geo->GetPMT(tubeID).GetPosition(2);

  //Top (OD || ID)
  if ( cylLoc == 5 || cylLoc == 0){
          square[0] = tube[0];
          square[1] = tube[1] + radius + height/2. + 1.;
  }
  //Bot (OD || ID)
  else if ( cylLoc == 3 || cylLoc == 2){
          square[0] = tube[0];
          square[1] = -(height/2. +radius +tube[1] + 1. );
  }
  //Barrel OD
  else {

          double zero[2] = {0,-1};
          double values[2] = {tube[0],tube[1]};
          double angle = Angles(zero, values);

          double length = angle*radius ;
          if (tube[0]<0) length *= -1;
          square[0] = length;
          square[1] = tube[2];
  }

}
// Convert the coordinates from a 2D unfolded plane back into a cylinder
void SquareToCylinder(double square[2], double cylinder[3], double radius, double height){

  // Find out where the PMT is in the tank

  //Top (OD || ID)
  if ( square[1] >= (height/2. + 1.) ){
          cylinder[0] = square[0];
          cylinder[1] = square[1] - (radius + height/2. + 1.);
          cylinder[2] = height/2.;
  }
  //Bot (OD || ID)
  else if ( square[1] <= - (height/2. + 1.) ){
          cylinder[0] = square[0];
          cylinder[1] = -(square[1] + height/2. +radius + 1.);
          cylinder[2] = -height/2.;
  }
  //Barrel OD
  else {

          double angle = square[0]/radius;
          cylinder[0] = radius*cos( angle - PI/2.);
          cylinder[1] = radius*sin( abs(angle) - PI/2.);
          cylinder[2] = square[1];

  }

}

int FindBin(int numBins, double low, double high, double val){
    double maxlength = high - low;
    double step = maxlength/ (double)numBins;

    double binNum = (val - low)/step;
    return (int) binNum;
}
// Returns an array which contains the 8 neighbours of any grid point in xbin,ybin coordinates
// e.g a[0]is the top right neighbour from the point
void GetNeighbours(int a[8][2], int i, int j){

  a[0][0] = i+1; a[0][1] =  j+1;
  a[1][0] = i+1; a[1][1] =  j-1;
  a[2][0] = i-1; a[2][1] =  j+1;
  a[3][0] = i-1; a[3][1] =  j-1;

  a[4][0] = i; a[4][1] = j+1;
  a[5][0] = i; a[5][1] = j-1;
  a[6][0] = i+1; a[6][1] = j;
  a[7][0] = i-1; a[7][1] = j;

}
// A function to find all of the peaks in the grid, this function looks for all of a point's neighbours to have a lower charge for the point to be considered a peak.
std::vector<section> FindPeaks(double array[GridXBin][GridYBin], double radius, double height){

  std::vector<section> peaks; // first create a vector of section objects

  for (int i = 0; i < GridXBin; i++){
    for (int j = 0; j < GridYBin; j++){
      int larger = 0; // counter to store how many neighbours are lower than the point
      double totcharge = array[i][j]; // initialise charge counter with the charge of the point
      int neighbours[8][2];
      GetNeighbours( neighbours, i, j ); 
      for (int n = 0; n < 8; n++){ // loop through all of the neighbouring sections

        if ( neighbours[n][0] < 0
            || neighbours[n][1] < 0
            || neighbours[n][0] > GridXBin -1
            || neighbours[n][1] > GridYBin -1
          ) continue;
          double val = array[ neighbours[n][0] ][ neighbours[n][1] ]; // find the charge of the neighbouring section
        if (val < array[i][j]) { // if the charge is lower then add it to our point and increase the counter
          totcharge += val;
          larger++;
        }

      } // end of loop over n
      if (larger == 8){ // if all of the neighbours have a lower charge than our point, then we add a peak
        section a;
        a.xbin = i;
        a.ybin = j;
        a.charge = array[i][j];
        a.totCharge = totcharge;
        peaks.push_back(a);
      }

    } // end loop over y (j index)
  }// end loop over x (i index)

  // sort out the clusters from highest to lowest charge
  std::sort(peaks.begin(), peaks.end(), CompSection);

  // find out the distance between clusters
  int i = 0;
  int j = (peaks.size()-1);
  //  std::cout << "PeakSize: " <<  peaks.size() << std::endl;
  if(peaks.size() > 0){

    while (i < (peaks.size()-1) ){
      //std::cout << "Mr hello: " << i << " "<< j <<  std::endl;
      double Val1[2]; double Val2[2];

      Val1[0] = -GridX + ( ((double)peaks[i].xbin + 0.5))*2*GridX/GridXBin;
      Val2[0] = -GridX + ( ((double)peaks[i+1].xbin + 0.5))*2*GridX/GridXBin;
      Val1[1] = -GridY + ( ((double)peaks[i].ybin + 0.5))*2*GridX/GridYBin;
      Val2[1] = -GridY + ( ((double)peaks[i+1].ybin + 0.5))*2*GridX/GridYBin;
      double cluster1[3];
      double cluster2[3];
      SquareToCylinder(Val1, cluster1, radius, height);
      SquareToCylinder(Val2, cluster2, radius, height);

      double difflength = sqrt( pow( (cluster2[0] - cluster1[0]  ) ,2 )  + pow( (cluster2[1] - cluster1[1]  ) ,2 ) + pow( (cluster2[2] - cluster1[2]  ) ,2 )  );

      peaks[i].x = cluster1[0];
      peaks[i].y = cluster1[1];
      peaks[i].z = cluster1[2];

      if (difflength < 1400 ) { // less than the distance of one box corner to corner
        // merge clusters
        peaks[i].charge += peaks[i+1].charge;
        peaks[i].totCharge += peaks[i+1].totCharge; // this double counts some charges but we can remove the double counted charge later
        peaks[i].x = (cluster1[0] + cluster2[0])/2;
        peaks[i].y = (cluster1[1] + cluster2[1])/2;
        peaks[i].z = (cluster1[2] + cluster2[2])/2;

        peaks.erase( peaks.begin() +i + 1 ); // removes the merged element

      }else if(peaks[i+1].totCharge/peaks[0].totCharge < 0.2 ){
        peaks.erase( peaks.begin() +i + 1 ); // removes the small element
      }
      else{i++;}

    }

  }

  return peaks;

}


void loadlibs(){

  char *wcsimdirenv;
  wcsimdirenv = getenv("WCSIMDIR");

  if (wcsimdirenv != NULL) {
    gSystem->Load("$WCSIMDIR/libWCSimRoot.so");
    gSystem->Load("$WCSIMDIR/libWCSimRoot.rootmap");
    gSystem->Load("$WCSIMDIR/libWCSimRootDict_rdict.pcm");
  } else {
    std::cout << "ERROR: WCSIMDIR environment variable has not been set." << std::endl;
  }


}  //end of loadlibs function

int main(int argc, char** argv){


  // load library
  load_library();

  char *inFileName ="../../wcsim.root";

  // open input file
  TFile *inFile = get_input_file(inFileName);

  const char *outFileName = "output.root";
  bool verbosity = 0;

  bool createCanvases = false;

  // Some nice formatting for text options
  std::cout << std::scientific; //all numbers appear in scientific notation
  std::cout << std::setprecision(2); //sets all numbers to output to no more than 2 D.P.
  std::cout << std::left; // Sets the text justification to the left
  const int txtW = 20; // Width of box 'holding' text
  const int numW = 10; // Width of box 'holding' numbers

  int modulo = 100;



	// Get a pointer to the tree from the input file
	TTree *wcsimT = (TTree*) inFile->Get("wcsimT");

	// Get the number of events in the tree
	long int nEvent = wcsimT->GetEntries();
	if (verbosity) { std::cout << "Number of events: "<< nEvent << std::endl;}

	// Create a WCSimRootEvent to put stuff from the tree in
	WCSimRootEvent *wcsimRoot = new WCSimRootEvent();
	WCSimRootEvent *wcsimRootID = new WCSimRootEvent();

	//ID Event readings

	TBranch *IDbranch = wcsimT->GetBranch("wcsimrootevent");
	IDbranch->SetAddress(&wcsimRootID);

	//Force Deletion to prevent memory leak
	wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	//OD Event readings and branch creation

	// Set the branch address for reading from the tree
	TBranch *branch = wcsimT->GetBranch("wcsimrootevent_OD");
	branch->SetAddress(&wcsimRoot);

	// Force deletion to prevent memory leak
	wcsimT->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);

	// Load the geometry tree (only 1 "event")
	TTree* input_geom_tree = (TTree*) inFile->Get("wcsimGeoT");
	WCSimRootGeom *geom = 0;
	input_geom_tree->SetBranchAddress("wcsimrootgeom", &geom);
	if (verbosity) {std::cout << "Geotree has " << input_geom_tree->GetEntries() << " entries." << std::endl;}
	input_geom_tree->GetEntry(0);

	// Start with the main trigger as it always exists and contains most of the info
	WCSimRootTrigger *wcsimTriggerID;
	WCSimRootTrigger *wcsimTriggerOD;

	// Create an output file
	TFile *outFile = new TFile(outFileName, "RECREATE");


	// geometry tree
	TTree geom_tree("geom_tree","geometry tree");
	Float_t detector_length, detector_radius, pmt_radius;
	Int_t number_of_pmts;
	Int_t number_of_pmts_OD;
	geom_tree.Branch("detector_length",&detector_length,"detector_length/F");
	geom_tree.Branch("detector_radius",&detector_radius,"detector_radius/F");
	geom_tree.Branch("pmt_radius",&pmt_radius,"pmt_radius/F");
	geom_tree.Branch("number_of_pmts",&number_of_pmts,"number_of_pmts/I");
	geom_tree.Branch("number_of_pmts_OD",&number_of_pmts_OD,"number_of_pmts_OD/I");
	detector_length = geom->GetWCCylLength();
	detector_radius = geom->GetWCCylRadius();
	pmt_radius = geom->GetWCPMTRadius();
	number_of_pmts = geom->GetWCNumPMT();
	number_of_pmts_OD = geom->GetODWCNumPMT();
	geom_tree.Fill();


  // all pmts tree
  TTree all_pmts_tree("all_pmts_tree","all pmts tree");
  Int_t pmt_number, pmt_location;
  Float_t pmt_ux, pmt_uy, pmt_uz, pmt_x, pmt_y, pmt_z;
  TTree all_pmts_OD_tree("all_pmts_OD_tree","all pmts OD tree");
  Int_t pmt_OD_number, pmt_OD_location;
  Float_t pmt_OD_ux, pmt_OD_uy, pmt_OD_uz, pmt_OD_x, pmt_OD_y, pmt_OD_z;
  all_pmts_tree.Branch("pmt_number",&pmt_number,"pmt_number/I");
  all_pmts_tree.Branch("pmt_location",&pmt_location,"pmt_location/I"); // 0 = top cap, 1 = wall, 2 = bottom cap
  all_pmts_tree.Branch("pmt_ux",&pmt_ux,"pmt_ux/F"); // (ux, uy, uz) direction of the PMT face
  all_pmts_tree.Branch("pmt_uy",&pmt_uy,"pmt_uy/F");
  all_pmts_tree.Branch("pmt_uz",&pmt_uz,"pmt_uz/F");
  all_pmts_tree.Branch("pmt_x",&pmt_x,"pmt_x/F"); // (x, y, z) of the center of the sphere which forms the PMT
  all_pmts_tree.Branch("pmt_y",&pmt_y,"pmt_y/F");
  all_pmts_tree.Branch("pmt_z",&pmt_z,"pmt_z/F");
  all_pmts_OD_tree.Branch("pmt_OD_number",&pmt_OD_number,"pmt_OD_number/I");
  all_pmts_OD_tree.Branch("pmt_OD_location",&pmt_OD_location,"pmt_OD_location/I"); // 0 = top cap, 1 = wall, 2 = bottom cap
  all_pmts_OD_tree.Branch("pmt_OD_ux",&pmt_OD_ux,"pmt_OD_ux/F"); // (ux, uy, uz) direction of the PMT face
  all_pmts_OD_tree.Branch("pmt_OD_uy",&pmt_OD_uy,"pmt_OD_uy/F");
  all_pmts_OD_tree.Branch("pmt_OD_uz",&pmt_OD_uz,"pmt_OD_uz/F");
  all_pmts_OD_tree.Branch("pmt_OD_x",&pmt_OD_x,"pmt_OD_x/F"); // (x, y, z) of the center of the sphere which forms the PMT
  all_pmts_OD_tree.Branch("pmt_OD_y",&pmt_OD_y,"pmt_OD_y/F");
  all_pmts_OD_tree.Branch("pmt_OD_z",&pmt_OD_z,"pmt_OD_z/F");
  for(int i=0; i<number_of_pmts; i++){
    WCSimRootPMT pmt = geom->GetPMT(i);
    pmt_number = pmt.GetTubeNo();
    pmt_location = pmt.GetCylLoc();
    pmt_ux = pmt.GetOrientation(0);
    pmt_uy = pmt.GetOrientation(1);
    pmt_uz = pmt.GetOrientation(2);
    pmt_x = pmt.GetPosition(0);
    pmt_y = pmt.GetPosition(1);
    pmt_z = pmt.GetPosition(2);
    all_pmts_tree.Fill();
  }
  for(int i=0; i<number_of_pmts_OD; i++){
    WCSimRootPMT pmt = geom->GetPMT(number_of_pmts+i);
    pmt_OD_number = pmt.GetTubeNo();
    pmt_OD_location = pmt.GetCylLoc();
    pmt_OD_ux = pmt.GetOrientation(0);
    pmt_OD_uy = pmt.GetOrientation(1);
    pmt_OD_uz = pmt.GetOrientation(2);
    pmt_OD_x = pmt.GetPosition(0);
    pmt_OD_y = pmt.GetPosition(1);
    pmt_OD_z = pmt.GetPosition(2);
    all_pmts_OD_tree.Fill();
  }



  // primary events tree
  TTree primary_events_tree("primary_events_tree","primary events tree");
  WCSimRootTrigger *wcsimrootevent;
  WCSimRootTrigger *wcsimrootevent_OD;

  int n_primary_events = wcsimT->GetEntries();
  int n_primary_events_OD = wcsimT->GetEntries();
  if( n_primary_events != n_primary_events_OD ){
    std::cout << " problem: n_primary_events " << n_primary_events << " n_primary_events_OD " << n_primary_events_OD << std::endl;
    exit(0);
  }
  int number_of_tracks;
  int number_of_raw_cherenkov_hits = 0;
  int number_of_digitized_cherenkov_hits;
  int number_of_raw_cherenkov_hits_OD = 0;
  int number_of_digitized_cherenkov_hits_OD;

  std::vector<Int_t> event_number, trigger_number, trigger_date, trigger_mode, v_trigger_vtxvol, trigger_vec_rec_number, trigger_jmu, trigger_jp, trigger_npar, trigger_ntrack, trigger_number_raw_hits, trigger_number_raw_hits_OD, trigger_number_digitized_hits, trigger_number_digitized_hits_OD, trigger_number_times, trigger_type, trigger_nvertex;
  std::vector<std::vector<Int_t> > trigger_vtxvol;
  std::vector<Float_t> v_trigger_vtx_x, v_trigger_vtx_y, v_trigger_vtx_z, trigger_sum_q;
  std::vector<std::vector<Float_t> > trigger_info, trigger_vtx_x, trigger_vtx_y, trigger_vtx_z;

  std::vector<Int_t> v_track_ipnu, v_track_parent_type, v_track_flag, v_track_start_volume, v_track_stop_volume, v_track_id;
  std::vector<std::vector<Int_t> > track_ipnu,   track_parent_type,   track_flag,   track_start_volume,   track_stop_volume,   track_id;

  std::vector<Float_t> v_track_ux, v_track_uy, v_track_uz, v_track_M, v_track_P, v_track_E, v_track_px, v_track_py, v_track_pz, v_track_stop_x, v_track_stop_y, v_track_stop_z, v_track_start_x, v_track_start_y, v_track_start_z, v_track_time;
  std::vector<std::vector<Float_t> > track_ux,   track_uy,   track_uz,   track_M,   track_P,   track_E,   track_px,   track_py,   track_pz,   track_stop_x,   track_stop_y,   track_stop_z,   track_start_x,   track_start_y,   track_start_z,   track_time;

  std::vector<Int_t> v_raw_hit_tube_id, v_raw_hit_tube_times_indexes, v_raw_hit_tube_pe;
  std::vector<std::vector<Int_t> > raw_hit_tube_id, raw_hit_tube_times_indexes, raw_hit_tube_pe;
  std::vector<Int_t> v_raw_hit_OD_tube_id, v_raw_hit_OD_tube_times_indexes, v_raw_hit_OD_tube_pe;
  std::vector<std::vector<Int_t> > raw_hit_OD_tube_id, raw_hit_OD_tube_times_indexes, raw_hit_OD_tube_pe;

  std::vector<Float_t> vv_raw_hit_times;
  std::vector<std::vector<Float_t> > v_raw_hit_times;
  std::vector<std::vector<std::vector<Float_t> > > raw_hit_times;
  std::vector<Int_t> vv_raw_hit_parent_ids;
  std::vector<std::vector<Int_t> > v_raw_hit_parent_ids;
  std::vector<std::vector<std::vector<Int_t> > > raw_hit_parent_ids;
  std::vector<Float_t> vv_raw_hit_OD_times;
  std::vector<std::vector<Float_t> > v_raw_hit_OD_times;
  std::vector<std::vector<std::vector<Float_t> > > raw_hit_OD_times;
  std::vector<Int_t> vv_raw_hit_OD_parent_ids;
  std::vector<std::vector<Int_t> > v_raw_hit_OD_parent_ids;
  std::vector<std::vector<std::vector<Int_t> > > raw_hit_OD_parent_ids;

  std::vector<Int_t> v_digitized_hit_tube_id;
  std::vector<std::vector<Int_t> > digitized_hit_tube_id;
  std::vector<Float_t> v_digitized_hit_Q, v_digitized_hit_time;
  std::vector<std::vector<Float_t> > digitized_hit_Q, digitized_hit_time;
  std::vector<Int_t> vv_digitized_hit_photon_ids;
  std::vector<std::vector<Int_t> > v_digitized_hit_photon_ids;
  std::vector<std::vector<std::vector<Int_t> > > digitized_hit_photon_ids;

  std::vector<Int_t> v_digitized_hit_OD_tube_id;
  std::vector<std::vector<Int_t> > digitized_hit_OD_tube_id;
  std::vector<Float_t> v_digitized_hit_OD_Q, v_digitized_hit_OD_time;
  std::vector<std::vector<Float_t> > digitized_hit_OD_Q, digitized_hit_OD_time;
  std::vector<Int_t> vv_digitized_hit_OD_photon_ids;
  std::vector<std::vector<Int_t> > v_digitized_hit_OD_photon_ids;
  std::vector<std::vector<std::vector<Int_t> > > digitized_hit_OD_photon_ids;

  primary_events_tree.Branch("event_number",&event_number);
  primary_events_tree.Branch("trigger_number",&trigger_number);
  primary_events_tree.Branch("trigger_date",&trigger_date);
  primary_events_tree.Branch("trigger_mode",&trigger_mode); // interaction mode
  primary_events_tree.Branch("trigger_vtxvol",&trigger_vtxvol); // volume of vertex
  primary_events_tree.Branch("trigger_vtx_x",&trigger_vtx_x); // interaction vertex
  primary_events_tree.Branch("trigger_vtx_y",&trigger_vtx_y);
  primary_events_tree.Branch("trigger_vtx_z",&trigger_vtx_z);
  primary_events_tree.Branch("trigger_vec_rec_number",&trigger_vec_rec_number); // info event number in inputvetcotfile
  primary_events_tree.Branch("trigger_jmu",&trigger_jmu); // index to muon
  primary_events_tree.Branch("trigger_jp",&trigger_jp); // index to proton
  primary_events_tree.Branch("trigger_npar",&trigger_npar); // number of final state particles
  primary_events_tree.Branch("trigger_ntrack",&trigger_ntrack);
  primary_events_tree.Branch("trigger_nvertex",&trigger_nvertex);
  primary_events_tree.Branch("trigger_number_raw_hits",&trigger_number_raw_hits); // Total number of tubes with hits
  primary_events_tree.Branch("trigger_number_raw_hits_OD",&trigger_number_raw_hits_OD); 
  primary_events_tree.Branch("trigger_number_digitized_hits",&trigger_number_digitized_hits); // Number of PMTs with digitized hits
  primary_events_tree.Branch("trigger_number_digitized_hits_OD",&trigger_number_digitized_hits_OD); // Number of PMTs with digitized hits
  primary_events_tree.Branch("trigger_sum_q",&trigger_sum_q); // sum of q(readout digitized pe) in event
  primary_events_tree.Branch("trigger_type",&trigger_type); // enumeration of trigger type (stored as int)
  primary_events_tree.Branch("trigger_info",&trigger_info); // info about why the trigger was passed
  primary_events_tree.Branch("trigger_number_times",&trigger_number_times);

  primary_events_tree.Branch("track_ipnu",&track_ipnu); // id of final state particle
  primary_events_tree.Branch("track_parent_type",&track_parent_type); // ID of parent of ith particle (0 if primary)
  primary_events_tree.Branch("track_ux",&track_ux); // direction of ith final state particle
  primary_events_tree.Branch("track_uy",&track_uy);
  primary_events_tree.Branch("track_uz",&track_uz);
  primary_events_tree.Branch("track_px",&track_px); // momentum of ith final state particle
  primary_events_tree.Branch("track_py",&track_py);
  primary_events_tree.Branch("track_pz",&track_pz);
  primary_events_tree.Branch("track_flag",&track_flag); // flag: -1 = incoming neutrino
  //       -2 = target
  //        1 = outgoing lepton
  //        2 = most energetic outgoing nucleon
  primary_events_tree.Branch("track_M",&track_M); // mass of ith final state particle
  primary_events_tree.Branch("track_P",&track_P); // momentum of ith final state particle
  primary_events_tree.Branch("track_E",&track_E); // energy of ith final state particle
  primary_events_tree.Branch("track_start_volume",&track_start_volume); // starting volume of ith final state particle
  primary_events_tree.Branch("track_stop_volume",&track_stop_volume); // stopping volume of ith final state particle
  primary_events_tree.Branch("track_stop_x",&track_stop_x); // stopping point of ith final state particle 
  primary_events_tree.Branch("track_stop_y",&track_stop_y);
  primary_events_tree.Branch("track_stop_z",&track_stop_z);
  primary_events_tree.Branch("track_start_x",&track_start_x); // starting point of ith final state particle
  primary_events_tree.Branch("track_start_y",&track_start_y);
  primary_events_tree.Branch("track_start_z",&track_start_z);
  primary_events_tree.Branch("track_time",&track_time); // creation time of ith final state particle
  primary_events_tree.Branch("track_id",&track_id);

  primary_events_tree.Branch("raw_hit_tube_id",&raw_hit_tube_id);
  primary_events_tree.Branch("raw_hit_tube_times_indexes",&raw_hit_tube_times_indexes); // The indexes of times recorded at each tube
  primary_events_tree.Branch("raw_hit_tube_pe",&raw_hit_tube_pe); // The totalPE recorded at each tube
  primary_events_tree.Branch("raw_hit_times",&raw_hit_times); // The time of each hit of each tube
  primary_events_tree.Branch("raw_hit_parent_ids",&raw_hit_parent_ids); // The parent id of each hit of each tube

  primary_events_tree.Branch("raw_hit_OD_tube_id",&raw_hit_OD_tube_id);
  primary_events_tree.Branch("raw_hit_OD_tube_times_indexes",&raw_hit_OD_tube_times_indexes); // The indexes of times recorded at each tube
  primary_events_tree.Branch("raw_hit_OD_tube_pe",&raw_hit_OD_tube_pe); // The totalPE recorded at each tube
  primary_events_tree.Branch("raw_hit_OD_times",&raw_hit_OD_times); // The time of each hit_OD of each tube
  primary_events_tree.Branch("raw_hit_OD_parent_ids",&raw_hit_OD_parent_ids); // The parent id of each hit_OD of each tube

  primary_events_tree.Branch("digitized_hit_tube_id",&digitized_hit_tube_id); // The readout tube ID
  primary_events_tree.Branch("digitized_hit_Q",&digitized_hit_Q); // The readout digitized pe
  primary_events_tree.Branch("digitized_hit_time",&digitized_hit_time);  // The readout digitized time
  primary_events_tree.Branch("digitized_hit_photon_ids",&digitized_hit_photon_ids);

  primary_events_tree.Branch("digitized_hit_OD_tube_id",&digitized_hit_OD_tube_id); // The readout tube ID
  primary_events_tree.Branch("digitized_hit_OD_Q",&digitized_hit_OD_Q); // The readout digitized pe
  primary_events_tree.Branch("digitized_hit_OD_time",&digitized_hit_OD_time);  // The readout digitized time
  primary_events_tree.Branch("digitized_hit_OD_photon_ids",&digitized_hit_OD_photon_ids);

  for(int ievent=0; ievent<n_primary_events; ievent++){ 

    // loop on primary events
    wcsimT->GetEvent(ievent); 

    // init trigger variables
    event_number.clear();
    trigger_number.clear();
    trigger_date.clear();
    trigger_mode.clear();
    trigger_vtxvol.clear();
    trigger_vtx_x.clear();
    trigger_vtx_y.clear();
    trigger_vtx_z.clear();
    trigger_vec_rec_number.clear();
    trigger_jmu.clear();
    trigger_jp.clear();
    trigger_npar.clear();
    trigger_ntrack.clear();
    trigger_nvertex.clear();
    trigger_number_raw_hits.clear();
    trigger_number_raw_hits_OD.clear();
    trigger_number_digitized_hits.clear();
    trigger_number_digitized_hits_OD.clear();
    trigger_number_times.clear();
    trigger_sum_q.clear();
    trigger_type.clear();
    trigger_info.clear();

    // init track variables
    track_ipnu.clear();
    track_parent_type.clear();
    track_ux.clear();
    track_uy.clear();
    track_uz.clear();
    track_px.clear();
    track_py.clear();
    track_pz.clear();
    track_flag.clear();
    track_M.clear();
    track_P.clear();
    track_E.clear();
    track_start_volume.clear();
    track_stop_volume.clear();
    track_stop_x.clear();
    track_stop_y.clear();
    track_stop_z.clear();
    track_start_x.clear();
    track_start_y.clear();
    track_start_z.clear();
    track_time.clear();
    track_id.clear();

    // init raw hits variables
    raw_hit_tube_id.clear();
    raw_hit_tube_times_indexes.clear();
    raw_hit_tube_pe.clear();
    raw_hit_times.clear();
    raw_hit_parent_ids.clear();

    raw_hit_OD_tube_id.clear();
    raw_hit_OD_tube_times_indexes.clear();
    raw_hit_OD_tube_pe.clear();
    raw_hit_OD_times.clear();
    raw_hit_OD_parent_ids.clear();

    // init digitized hits variables
    digitized_hit_tube_id.clear();
    digitized_hit_Q.clear();
    digitized_hit_time.clear();
    digitized_hit_photon_ids.clear();

    digitized_hit_OD_tube_id.clear();
    digitized_hit_OD_Q.clear();
    digitized_hit_OD_time.clear();
    digitized_hit_OD_photon_ids.clear();
  
    if( wcsimRootID->GetNumberOfEvents() != wcsimRoot->GetNumberOfEvents() ){
      std::cout << " problem: n events " << wcsimRootID->GetNumberOfEvents() << " OD " << wcsimRoot->GetNumberOfEvents() << std::endl;
      exit(0);
    }

    for(int itrigger=0; itrigger<wcsimRootID->GetNumberOfEvents(); itrigger++){
      // loop on triggers in the event

      if( itrigger > 0 ) continue; // seg faults

      wcsimrootevent = wcsimRootID->GetTrigger(itrigger);
      wcsimrootevent_OD = wcsimRoot->GetTrigger(itrigger);

      event_number.push_back(wcsimrootevent->GetHeader()->GetEvtNum());
      trigger_number.push_back(itrigger); 
      trigger_date.push_back(wcsimrootevent->GetHeader()->GetDate());
      trigger_mode.push_back(wcsimrootevent->GetMode());

      // init trigger vertex variables
      v_trigger_vtxvol.clear();
      v_trigger_vtx_x.clear();
      v_trigger_vtx_y.clear();
      v_trigger_vtx_z.clear();

      int n_vertexes = 1;
      trigger_nvertex.push_back(n_vertexes);
      v_trigger_vtxvol.push_back(wcsimrootevent->GetVtxvol());
      v_trigger_vtx_x.push_back(wcsimrootevent->GetVtx(0));
      v_trigger_vtx_y.push_back(wcsimrootevent->GetVtx(1));
      v_trigger_vtx_z.push_back(wcsimrootevent->GetVtx(2));

      trigger_vtxvol.push_back(v_trigger_vtxvol);
      trigger_vtx_x.push_back(v_trigger_vtx_x);
      trigger_vtx_y.push_back(v_trigger_vtx_y);
      trigger_vtx_z.push_back(v_trigger_vtx_z);

      trigger_vec_rec_number.push_back(wcsimrootevent->GetVecRecNumber());
      trigger_jmu.push_back(wcsimrootevent->GetJmu());
      trigger_jp.push_back(wcsimrootevent->GetJp());
      trigger_npar.push_back(wcsimrootevent->GetNpar());
      number_of_tracks = wcsimrootevent->GetNtrack();
      trigger_ntrack.push_back(number_of_tracks);
      number_of_raw_cherenkov_hits = wcsimrootevent->GetNumTubesHit();
      trigger_number_raw_hits.push_back(number_of_raw_cherenkov_hits);
      number_of_raw_cherenkov_hits_OD = wcsimrootevent_OD->GetNumTubesHit();
      trigger_number_raw_hits_OD.push_back(number_of_raw_cherenkov_hits_OD);
      number_of_digitized_cherenkov_hits = wcsimrootevent->GetNcherenkovdigihits();
      number_of_digitized_cherenkov_hits_OD = wcsimrootevent_OD->GetNcherenkovdigihits();
      trigger_number_digitized_hits.push_back(number_of_digitized_cherenkov_hits);
      trigger_number_digitized_hits_OD.push_back(number_of_digitized_cherenkov_hits_OD);
      trigger_number_times.push_back(wcsimrootevent->GetNcherenkovhittimes());
      trigger_sum_q.push_back(wcsimrootevent->GetSumQ());

      trigger_type.push_back(wcsimrootevent->GetTriggerType());
      trigger_info.push_back(wcsimrootevent->GetTriggerInfo());

      // init trigger track variables
      v_track_ipnu.clear();
      v_track_parent_type.clear();
      v_track_flag.clear();
      v_track_start_volume.clear();
      v_track_stop_volume.clear();
      v_track_id.clear();
      v_track_ux.clear();
      v_track_uy.clear();
      v_track_uz.clear();
      v_track_M.clear();
      v_track_P.clear();
      v_track_E.clear();
      v_track_px.clear();
      v_track_py.clear();
      v_track_pz.clear();
      v_track_stop_x.clear();
      v_track_stop_y.clear();
      v_track_stop_z.clear();
      v_track_start_x.clear();  
      v_track_start_y.clear();
      v_track_start_z.clear();
      v_track_time.clear();

      for(int itrack=0; itrack<number_of_tracks; itrack++){

	// loop on tracks in the trigger
	TObject *element = (wcsimrootevent->GetTracks())->At(itrack);
	WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);

	v_track_ipnu.push_back(wcsimroottrack->GetIpnu());
	v_track_parent_type.push_back(wcsimroottrack->GetParenttype());
	v_track_ux.push_back(wcsimroottrack->GetDir(0));
	v_track_uy.push_back(wcsimroottrack->GetDir(1));
	v_track_uz.push_back(wcsimroottrack->GetDir(2));
	v_track_px.push_back(wcsimroottrack->GetPdir(0));
	v_track_py.push_back(wcsimroottrack->GetPdir(1));
	v_track_pz.push_back(wcsimroottrack->GetPdir(2));
	v_track_flag.push_back(wcsimroottrack->GetFlag());
	v_track_M.push_back(wcsimroottrack->GetM());
	v_track_P.push_back(wcsimroottrack->GetP());
	v_track_E.push_back(wcsimroottrack->GetE());
	v_track_start_volume.push_back(wcsimroottrack->GetStartvol());
	v_track_stop_volume.push_back(wcsimroottrack->GetStopvol());
	v_track_stop_x.push_back(wcsimroottrack->GetStop(0));
	v_track_stop_y.push_back(wcsimroottrack->GetStop(1));
	v_track_stop_z.push_back(wcsimroottrack->GetStop(2));
	v_track_start_x.push_back(wcsimroottrack->GetStart(0));
	v_track_start_y.push_back(wcsimroottrack->GetStart(1));
	v_track_start_z.push_back(wcsimroottrack->GetStart(2));
	v_track_time.push_back(wcsimroottrack->GetTime());
	v_track_id.push_back(wcsimroottrack->GetId());
      }
      
      track_ipnu.push_back(v_track_ipnu);
      track_parent_type.push_back(v_track_parent_type);
      track_ux.push_back(v_track_ux);
      track_uy.push_back(v_track_uy);
      track_uz.push_back(v_track_uz);
      track_px.push_back(v_track_px);
      track_py.push_back(v_track_py);
      track_pz.push_back(v_track_pz);
      track_flag.push_back(v_track_flag);
      track_M.push_back(v_track_M);
      track_P.push_back(v_track_P);
      track_E.push_back(v_track_E);
      track_start_volume.push_back(v_track_start_volume);
      track_stop_volume.push_back(v_track_stop_volume);
      track_stop_x.push_back(v_track_stop_x);
      track_stop_y.push_back(v_track_stop_y);
      track_stop_z.push_back(v_track_stop_z);
      track_start_x.push_back(v_track_start_x);
      track_start_y.push_back(v_track_start_y);
      track_start_z.push_back(v_track_start_z);
      track_time.push_back(v_track_time);
      track_id.push_back(v_track_id);
        
      // init trigger raw hits variables
      v_raw_hit_tube_id.clear();
      v_raw_hit_tube_times_indexes.clear();
      v_raw_hit_tube_pe.clear();
      v_raw_hit_times.clear();
      v_raw_hit_parent_ids.clear();

      v_raw_hit_OD_tube_id.clear();
      v_raw_hit_OD_tube_times_indexes.clear();
      v_raw_hit_OD_tube_pe.clear();
      v_raw_hit_OD_times.clear();
      v_raw_hit_OD_parent_ids.clear();

      // init trigger digitized hits variables
      v_digitized_hit_tube_id.clear();
      v_digitized_hit_Q.clear();
      v_digitized_hit_time.clear();
      v_digitized_hit_photon_ids.clear();

      v_digitized_hit_OD_tube_id.clear();
      v_digitized_hit_OD_Q.clear();
      v_digitized_hit_OD_time.clear();
      v_digitized_hit_OD_photon_ids.clear();

      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();
      TClonesArray *timeArray_OD = wcsimrootevent_OD->GetCherenkovHitTimes();
    
      for(int irawhit=0; irawhit<number_of_raw_cherenkov_hits; irawhit++){

	// loop on raw hits in the trigger
	TObject *Hit = (wcsimrootevent->GetCherenkovHits())->At(irawhit);
	WCSimRootCherenkovHit *wcsimrootcherenkovhit = dynamic_cast<WCSimRootCherenkovHit*>(Hit);

	int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	v_raw_hit_tube_id.push_back(tubeNumber);
	v_raw_hit_tube_times_indexes.push_back(timeArrayIndex);
	v_raw_hit_tube_pe.push_back(peForTube);

	// extract the PMT object
	//WCSimRootPMT pmt   = geom->GetPMT(tubeNumber-1);

	// init trigger raw hit time variables
	vv_raw_hit_times.clear();
	vv_raw_hit_parent_ids.clear();

	for (int irawhittimes = timeArrayIndex; irawhittimes < timeArrayIndex + peForTube; irawhittimes++)
	  {
	    WCSimRootCherenkovHitTime * HitTime = dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(irawhittimes));
	    
	    vv_raw_hit_times.push_back(HitTime->GetTruetime());
	    vv_raw_hit_parent_ids.push_back(HitTime->GetParentID());
	  }
	
	v_raw_hit_times.push_back(vv_raw_hit_times);
	v_raw_hit_parent_ids.push_back(vv_raw_hit_parent_ids);

      } // End of loop over Cherenkov hits

      raw_hit_tube_id.push_back(v_raw_hit_tube_id);
      raw_hit_tube_times_indexes.push_back(v_raw_hit_tube_times_indexes);
      raw_hit_tube_pe.push_back(v_raw_hit_tube_pe);
      raw_hit_times.push_back(v_raw_hit_times);
      raw_hit_parent_ids.push_back(v_raw_hit_parent_ids);


      for(int irawhit=0; irawhit<number_of_raw_cherenkov_hits_OD; irawhit++){

	// loop on raw hits in the trigger
	TObject *Hit = (wcsimrootevent_OD->GetCherenkovHits())->At(irawhit);
	WCSimRootCherenkovHit *wcsimrootcherenkovhit = dynamic_cast<WCSimRootCherenkovHit*>(Hit);

	int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	v_raw_hit_OD_tube_id.push_back(tubeNumber);
	v_raw_hit_OD_tube_times_indexes.push_back(timeArrayIndex);
	v_raw_hit_OD_tube_pe.push_back(peForTube);

	// extract the PMT object
	//WCSimRootPMT pmt   = geom->GetPMT(tubeNumber-1);

	// init trigger raw hit time variables
	vv_raw_hit_OD_times.clear();
	vv_raw_hit_OD_parent_ids.clear();

	for (int irawhittimes = timeArrayIndex; irawhittimes < timeArrayIndex + peForTube; irawhittimes++)
	  {
	    WCSimRootCherenkovHitTime * HitTime = dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray_OD->At(irawhittimes));
	    
	    vv_raw_hit_OD_times.push_back(HitTime->GetTruetime());
	    vv_raw_hit_OD_parent_ids.push_back(HitTime->GetParentID());
	  }
	
	v_raw_hit_OD_times.push_back(vv_raw_hit_OD_times);
	v_raw_hit_OD_parent_ids.push_back(vv_raw_hit_OD_parent_ids);

      } // End of loop over Cherenkov hits

      raw_hit_OD_tube_id.push_back(v_raw_hit_OD_tube_id);
      raw_hit_OD_tube_times_indexes.push_back(v_raw_hit_OD_tube_times_indexes);
      raw_hit_OD_tube_pe.push_back(v_raw_hit_OD_tube_pe);
      raw_hit_OD_times.push_back(v_raw_hit_OD_times);
      raw_hit_OD_parent_ids.push_back(v_raw_hit_OD_parent_ids);

      // init trigger digitized hits variables
      v_digitized_hit_tube_id.clear();
      v_digitized_hit_Q.clear();
      v_digitized_hit_time.clear();

      for(int idigitizedhit=0; idigitizedhit<number_of_digitized_cherenkov_hits; idigitizedhit++){

	// loop on digitized hits in the trigger
	TObject *Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(idigitizedhit);
    	WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);

	v_digitized_hit_tube_id.push_back(wcsimrootcherenkovdigihit->GetTubeId());
	v_digitized_hit_Q.push_back(wcsimrootcherenkovdigihit->GetQ());
	v_digitized_hit_time.push_back(wcsimrootcherenkovdigihit->GetT());

	// init trigger raw hit time variables
	vv_digitized_hit_photon_ids.clear();
	wcsimrootcherenkovdigihit->GetPhotonIds();
	for(unsigned int iph=0; iph<(wcsimrootcherenkovdigihit->GetPhotonIds()).size(); iph++)
	  vv_digitized_hit_photon_ids.push_back((wcsimrootcherenkovdigihit->GetPhotonIds()).at(iph));

	v_digitized_hit_photon_ids.push_back(vv_digitized_hit_photon_ids);

	//std::clog << " ievent " << ievent << " itrigger " << itrigger << " idigitizedhit " << idigitizedhit << " tubeid " << wcsimrootcherenkovdigihit->GetTubeId() << " Q " << wcsimrootcherenkovdigihit->GetQ() << " time " << wcsimrootcherenkovdigihit->GetT() << std::endl;

      }

      digitized_hit_tube_id.push_back(v_digitized_hit_tube_id);
      digitized_hit_Q.push_back(v_digitized_hit_Q);
      digitized_hit_time.push_back(v_digitized_hit_time);
      digitized_hit_photon_ids.push_back(v_digitized_hit_photon_ids);


      // init trigger digitized hits variables
      v_digitized_hit_OD_tube_id.clear();
      v_digitized_hit_OD_Q.clear();
      v_digitized_hit_OD_time.clear();

      for(int idigitizedhit=0; idigitizedhit<number_of_digitized_cherenkov_hits_OD; idigitizedhit++){

	// loop on digitized hits in the trigger
	TObject *Hit = (wcsimrootevent_OD->GetCherenkovDigiHits())->At(idigitizedhit);
    	WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);

	v_digitized_hit_OD_tube_id.push_back(wcsimrootcherenkovdigihit->GetTubeId());
	v_digitized_hit_OD_Q.push_back(wcsimrootcherenkovdigihit->GetQ());
	v_digitized_hit_OD_time.push_back(wcsimrootcherenkovdigihit->GetT());

	// init trigger raw hit time variables
	vv_digitized_hit_OD_photon_ids.clear();
	wcsimrootcherenkovdigihit->GetPhotonIds();
	for(unsigned int iph=0; iph<(wcsimrootcherenkovdigihit->GetPhotonIds()).size(); iph++)
	  vv_digitized_hit_OD_photon_ids.push_back((wcsimrootcherenkovdigihit->GetPhotonIds()).at(iph));

	v_digitized_hit_OD_photon_ids.push_back(vv_digitized_hit_OD_photon_ids);

      }


      digitized_hit_OD_tube_id.push_back(v_digitized_hit_OD_tube_id);
      digitized_hit_OD_Q.push_back(v_digitized_hit_OD_Q);
      digitized_hit_OD_time.push_back(v_digitized_hit_OD_time);
      digitized_hit_OD_photon_ids.push_back(v_digitized_hit_OD_photon_ids);

      if( ievent%modulo == 0 )
	std::clog << " ievent " << ievent << " of " << n_primary_events << " itrigger " << itrigger << " number_of_raw_cherenkov_hits " << number_of_raw_cherenkov_hits << " number_of_raw_cherenkov_hits_OD " << number_of_raw_cherenkov_hits_OD << std::endl;
    
    
    }

    primary_events_tree.Fill();

  }


	TTree *outBranch = new TTree("simulation", "simulation");

  TH1D *clusterHist = new TH1D("clusterHist", "clusterHist", 21, 0, 20);

	// Detector Geometry Details

	int MAXPMT = geom->GetWCNumPMT(); //Get the maximum number of PMTs in the ID
	int MAXPMTA = geom->GetODWCNumPMT(); //Get the maximum number of PMTs in the OD

	bool idOn = true; // Boolean to keep track of whether or not the ID was constructed (sometimes turned off for speed) true = on , false = off
	bool odOn = true; // Boolean to keep track of whether or not the ID was constructed (sometimes turned off for speed) true = on , false = off
	if (MAXPMT == 0 ) {idOn = false;}
	if (MAXPMTA == 0 ) {odOn = false;}

	double RadiusID = 0;
	double RadiusOD = 0;
	double HeightID = 0;
	double HeightOD = 0;

	// Find a barrel pmt
	bool barrelID = false; // Boolean to see if a barrel ID pmt has been found
	bool barrelOD = false;  // Boolean to see if a barrel OD pmt has been found
	bool capID = false;  // Boolean to see if a cap ID pmt has been found
	bool capOD = false;  // Boolean to see if a cap OD pmt has been found
	int barrelPMTID = -1; // PMT number of the barrel ID PMT
	int barrelPMTOD = -1; // PMT number of the barrel OD PMT
	int capPMTID = -1; // PMT number of the cap ID PMT
	int capPMTOD = -1; // PMT number of the cap OD PMT
	int pmtCount = 0; // Number used to count through all of the PMTs

	if (!idOn) {barrelID = true; capID = true;} // If no ID constructed, don't look for ID PMTs
	if (!odOn) {barrelOD = true; capOD = true;} // If no OD constructed, don't look for OD PMTs

	while (!barrelID || !barrelOD || !capID || !capOD ){ // Loop to look for barrel and cap PMTs to work out the radius and height respectively.

		//std::cout <<"Looking for barrel ID PMT: " << geom->GetPMT(pmtCount).GetCylLoc()<<std::endl;
		if ( !barrelID && ( geom->GetPMT(pmtCount).GetCylLoc() == 1)  ) {barrelID = true; barrelPMTID = pmtCount; }
		if ( !barrelOD && (geom->GetPMT(pmtCount).GetCylLoc() == 4 )  ) {barrelOD = true; barrelPMTOD = pmtCount; }
		if ( !capID && (geom->GetPMT(pmtCount).GetCylLoc() == 0 || geom->GetPMT(pmtCount).GetCylLoc() == 2)  ) {capID = true; capPMTID = pmtCount; }
		if ( !capOD && (geom->GetPMT(pmtCount).GetCylLoc() == 3 || geom->GetPMT(pmtCount).GetCylLoc() == 5)  ) {capOD = true; capPMTOD = pmtCount; }
		pmtCount++;
		//pmtCount+= 10; // Can speed up this process by checking PMTs in multiples higher than 1
	}


	if (idOn) { // If ID is on, check the PMTs are correct and set the height and radius.
		if (checkPMT(barrelPMTID, 0, MAXPMT -1 ) && checkPMT(capPMTID, 0, MAXPMT -1) ) {
			// Set the radius and height of the ID using the PMTs' positions.
			RadiusID = sqrt( pow( geom->GetPMT(barrelPMTID).GetPosition(0),2) + pow(geom->GetPMT(barrelPMTID).GetPosition(1),2) );
			HeightID = 2*(abs(geom->GetPMT(capPMTID).GetPosition(2)));

		}
		else {
			std::cerr << "Can not understand the tank geometry. Exiting..." << std::endl;
			exit(1);
		}
	}

	if (odOn) { // If OD is on, check the PMTs are correct and set the height and radius.
		if (checkPMT(barrelPMTOD, MAXPMT, MAXPMT + MAXPMTA-1) && checkPMT(capPMTOD, MAXPMT, MAXPMT + MAXPMTA - 1) ) {
			// Set the radius and height of the ID and OD using the PMTs' positions.
			RadiusOD = sqrt( pow( geom->GetPMT(barrelPMTOD).GetPosition(0),2) + pow(geom->GetPMT(barrelPMTOD).GetPosition(1),2) );
			HeightOD = 2*(abs(geom->GetPMT(capPMTOD).GetPosition(2)));

		}
		else {
			std::cerr << "Can not understand the tank geometry. Exiting..." << std::endl;
			exit(1);
		}
	}


	if (idOn) {
		std::cout << "Barrel Radius (ID) is: " << RadiusID <<std::endl;
		std::cout << "Barrel Height (ID) " <<  HeightID <<std::endl;
	}
	if (odOn) {
		std::cout << "Barrel Radius (OD) is: " << RadiusOD <<std::endl;
		std::cout << "Barrel Height (OD) " <<  HeightOD <<std::endl;
	}


	// OD Event Analysis
  for (int ev = 0; ev < nEvent; ev++){ // Loop over events
  //for (int ev = 21; ev < 22; ev++){ // Loop over events
	  wcsimT->GetEntry(ev);
	  wcsimTriggerOD = wcsimRoot->GetTrigger(0);
	  int numTriggers = wcsimRoot->GetNumberOfEvents();
	  int numSubTriggers = wcsimRoot->GetNumberOfSubEvents();

	  for (int nTrig = 0; nTrig < numTriggers; nTrig++){

	    wcsimTriggerOD = wcsimRoot->GetTrigger(nTrig);
	    int numTracks = wcsimTriggerOD->GetNtrack();

	    if ( numTracks != 0){
        WCSimRootTrack * trackOD = (WCSimRootTrack*) wcsimTriggerOD->GetTracks()->At(0);

        double tankArray[GridXBin][GridYBin] = {}; // split tank into GridXBin by GridYBin

	      double vtxX = wcsimTriggerOD->GetVtx(0);
	      double vtxY = wcsimTriggerOD->GetVtx(1);
	      double vtxZ = wcsimTriggerOD->GetVtx(2);
	      double dirX = trackOD->GetDir(0);
	      double dirY = trackOD->GetDir(1);
	      double dirZ = trackOD->GetDir(2);

	      double energy = trackOD->GetE();
	      double rawHits = 0;
	      double digiHits = 0;
        double largestHit = -1;
        double currentDist = 9999999;
        int largestHitID = -1;
        int largestHitID2 = -1;

	      int numPMTsHit = wcsimTriggerOD->GetNcherenkovhits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
	      int numPMTsDigiHit = wcsimTriggerOD->GetNcherenkovdigihits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
        if (verbosity) {std::cout << "Number of OD Hits: " << numPMTsDigiHit << std::endl;}
        if (verbosity) {std::cout << "Number triggers: " << nTrig << std::endl;}

    		for (int i = 0; i < numPMTsDigiHit; i++){

    		  WCSimRootCherenkovDigiHit *cherenkovDigiHit = (WCSimRootCherenkovDigiHit*) wcsimTriggerOD->GetCherenkovDigiHits()->At(i);
    		  double tmpDigiHits = cherenkovDigiHit->GetQ();
          int tubeID = cherenkovDigiHit->GetTubeId();
          int ODtubeID = tubeID + MAXPMT - 1;

          double tube[3];
          int cylLoc = geom->GetPMT(ODtubeID).GetCylLoc();
          tube[0] = geom->GetPMT(ODtubeID).GetPosition(0);
          tube[1] = geom->GetPMT(ODtubeID).GetPosition(1);
          tube[2] = geom->GetPMT(ODtubeID).GetPosition(2);

          double square[2] = {};
          double cylinder[3] = {};
          CylinderToSquare(square, ODtubeID, tmpDigiHits, geom, RadiusOD, HeightOD);
          tankArray[FindBin(GridXBin, -GridX, GridX, square[0])][FindBin(GridYBin, -GridY, GridY, square[1])] += tmpDigiHits;

    		} // End of loop over digi hits


        std::vector<section> peaks = FindPeaks(tankArray, RadiusOD, HeightOD);
        if (verbosity) {std::cout << "Number of Peaks found: " << peaks.size() << std::endl;}
        clusterHist->Fill(peaks.size());

	    }

	  } // End of loop over triggers

	}


	if (createCanvases) {
    TCanvas *c1 = new TCanvas("c1", "c1");
    clusterHist->Draw();
	}


  clusterHist->Write();
	outFile->Write(); // Write all of objects to the output file.
	outFile->Close(); // Close the output file.

  return 1;

}


void load_library(){

#if !defined(__MAKECINT__)
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.rootmap");
    gSystem->Load("${WCSIMDIR}/src/WCSimRootDict_rdict.pcm");
  }else{
    gSystem->Load("../../libWCSimRoot.so");
    gSystem->Load("../../libWCSimRoot.rootmap");
    gSystem->Load("../../src/WCSimRootDict_rdict.pcm");
  }
#endif
  return;

}

TString create_filename(const char * prefix, TString& filename_string)
{
  //std::cout << "Creating filename from prefix " << prefix << " and filename_string " << filename_string << std::endl;                                                                                                                                                  
  TString prefix_string(prefix);
  TString outfilename = prefix_string + filename_string;
  return outfilename;
}

TFile * get_input_file(char *filename){

  TFile *f;
  // Open the file                                                                                                            
  f = new TFile(filename,"read");

  if (!f->IsOpen()){
    std::cerr << "Error, could not open input file: " << filename << std::endl;
    exit(0);
  }
  
  return f;

}



