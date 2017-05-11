
#ifndef LIBRARY_DAQ_H
#define LIBRARY_DAQ_H

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <thrust/extrema.h>
#include <limits>
#include <limits.h>

/////////////////////////////
// define global variables //
/////////////////////////////
/// parameters
double distance_between_vertices; // linear distance between test vertices
double wall_like_distance; // distance from wall (in units of distance_between_vertices) to define wall-like events
unsigned int time_step_size; // time binning for the trigger
__constant__ unsigned int constant_time_step_size; 
unsigned int water_like_threshold_number_of_pmts; // number of pmts above which a trigger is possible for water-like events
unsigned int wall_like_threshold_number_of_pmts; // number of pmts above which a trigger is possible for wall-like events
unsigned int nhits_threshold_min, nhits_threshold_max;
double coalesce_time; // time such that if two triggers are closer than this they are coalesced into a single trigger
double trigger_gate_up; // duration to be saved after the trigger time
double trigger_gate_down; // duration to be saved before the trigger time
unsigned int max_n_hits_per_job; // max n of pmt hits per job
double dark_rate;
/// detector
double detector_height; // detector height
double detector_radius; // detector radius
/// pmts
unsigned int n_PMTs; // number of pmts in the detector
__constant__ unsigned int constant_n_PMTs;
double * PMT_x, *PMT_y, *PMT_z; // coordinates of the pmts in the detector
/// vertices
unsigned int n_test_vertices; // number of test vertices
unsigned int n_water_like_test_vertices; // number of test vertices
__constant__ unsigned int constant_n_test_vertices;
__constant__ unsigned int constant_n_water_like_test_vertices;
double * vertex_x, * vertex_y, * vertex_z; // coordinates of test vertices
/// threads
unsigned int number_of_kernel_blocks;  // number of cores to be used
dim3 number_of_kernel_blocks_3d;
unsigned int number_of_threads_per_block; // number of threads per core to be used
dim3 number_of_threads_per_block_3d;
unsigned int grid_size;  // grid = (n cores) X (n threads / core)
/// hits
double time_offset;  // ns, offset to make times positive
__constant__ double constant_time_offset;
unsigned int n_time_bins; // number of time bins 
__constant__ unsigned int constant_n_time_bins;
unsigned int n_direction_bins_theta; // number of direction bins 
__constant__ unsigned int constant_n_direction_bins_theta;
unsigned int n_direction_bins_phi; // number of direction bins 
__constant__ unsigned int constant_n_direction_bins_phi;
unsigned int n_direction_bins; // number of direction bins 
__constant__ unsigned int constant_n_direction_bins;
unsigned int n_hits; // number of input hits from the detector
__constant__ unsigned int constant_n_hits;
unsigned int * host_ids; // pmt id of a hit
unsigned int *device_ids;
texture<unsigned int, 1, cudaReadModeElementType> tex_ids;
unsigned int * host_times;  // time of a hit
unsigned int *device_times;
texture<unsigned int, 1, cudaReadModeElementType> tex_times;
// corrected tim bin of each hit (for each vertex)
unsigned int * host_time_bin_of_hit;
unsigned int * device_time_bin_of_hit;
// npmts per time bin
unsigned int * device_n_pmts_per_time_bin; // number of active pmts in a time bin
unsigned int * host_n_pmts_per_time_bin;
unsigned int * device_n_pmts_nhits; // number of active pmts
unsigned int * host_n_pmts_nhits;
unsigned int * device_n_pmts_per_time_bin_and_direction_bin; // number of active pmts in a time bin and direction bin
//unsigned int * device_time_nhits; // trigger time
//unsigned int * host_time_nhits;
// tof
double speed_light_water;
double cerenkov_angle_water;
double twopi;
bool cylindrical_grid;
float *device_times_of_flight; // time of flight between a vertex and a pmt
float *host_times_of_flight;
bool *device_directions_for_vertex_and_pmt; // test directions for vertex and pmt
bool *host_directions_for_vertex_and_pmt;
texture<float, 1, cudaReadModeElementType> tex_times_of_flight;
//texture<bool, 1, cudaReadModeElementType> tex_directions_for_vertex_and_pmt;
// triggers
std::vector<std::pair<unsigned int,unsigned int> > candidate_trigger_pair_vertex_time;  // pair = (v, t) = (a vertex, a time at the end of the 2nd of two coalesced bins)
std::vector<unsigned int> candidate_trigger_npmts_in_time_bin; // npmts in time bin
std::vector<std::pair<unsigned int,unsigned int> > trigger_pair_vertex_time;
std::vector<unsigned int> trigger_npmts_in_time_bin;
std::vector<std::pair<unsigned int,unsigned int> > final_trigger_pair_vertex_time;
std::vector<double> output_trigger_information;
// C timing
struct timeval t0;
struct timeval t1;
// CUDA timing
cudaEvent_t start, stop, total_start, total_stop;
// make output txt file for plotting?
bool output_txt;
unsigned int correct_mode;
// find candidates
unsigned int * host_max_number_of_pmts_in_time_bin;
unsigned int * device_max_number_of_pmts_in_time_bin;
unsigned int *  host_vertex_with_max_n_pmts;
unsigned int *  device_vertex_with_max_n_pmts;
// gpu properties
int max_n_threads_per_block;
int max_n_blocks;
// verbosity level
bool use_verbose;
bool use_timing;
// files
std::string event_file;
std::string detector_file;
std::string pmts_file;
std::string output_file;
std::string event_file_base;
std::string event_file_suffix;
std::string output_file_base;
float elapsed_parameters, elapsed_pmts, elapsed_detector, elapsed_vertices,
  elapsed_threads, elapsed_tof, elapsed_directions, elapsed_memory_tofs_dev, elapsed_memory_directions_dev, elapsed_memory_candidates_host, elapsed_tofs_copy_dev,  elapsed_directions_copy_dev,
  elapsed_input, elapsed_memory_dev, elapsed_copy_dev, elapsed_kernel_correct_times_and_get_n_pmts_per_time_bin, 
  elapsed_threads_candidates, elapsed_candidates_memory_dev, elapsed_candidates_kernel,
  elapsed_candidates_copy_host, choose_candidates, elapsed_coalesce, elapsed_gates, elapsed_free, elapsed_total,
  elapsed_tofs_free, elapsed_reset, elapsed_write_output;
unsigned int greatest_divisor;
unsigned int the_max_time;
unsigned int nhits_window;
int n_events;


__global__ void kernel_find_vertex_with_max_npmts_in_timebin(unsigned int * np, unsigned int * mnp, unsigned int * vmnp);
__global__ void kernel_find_vertex_with_max_npmts_in_timebin_and_directionbin(unsigned int * np, unsigned int * mnp, unsigned int * vmnp);

unsigned int read_number_of_input_hits();
bool read_input();
void print_parameters();
void print_parameters_2d();
void print_input();
void print_times_of_flight();
void print_directions();
void print_pmts();
bool read_detector();
bool read_the_pmts();
bool read_the_detector();
void make_test_vertices();
bool setup_threads_for_tof();
bool setup_threads_for_tof_biparallel();
bool setup_threads_for_tof_2d(unsigned int A, unsigned int B);
bool setup_threads_to_find_candidates();
bool setup_threads_nhits();
bool read_the_input();
void allocate_tofs_memory_on_device();
void allocate_directions_memory_on_device();
void allocate_correct_memory_on_device();
void allocate_correct_memory_on_device_nhits();
void allocate_candidates_memory_on_host();
void allocate_candidates_memory_on_device();
void make_table_of_tofs();
void make_table_of_directions();
void fill_correct_memory_on_device();
void fill_tofs_memory_on_device();
void fill_directions_memory_on_device();
void fill_tofs_memory_on_device_nhits();
void coalesce_triggers();
void separate_triggers_into_gates();
float timedifference_msec(struct timeval t0, struct timeval t1);
void start_c_clock();
double stop_c_clock();
void start_cuda_clock();
double stop_cuda_clock();
void start_total_cuda_clock();
double stop_total_cuda_clock();
unsigned int get_distance_index(unsigned int pmt_id, unsigned int vertex_block);
unsigned int get_time_index(unsigned int hit_index, unsigned int vertex_block);
unsigned int get_direction_index_at_pmt(unsigned int pmt_id, unsigned int vertex_index, unsigned int direction_index);
unsigned int get_direction_index_at_angles(unsigned int iphi, unsigned int itheta);
unsigned int get_direction_index_at_time(unsigned int time_bin, unsigned int vertex_index, unsigned int direction_index);
__device__ unsigned int device_get_distance_index(unsigned int pmt_id, unsigned int vertex_block);
__device__ unsigned int device_get_time_index(unsigned int hit_index, unsigned int vertex_block);
__device__ unsigned int device_get_direction_index_at_pmt(unsigned int pmt_id, unsigned int vertex_index, unsigned int direction_index);
__device__ unsigned int device_get_direction_index_at_angles(unsigned int iphi, unsigned int itheta);
__device__ unsigned int device_get_direction_index_at_time(unsigned int time_bin, unsigned int vertex_index, unsigned int direction_index);
void print_gpu_properties();
unsigned int read_number_of_pmts();
bool read_pmts();
void free_event_memories();
void free_event_memories_nhits();
void free_global_memories();
void copy_candidates_from_device_to_host();
void choose_candidates_above_threshold();
bool set_input_file_for_event(int n);
void set_output_file();
void set_output_file_nhits(unsigned int threshold);
void write_output();
void write_output_nhits(unsigned int n);
void initialize_output();
void initialize_output_nhits();
float read_value_from_file(std::string paramname, std::string filename);
void read_user_parameters();
void read_user_parameters_nhits();
void check_cudamalloc_float(unsigned int size);
void check_cudamalloc_int(unsigned int size);
void check_cudamalloc_unsigned_int(unsigned int size);
void check_cudamalloc_bool(unsigned int size);
void setup_threads_for_histo(unsigned int n);
unsigned int find_greatest_divisor(unsigned int n, unsigned int max);
void setup_threads_for_histo();
void setup_threads_for_histo_iterated(bool last);
void setup_threads_for_histo_per(unsigned int n);


unsigned int read_number_of_input_hits(){

  FILE *f=fopen(event_file.c_str(), "r");
  if (f == NULL){
    printf(" cannot read input file \n");
    fclose(f);
    return 0;
  }

  unsigned int n_hits = 0;

  for (char c = getc(f); c != EOF; c = getc(f))
    if (c == '\n')
      n_hits ++;

  fclose(f);
  return n_hits;

}

bool read_input(){

  FILE *f=fopen(event_file.c_str(), "r");

  int time;
  double timef;
  unsigned int id;
  int min = INT_MAX;
  int max = INT_MIN;
  for( unsigned int i=0; i<n_hits; i++){
    if( fscanf(f, "%d %lf", &id, &timef) != 2 ){
      printf(" problem scanning hit %d from input \n", i);
      fclose(f);
      return false;
    }
    time = int(floor(timef));
    host_times[i] = time;
    host_ids[i] = id;
    if( time > max ) max = time;
    if( time < min ) min = time;
  }

  if( min < 0 ){
    time_offset -= min;
  }


  the_max_time = max;

  fclose(f);
  return true;

}


bool read_detector(){

  FILE *f=fopen(detector_file.c_str(), "r");
  double pmt_radius;
  if( fscanf(f, "%lf %lf %lf", &detector_height, &detector_radius, &pmt_radius) != 3 ){
    printf(" problem scanning detector \n");
    fclose(f);
    return false;
  }

  fclose(f);
  return true;

}



void print_parameters(){

  printf(" n_test_vertices = %d \n", n_test_vertices);
  printf(" n_water_like_test_vertices = %d \n", n_water_like_test_vertices);
  printf(" n_PMTs = %d \n", n_PMTs);
  printf(" number_of_kernel_blocks = %d \n", number_of_kernel_blocks);
  printf(" number_of_threads_per_block = %d \n", number_of_threads_per_block);
  printf(" grid size = %d -> %d \n", number_of_kernel_blocks*number_of_threads_per_block, grid_size);

}

void print_parameters_2d(){

  printf(" n_test_vertices = %d \n", n_test_vertices);
  printf(" n_water_like_test_vertices = %d \n", n_water_like_test_vertices);
  printf(" n_PMTs = %d \n", n_PMTs);
  printf(" number_of_kernel_blocks = (%d, %d) = %d \n", number_of_kernel_blocks_3d.x, number_of_kernel_blocks_3d.y, number_of_kernel_blocks_3d.x * number_of_kernel_blocks_3d.y);
  printf(" number_of_threads_per_block = (%d, %d) = %d \n", number_of_threads_per_block_3d.x, number_of_threads_per_block_3d.y, number_of_threads_per_block_3d.x * number_of_threads_per_block_3d.y);
  printf(" grid size = %d -> %d \n", number_of_kernel_blocks_3d.x*number_of_kernel_blocks_3d.y*number_of_threads_per_block_3d.x*number_of_threads_per_block_3d.y, grid_size);

}

void print_input(){

  for(unsigned int i=0; i<n_hits; i++)
    printf(" hit %d time %d id %d \n", i, host_times[i], host_ids[i]);

}

void print_pmts(){

  for(unsigned int i=0; i<n_PMTs; i++)
    printf(" pmt %d x %f y %f z %f  \n", i, PMT_x[i], PMT_y[i], PMT_z[i]);

}

void print_times_of_flight(){

  printf(" times_of_flight: (vertex, PMT) \n");
  unsigned int distance_index;
  for(unsigned int iv=0; iv<n_test_vertices; iv++){
    printf(" ( ");
    for(unsigned int ip=0; ip<n_PMTs; ip++){
      distance_index = get_distance_index(ip + 1, n_PMTs*iv);
      printf(" %f ", host_times_of_flight[distance_index]);
    }
    printf(" ) \n");
  }
}


void print_directions(){

  printf(" directions: (vertex, PMT) \n");
  for(unsigned int iv=0; iv<n_test_vertices; iv++){
    printf(" [ ");
    for(unsigned int ip=0; ip<n_PMTs; ip++){
      printf(" ( ");
      for(unsigned int id=0; id<n_direction_bins; id++){
	printf("%d ", host_directions_for_vertex_and_pmt[get_direction_index_at_pmt(ip, iv, id)]);
      }
      printf(" ) ");
    }
    printf(" ] \n");
  }
}


bool read_the_pmts(){

  printf(" --- read pmts \n");
  n_PMTs = read_number_of_pmts();
  if( !n_PMTs ) return false;
  printf(" detector contains %d PMTs \n", n_PMTs);
  PMT_x = (double *)malloc(n_PMTs*sizeof(double));
  PMT_y = (double *)malloc(n_PMTs*sizeof(double));
  PMT_z = (double *)malloc(n_PMTs*sizeof(double));
  if( !read_pmts() ) return false;
  //print_pmts();
  return true;

}

bool read_the_detector(){

  printf(" --- read detector \n");
  if( !read_detector() ) return false;
  printf(" detector height %f cm, radius %f cm \n", detector_height, detector_radius);
  return true;

}

void make_test_vertices(){

  printf(" --- make test vertices \n");
  float semiheight = detector_height/2.;
  n_test_vertices = 0;


  if( !cylindrical_grid ){

    // 1: count number of test vertices
    for(int i=-1*semiheight; i <= semiheight; i+=distance_between_vertices) {
      for(int j=-1*detector_radius; j<=detector_radius; j+=distance_between_vertices) {
	for(int k=-1*detector_radius; k<=detector_radius; k+=distance_between_vertices) {
	  if(pow(j,2)+pow(k,2) > pow(detector_radius,2))
	    continue;
	  n_test_vertices++;
	}
      }
    }
    vertex_x = (double *)malloc(n_test_vertices*sizeof(double));
    vertex_y = (double *)malloc(n_test_vertices*sizeof(double));
    vertex_z = (double *)malloc(n_test_vertices*sizeof(double));

    // 2: assign coordinates to test vertices
    // water-like events
    n_test_vertices = 0;
    for(int i=-1*semiheight; i <= semiheight; i+=distance_between_vertices) {
      for(int j=-1*detector_radius; j<=detector_radius; j+=distance_between_vertices) {
	for(int k=-1*detector_radius; k<=detector_radius; k+=distance_between_vertices) {

	
	  if( 
	     // skip endcap region
	     abs(i) > semiheight - wall_like_distance*distance_between_vertices ||
	     // skip sidewall region
	     pow(j,2)+pow(k,2) > pow(detector_radius - wall_like_distance*distance_between_vertices,2)
	      ) continue;
	
	  vertex_x[n_test_vertices] = j*1.;
	  vertex_y[n_test_vertices] = k*1.;
	  vertex_z[n_test_vertices] = i*1.;
	  n_test_vertices++;
	}
      }
    }
    n_water_like_test_vertices = n_test_vertices;

    // wall-like events
    for(int i=-1*semiheight; i <= semiheight; i+=distance_between_vertices) {
      for(int j=-1*detector_radius; j<=detector_radius; j+=distance_between_vertices) {
	for(int k=-1*detector_radius; k<=detector_radius; k+=distance_between_vertices) {

	  if( 
	     abs(i) > semiheight - wall_like_distance*distance_between_vertices ||
	     pow(j,2)+pow(k,2) > pow(detector_radius - wall_like_distance*distance_between_vertices,2)
	      ){

	    if(pow(j,2)+pow(k,2) > pow(detector_radius,2)) continue;
	  
	    vertex_x[n_test_vertices] = j*1.;
	    vertex_y[n_test_vertices] = k*1.;
	    vertex_z[n_test_vertices] = i*1.;
	    n_test_vertices++;
	  }
	}
      }
    }

  }else{ // cylindrical grid
  
    int n_vertical = detector_height/distance_between_vertices;
    double distance_vertical = detector_height/n_vertical;
    int n_radial = 2.*detector_radius/distance_between_vertices;
    double distance_radial = 2.*detector_radius/n_radial;
    int n_angular;
    double distance_angular;
    
    printf(" distance_between_vertices %f, distance_vertical %f, distance_radial %f \n",
	   distance_between_vertices, distance_vertical, distance_radial);
    
    double the_r, the_z, the_phi;
    bool first = false; // true: add extra layer near wall
                       // false: regular spacing

    bool add_extra_layer = first;
    
    // 1: count number of test vertices
    the_r = detector_radius;
    while( the_r >= 0. ){
      n_angular = twopi*the_r / distance_between_vertices;
      distance_angular = twopi/n_angular;
      
      the_z = -semiheight;
      
      while( the_z <= semiheight){
	
	the_phi = 0.;
	while( the_phi < twopi - distance_angular/2. ){
	  
	  n_test_vertices ++;
	  
	  if( the_r == 0. ) break;
	  the_phi += distance_angular;
	}


	if( add_extra_layer ){
	  if( the_z + semiheight < 0.3*distance_vertical ) // only true at bottom endcap
	    the_z += distance_vertical/2.;
	  else if( semiheight - the_z < 0.7*distance_vertical ) // only true near top endcap
	    the_z += distance_vertical/2.;
	  else
	    the_z += distance_vertical;
	}else{
	  the_z += distance_vertical;
	}
      }
      if( first ){
	the_r -= distance_radial/2.;
	first = false;
      }
      else{
	the_r -= distance_radial;
      }
    }

    vertex_x = (double *)malloc(n_test_vertices*sizeof(double));
    vertex_y = (double *)malloc(n_test_vertices*sizeof(double));
    vertex_z = (double *)malloc(n_test_vertices*sizeof(double));

    first = add_extra_layer;
    // 2: assign coordinates to test vertices
    // water-like events
    n_test_vertices = 0;

    the_r = detector_radius;
    while( the_r >= 0. ){

      // skip sidewall region
      if(the_r <= detector_radius - wall_like_distance*distance_between_vertices ){

	n_angular = twopi*the_r / distance_between_vertices;
	distance_angular = twopi/n_angular;
	
	the_z = -semiheight;
	
	while( the_z <= semiheight){
	  
	  // skip endcap region
	  if( fabs(the_z) <= semiheight - wall_like_distance*distance_between_vertices ){

	    the_phi = 0.;
	    while( the_phi < twopi - distance_angular/2. ){
	      
	      vertex_x[n_test_vertices] = the_r*cos(the_phi);
	      vertex_y[n_test_vertices] = the_r*sin(the_phi);
	      vertex_z[n_test_vertices] = the_z;
	      n_test_vertices ++;
	      
	      if( the_r == 0. ) break;
	      the_phi += distance_angular;
	    }
	  }

	  if( add_extra_layer ){
	    if( the_z + semiheight < 0.3*distance_vertical ) // only true at bottom endcap
	      the_z += distance_vertical/2.;
	    else if( semiheight - the_z < 0.7*distance_vertical ) // only true near top endcap
	      the_z += distance_vertical/2.;
	    else
	      the_z += distance_vertical;
	  }else{
	    the_z += distance_vertical;
	  }
	}

      }
      if( first ){
	the_r -= distance_radial/2.;
	first = false;
      }
      else{
	the_r -= distance_radial;
      }
    }


    n_water_like_test_vertices = n_test_vertices;

    first = add_extra_layer;
    // wall-like events
    the_r = detector_radius;
    while( the_r >= 0. ){

      n_angular = twopi*the_r / distance_between_vertices;
      distance_angular = twopi/n_angular;
      
      the_z = -semiheight;
      
      while( the_z <= semiheight){
	
	if( fabs(the_z) > semiheight - wall_like_distance*distance_between_vertices ||
	    the_r > detector_radius - wall_like_distance*distance_between_vertices ){
	  
	  the_phi = 0.;
	  while( the_phi < twopi - distance_angular/2. ){
	    
	    vertex_x[n_test_vertices] = the_r*cos(the_phi);
	    vertex_y[n_test_vertices] = the_r*sin(the_phi);
	    vertex_z[n_test_vertices] = the_z;
	    n_test_vertices ++;
	    
	    if( the_r == 0. ) break;
	    the_phi += distance_angular;
	  }
	}
	if( add_extra_layer ){
	  if( the_z + semiheight < 0.3*distance_vertical ) // only true at bottom endcap
	    the_z += distance_vertical/2.;
	  else if( semiheight - the_z < 0.7*distance_vertical ) // only true near top endcap
	    the_z += distance_vertical/2.;
	  else
	    the_z += distance_vertical;
	}else{
	  the_z += distance_vertical;
	}
      }
      if( first ){
	the_r -= distance_radial/2.;
	first = false;
      }
      else{
	the_r -= distance_radial;
      }
    }
    

  }

  return;

}

bool setup_threads_for_tof(){

  grid_size = n_test_vertices;

  number_of_kernel_blocks = grid_size / max_n_threads_per_block + 1;
  number_of_threads_per_block = ( number_of_kernel_blocks > 1 ? max_n_threads_per_block : grid_size);

  print_parameters();

  if( number_of_threads_per_block > max_n_threads_per_block ){
    printf(" warning: number_of_threads_per_block = %d cannot exceed max value %d \n", number_of_threads_per_block, max_n_threads_per_block );
    return false;
  }

  if( number_of_kernel_blocks > max_n_blocks ){
    printf(" warning: number_of_kernel_blocks = %d cannot exceed max value %d \n", number_of_kernel_blocks, max_n_blocks );
    return false;
  }

  return true;
}


bool setup_threads_for_tof_biparallel(){

  grid_size = n_test_vertices * n_hits;

  number_of_kernel_blocks = grid_size / max_n_threads_per_block + 1;
  number_of_threads_per_block = ( number_of_kernel_blocks > 1 ? max_n_threads_per_block : grid_size);

  print_parameters();

  if( number_of_threads_per_block > max_n_threads_per_block ){
    printf(" --------------------- warning: number_of_threads_per_block = %d cannot exceed max value %d \n", number_of_threads_per_block, max_n_threads_per_block );
    return false;
  }

  if( number_of_kernel_blocks > max_n_blocks ){
    printf(" warning: number_of_kernel_blocks = %d cannot exceed max value %d \n", number_of_kernel_blocks, max_n_blocks );
    return false;
  }

  return true;

}

bool setup_threads_for_tof_2d(unsigned int A, unsigned int B){

  if( std::numeric_limits<unsigned int>::max() / B  < A ){
    printf(" --------------------- warning: B = %d times A = %d cannot exceed max value %u \n", B, A, std::numeric_limits<unsigned int>::max() );
    return false;
  }

  grid_size = A * B;
  unsigned int max_n_threads_per_block_2d = sqrt(max_n_threads_per_block);

  number_of_kernel_blocks_3d.x = A / max_n_threads_per_block_2d + 1;
  number_of_kernel_blocks_3d.y = B / max_n_threads_per_block_2d + 1;

  number_of_threads_per_block_3d.x = ( number_of_kernel_blocks_3d.x > 1 ? max_n_threads_per_block_2d : A);
  number_of_threads_per_block_3d.y = ( number_of_kernel_blocks_3d.y > 1 ? max_n_threads_per_block_2d : B);

  print_parameters_2d();

  if( number_of_threads_per_block_3d.x > max_n_threads_per_block_2d ){
    printf(" --------------------- warning: number_of_threads_per_block x = %d cannot exceed max value %d \n", number_of_threads_per_block_3d.x, max_n_threads_per_block_2d );
    return false;
  }

  if( number_of_threads_per_block_3d.y > max_n_threads_per_block_2d ){
    printf(" --------------------- warning: number_of_threads_per_block y = %d cannot exceed max value %d \n", number_of_threads_per_block_3d.y, max_n_threads_per_block_2d );
    return false;
  }

  if( number_of_kernel_blocks_3d.x > max_n_blocks ){
    printf(" warning: number_of_kernel_blocks x = %d cannot exceed max value %d \n", number_of_kernel_blocks_3d.x, max_n_blocks );
    return false;
  }

  if( number_of_kernel_blocks_3d.y > max_n_blocks ){
    printf(" warning: number_of_kernel_blocks y = %d cannot exceed max value %d \n", number_of_kernel_blocks_3d.y, max_n_blocks );
    return false;
  }

  if( std::numeric_limits<int>::max() / (number_of_kernel_blocks_3d.x*number_of_kernel_blocks_3d.y)  < number_of_threads_per_block_3d.x*number_of_threads_per_block_3d.y ){
    printf(" --------------------- warning: grid size cannot exceed max value %u \n", std::numeric_limits<int>::max() );
    return false;
  }


  return true;

}

bool setup_threads_to_find_candidates(){

  number_of_kernel_blocks = n_time_bins / max_n_threads_per_block + 1;
  number_of_threads_per_block = ( number_of_kernel_blocks > 1 ? max_n_threads_per_block : n_time_bins);

  if( number_of_threads_per_block > max_n_threads_per_block ){
    printf(" warning: number_of_threads_per_block = %d cannot exceed max value %d \n", number_of_threads_per_block, max_n_threads_per_block );
    return false;
  }

  return true;
}

bool setup_threads_nhits(){

  number_of_kernel_blocks_3d.x = 100;
  number_of_kernel_blocks_3d.y = 1;

  number_of_threads_per_block_3d.x = 1024;
  number_of_threads_per_block_3d.y = 1;

  print_parameters_2d();

  return true;

}



bool read_the_input(){

  printf(" --- read input \n");
  n_hits = read_number_of_input_hits();
  if( !n_hits ) return false;
  printf(" input contains %d hits \n", n_hits);
  host_ids = (unsigned int *)malloc(n_hits*sizeof(unsigned int));
  host_times = (unsigned int *)malloc(n_hits*sizeof(unsigned int));
  if( !read_input() ) return false;
  //time_offset = 600.; // set to constant to match trevor running
  n_time_bins = int(floor((the_max_time + time_offset)/time_step_size))+1; // floor returns the integer below
  printf(" input max_time %d, n_time_bins %d \n", the_max_time, n_time_bins);
  printf(" time_offset = %f ns \n", time_offset);
  //print_input();

  checkCudaErrors( cudaMemcpyToSymbol(constant_n_time_bins, &n_time_bins, sizeof(n_time_bins)) );
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_hits, &n_hits, sizeof(n_hits)) );

  return true;
}

void allocate_tofs_memory_on_device(){

  printf(" --- allocate memory tofs \n");
  check_cudamalloc_float(n_test_vertices*n_PMTs);
  checkCudaErrors(cudaMalloc((void **)&device_times_of_flight, n_test_vertices*n_PMTs*sizeof(float)));
  /*
  if( n_hits*n_test_vertices > available_memory ){
    printf(" cannot allocate vector of %d, available_memory %d \n", n_hits*n_test_vertices, available_memory);
    return 0;
  }
  */


  if( correct_mode == 1 ){
    
    unsigned int max = max_n_threads_per_block;
    greatest_divisor = find_greatest_divisor ( n_test_vertices , max);
    printf(" greatest divisor of %d below %d is %d \n", n_test_vertices, max, greatest_divisor);
    
  }

  return;

}

void allocate_directions_memory_on_device(){

  printf(" --- allocate memory directions \n");
  check_cudamalloc_bool( n_test_vertices*n_direction_bins*n_PMTs);
  checkCudaErrors(cudaMalloc((void **)&device_directions_for_vertex_and_pmt, n_test_vertices*n_direction_bins*n_PMTs*sizeof(bool)));
  /*
  if( n_hits*n_test_vertices > available_memory ){
    printf(" cannot allocate vector of %d, available_memory %d \n", n_hits*n_test_vertices, available_memory);
    return 0;
  }
  */

  return;

}

void allocate_correct_memory_on_device(){

  printf(" --- allocate memory \n");
  /*
  if( n_hits > available_memory ){
    printf(" cannot allocate vector of %d, available_memory %d \n", n_hits, available_memory);
    return 0;
  }
  */
  check_cudamalloc_unsigned_int(n_hits);
  checkCudaErrors(cudaMalloc((void **)&device_ids, n_hits*sizeof(unsigned int)));

  check_cudamalloc_unsigned_int(n_hits);
  checkCudaErrors(cudaMalloc((void **)&device_times, n_hits*sizeof(unsigned int)));
  /*
  if( n_test_vertices*n_PMTs > available_memory ){
    printf(" cannot allocate vector of %d, available_memory %d \n", n_test_vertices*n_PMTs, available_memory);
    return 0;
  }
  */

  if( correct_mode != 9 ){
    check_cudamalloc_unsigned_int(n_time_bins*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_n_pmts_per_time_bin, n_time_bins*n_test_vertices*sizeof(unsigned int)));
  }

  if( correct_mode == 0 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));
  } else if( correct_mode == 1 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));

    check_cudamalloc_unsigned_int(n_hits*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_time_bin_of_hit, n_hits*n_test_vertices*sizeof(unsigned int)));
    //    checkCudaErrors(cudaMemset(device_time_bin_of_hit, 0, n_hits*n_test_vertices*sizeof(unsigned int)));
  }
  else if( correct_mode == 2 ){
    check_cudamalloc_unsigned_int(n_hits*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_time_bin_of_hit, n_hits*n_test_vertices*sizeof(unsigned int)));
    //checkCudaErrors(cudaMemset(device_time_bin_of_hit, 0, n_hits*n_test_vertices*sizeof(unsigned int)));

    host_time_bin_of_hit = (unsigned int *)malloc(n_hits*n_test_vertices*sizeof(unsigned int));

    host_n_pmts_per_time_bin = (unsigned int *)malloc(n_time_bins*n_test_vertices*sizeof(unsigned int));
    memset(host_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int));
  } else if( correct_mode == 3 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));

    check_cudamalloc_unsigned_int(n_hits*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_time_bin_of_hit, n_hits*n_test_vertices*sizeof(unsigned int)));
    //checkCudaErrors(cudaMemset(device_time_bin_of_hit, 0, n_hits*n_test_vertices*sizeof(unsigned int)));
  } else if( correct_mode == 4 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));

    check_cudamalloc_unsigned_int(n_hits*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_time_bin_of_hit, n_hits*n_test_vertices*sizeof(unsigned int)));
  } else if( correct_mode == 5 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));

    check_cudamalloc_unsigned_int(n_hits*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_time_bin_of_hit, n_hits*n_test_vertices*sizeof(unsigned int)));
  } else if( correct_mode == 6 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));

    check_cudamalloc_unsigned_int(n_hits*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_time_bin_of_hit, n_hits*n_test_vertices*sizeof(unsigned int)));
  } else if( correct_mode == 7 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));

    check_cudamalloc_unsigned_int(n_hits*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_time_bin_of_hit, n_hits*n_test_vertices*sizeof(unsigned int)));
  } else if( correct_mode == 8 ){
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin, 0, n_time_bins*n_test_vertices*sizeof(unsigned int)));
  } else if( correct_mode == 9 ){

    check_cudamalloc_unsigned_int(n_time_bins*n_direction_bins*n_test_vertices);
    checkCudaErrors(cudaMalloc((void **)&device_n_pmts_per_time_bin_and_direction_bin, n_time_bins*n_direction_bins*n_test_vertices*sizeof(unsigned int)));
    checkCudaErrors(cudaMemset(device_n_pmts_per_time_bin_and_direction_bin, 0, n_time_bins*n_direction_bins*n_test_vertices*sizeof(unsigned int)));
  }


  return;

}


void allocate_correct_memory_on_device_nhits(){

  printf(" --- allocate memory \n");
  /*
  if( n_hits > available_memory ){
    printf(" cannot allocate vector of %d, available_memory %d \n", n_hits, available_memory);
    return 0;
  }
  */
  check_cudamalloc_unsigned_int(n_hits);
  checkCudaErrors(cudaMalloc((void **)&device_ids, n_hits*sizeof(unsigned int)));

  check_cudamalloc_unsigned_int(n_hits);
  checkCudaErrors(cudaMalloc((void **)&device_times, n_hits*sizeof(unsigned int)));
  /*
  if( n_test_vertices*n_PMTs > available_memory ){
    printf(" cannot allocate vector of %d, available_memory %d \n", n_test_vertices*n_PMTs, available_memory);
    return 0;
  }
  */

  check_cudamalloc_unsigned_int(1);
  checkCudaErrors(cudaMalloc((void **)&device_n_pmts_nhits, 1*sizeof(unsigned int)));
  //  checkCudaErrors(cudaMalloc((void **)&device_time_nhits, (nhits_window/time_step_size + 1)*sizeof(unsigned int)));

  host_n_pmts_nhits = (unsigned int *)malloc(1*sizeof(unsigned int));
  //  host_time_nhits = (unsigned int *)malloc((nhits_window/time_step_size + 1)*sizeof(unsigned int));

  return;

}

void allocate_candidates_memory_on_host(){

  printf(" --- allocate candidates memory on host \n");

  host_max_number_of_pmts_in_time_bin = (unsigned int *)malloc(n_time_bins*sizeof(unsigned int));
  host_vertex_with_max_n_pmts = (unsigned int *)malloc(n_time_bins*sizeof(unsigned int));

  return;

}

void allocate_candidates_memory_on_device(){

  printf(" --- allocate candidates memory on device \n");

  check_cudamalloc_unsigned_int(n_time_bins);
  checkCudaErrors(cudaMalloc((void **)&device_max_number_of_pmts_in_time_bin, n_time_bins*sizeof(unsigned int)));

  check_cudamalloc_unsigned_int(n_time_bins);
  checkCudaErrors(cudaMalloc((void **)&device_vertex_with_max_n_pmts, n_time_bins*sizeof(unsigned int)));

  return;

}

void make_table_of_tofs(){

  printf(" --- fill times_of_flight \n");
  host_times_of_flight = (float*)malloc(n_test_vertices*n_PMTs * sizeof(double));
  printf(" speed_light_water %f \n", speed_light_water);
  unsigned int distance_index;
  time_offset = 0.;
  for(unsigned int ip=0; ip<n_PMTs; ip++){
    for(unsigned int iv=0; iv<n_test_vertices; iv++){
      distance_index = get_distance_index(ip + 1, n_PMTs*iv);
      host_times_of_flight[distance_index] = sqrt(pow(vertex_x[iv] - PMT_x[ip],2) + pow(vertex_y[iv] - PMT_y[ip],2) + pow(vertex_z[iv] - PMT_z[ip],2))/speed_light_water;
      if( host_times_of_flight[distance_index] > time_offset )
	time_offset = host_times_of_flight[distance_index];

    }
  }
  //print_times_of_flight();

  return;
}


void make_table_of_directions(){

  printf(" --- fill directions \n");
  printf(" cerenkov_angle_water %f \n", cerenkov_angle_water);
  host_directions_for_vertex_and_pmt = (bool*)malloc(n_test_vertices*n_PMTs*n_direction_bins * sizeof(bool));
  float dx, dy, dz, dr, phi, cos_theta, sin_theta;
  float phi2, cos_theta2, angle;
  unsigned int dir_index_at_angles;
  unsigned int dir_index_at_pmt;
  for(unsigned int ip=0; ip<n_PMTs; ip++){
    for(unsigned int iv=0; iv<n_test_vertices; iv++){
      dx = PMT_x[ip] - vertex_x[iv];
      dy = PMT_y[ip] - vertex_y[iv];
      dz = PMT_z[ip] - vertex_z[iv];
      dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
      phi = atan2(dy,dx);
      // light direction
      cos_theta = dz/dr;
      sin_theta = sqrt(1. - pow(cos_theta,2));
      // particle direction
      for(unsigned int itheta = 0; itheta < n_direction_bins_theta; itheta++){
	cos_theta2 = -1. + 2.*itheta/(n_direction_bins_theta - 1);
	for(unsigned int iphi = 0; iphi < n_direction_bins_phi; iphi++){
	  phi2 = 0. + twopi*iphi/n_direction_bins_phi;

	  if( (itheta == 0 || itheta + 1 == n_direction_bins_theta ) && iphi != 0 ) break;

	  // angle between light direction and particle direction
	  angle = acos( sin_theta*sqrt(1 - pow(cos_theta2,2)) * cos(phi - phi2) + cos_theta*cos_theta2 );

	  dir_index_at_angles = get_direction_index_at_angles(iphi, itheta);
	  dir_index_at_pmt = get_direction_index_at_pmt(ip, iv, dir_index_at_angles);

	  //printf(" phi %f ctheta %f phi' %f ctheta' %f angle %f dir_index_at_angles %d dir_index_at_pmt %d \n", 
	  //	 phi, cos_theta, phi2, cos_theta2, angle, dir_index_at_angles, dir_index_at_pmt);

	  host_directions_for_vertex_and_pmt[dir_index_at_pmt] 
	    = (bool)(fabs(angle - cerenkov_angle_water) < twopi/(2.*n_direction_bins_phi));
	}
      }
    }
  }
  //print_directions();

  return;
}


void fill_tofs_memory_on_device(){

  printf(" --- copy tofs from host to device \n");
  checkCudaErrors(cudaMemcpy(device_times_of_flight,
			     host_times_of_flight,
			     n_test_vertices*n_PMTs*sizeof(float),
			     cudaMemcpyHostToDevice));
  checkCudaErrors( cudaMemcpyToSymbol(constant_time_step_size, &time_step_size, sizeof(time_step_size)) );
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_test_vertices, &n_test_vertices, sizeof(n_test_vertices)) );
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_water_like_test_vertices, &n_water_like_test_vertices, sizeof(n_water_like_test_vertices)) );
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_PMTs, &n_PMTs, sizeof(n_PMTs)) );

  // Bind the array to the texture
  checkCudaErrors(cudaBindTexture(0,tex_times_of_flight, device_times_of_flight, n_test_vertices*n_PMTs*sizeof(float)));
  


  return;
}


void fill_directions_memory_on_device(){

  printf(" --- copy directions from host to device \n");
  checkCudaErrors(cudaMemcpy(device_directions_for_vertex_and_pmt,
			     host_directions_for_vertex_and_pmt,
			     n_test_vertices*n_PMTs*n_direction_bins*sizeof(bool),
			     cudaMemcpyHostToDevice));
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_direction_bins_theta, &n_direction_bins_theta, sizeof(n_direction_bins_theta)) );
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_direction_bins_phi, &n_direction_bins_phi, sizeof(n_direction_bins_phi)) );
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_direction_bins, &n_direction_bins, sizeof(n_direction_bins)) );

  // Bind the array to the texture
  //  checkCudaErrors(cudaBindTexture(0,tex_directions_for_vertex_and_pmt, device_directions_for_vertex_and_pmt, n_test_vertices*n_PMTs*n_direction_bins_theta*n_direction_bins_theta*sizeof(bool)));
  


  return;
}


void fill_tofs_memory_on_device_nhits(){

  printf(" --- copy tofs from host to device \n");
  checkCudaErrors( cudaMemcpyToSymbol(constant_time_step_size, &time_step_size, sizeof(time_step_size)) );
  checkCudaErrors( cudaMemcpyToSymbol(constant_n_PMTs, &n_PMTs, sizeof(n_PMTs)) );


  return;
}


void fill_correct_memory_on_device(){

  printf(" --- copy from host to device \n");
  checkCudaErrors(cudaMemcpy(device_ids,
			     host_ids,
			     n_hits*sizeof(unsigned int),
			     cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(device_times,
			     host_times,
			     n_hits*sizeof(unsigned int),
			     cudaMemcpyHostToDevice));
  checkCudaErrors( cudaMemcpyToSymbol(constant_time_offset, &time_offset, sizeof(time_offset)) );

  checkCudaErrors(cudaBindTexture(0,tex_ids, device_ids, n_hits*sizeof(unsigned int)));
  checkCudaErrors(cudaBindTexture(0,tex_times, device_times, n_hits*sizeof(unsigned int)));


  return;
}





unsigned int read_number_of_pmts(){

  FILE *f=fopen(pmts_file.c_str(), "r");
  if (f == NULL){
    printf(" cannot read pmts file \n");
    fclose(f);
    return 0;
  }

  unsigned int n_pmts = 0;

  for (char c = getc(f); c != EOF; c = getc(f))
    if (c == '\n')
      n_pmts ++;

  fclose(f);
  return n_pmts;

}

bool read_pmts(){

  FILE *f=fopen(pmts_file.c_str(), "r");

  double x, y, z;
  unsigned int id;
  for( unsigned int i=0; i<n_PMTs; i++){
    if( fscanf(f, "%d %lf %lf %lf", &id, &x, &y, &z) != 4 ){
      printf(" problem scanning pmt %d \n", i);
      fclose(f);
      return false;
    }
    PMT_x[id-1] = x;
    PMT_y[id-1] = y;
    PMT_z[id-1] = z;
  }

  fclose(f);


  return true;

}


void coalesce_triggers(){

  trigger_pair_vertex_time.clear();
  trigger_npmts_in_time_bin.clear();

  unsigned int vertex_index, time_upper, number_of_pmts_in_time_bin;
  unsigned int max_pmt=0,max_vertex_index=0,max_time=0;
  bool first_trigger, last_trigger, coalesce_triggers;
  unsigned int trigger_index;
  for(std::vector<std::pair<unsigned int,unsigned int> >::const_iterator itrigger=candidate_trigger_pair_vertex_time.begin(); itrigger != candidate_trigger_pair_vertex_time.end(); ++itrigger){

    vertex_index =      itrigger->first;
    time_upper = itrigger->second;
    trigger_index = itrigger - candidate_trigger_pair_vertex_time.begin();
    number_of_pmts_in_time_bin = candidate_trigger_npmts_in_time_bin.at(trigger_index);

    first_trigger = (trigger_index == 0);
    last_trigger = (trigger_index == candidate_trigger_pair_vertex_time.size()-1);
       
    if( first_trigger ){
      max_pmt = number_of_pmts_in_time_bin;
      max_vertex_index = vertex_index;
      max_time = time_upper;
    }
    else{
      coalesce_triggers = (std::abs((int)(max_time - time_upper)) < coalesce_time/time_step_size);

      if( coalesce_triggers ){
	if( number_of_pmts_in_time_bin >= max_pmt) {
	  max_pmt = number_of_pmts_in_time_bin;
	  max_vertex_index = vertex_index;
	  max_time = time_upper;
	}
      }else{
	trigger_pair_vertex_time.push_back(std::make_pair(max_vertex_index,max_time));
	trigger_npmts_in_time_bin.push_back(max_pmt);
	max_pmt = number_of_pmts_in_time_bin;
	max_vertex_index = vertex_index;
	max_time = time_upper;     
      }
    }

    if(last_trigger){
      trigger_pair_vertex_time.push_back(std::make_pair(max_vertex_index,max_time));
      trigger_npmts_in_time_bin.push_back(max_pmt);
    }
     
  }

  for(std::vector<std::pair<unsigned int,unsigned int> >::const_iterator itrigger=trigger_pair_vertex_time.begin(); itrigger != trigger_pair_vertex_time.end(); ++itrigger)
    printf(" coalesced trigger timebin %d vertex index %d \n", itrigger->first, itrigger->second);

  return;

}


void separate_triggers_into_gates(){

  final_trigger_pair_vertex_time.clear();
  unsigned int trigger_index;

  unsigned int time_start=0;
  for(std::vector<std::pair<unsigned int,unsigned int> >::const_iterator itrigger=trigger_pair_vertex_time.begin(); itrigger != trigger_pair_vertex_time.end(); ++itrigger){
    //once a trigger is found, we must jump in the future before searching for the next
    if(itrigger->second > time_start) {
      unsigned int triggertime = itrigger->second*time_step_size - time_offset;
      final_trigger_pair_vertex_time.push_back(std::make_pair(itrigger->first,triggertime));
      time_start = triggertime + trigger_gate_up;
      trigger_index = itrigger - trigger_pair_vertex_time.begin();
      output_trigger_information.clear();
      output_trigger_information.push_back(vertex_x[itrigger->first]);
      output_trigger_information.push_back(vertex_y[itrigger->first]);
      output_trigger_information.push_back(vertex_z[itrigger->first]);
      output_trigger_information.push_back(trigger_npmts_in_time_bin.at(trigger_index));
      output_trigger_information.push_back(triggertime);

      printf(" triggertime: %d, npmts: %d, x: %f, y: %f, z: %f \n", triggertime, trigger_npmts_in_time_bin.at(trigger_index), vertex_x[itrigger->first], vertex_y[itrigger->first], vertex_z[itrigger->first]);

      /* if( output_txt ){ */
      /* 	FILE *of=fopen(output_file.c_str(), "w"); */

      /* 	unsigned int distance_index; */
      /* 	double tof; */
      /* 	double corrected_time; */

      /* 	for(unsigned int i=0; i<n_hits; i++){ */

      /* 	  distance_index = get_distance_index(host_ids[i], n_PMTs*(itrigger->first)); */
      /* 	  tof = host_times_of_flight[distance_index]; */

      /* 	  corrected_time = host_times[i]-tof; */

      /* 	  //fprintf(of, " %d %d %f \n", host_ids[i], host_times[i], corrected_time); */
      /* 	  fprintf(of, " %d %f \n", host_ids[i], corrected_time); */
      /* 	} */

      /* 	fclose(of); */
      /* } */

    }
  }


  return;
}


float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}



void start_c_clock(){
  gettimeofday(&t0,0);

}
double stop_c_clock(){
  gettimeofday(&t1,0);
  return timedifference_msec(t0, t1);
}
void start_cuda_clock(){
  cudaEventRecord(start);

}
double stop_cuda_clock(){
  float milli;
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&milli, start, stop);
  return milli;
}
void start_total_cuda_clock(){
  cudaEventRecord(total_start);

}
double stop_total_cuda_clock(){
  float milli;
  cudaEventRecord(total_stop);
  cudaEventSynchronize(total_stop);
  cudaEventElapsedTime(&milli, total_start, total_stop);
  return milli;
}

unsigned int get_distance_index(unsigned int pmt_id, unsigned int vertex_block){
  // block = (npmts) * (vertex index)

  return pmt_id - 1 + vertex_block;

}

unsigned int get_time_index(unsigned int hit_index, unsigned int vertex_block){
  // block = (n time bins) * (vertex index)

  return hit_index + vertex_block;

}

 unsigned int get_direction_index_at_angles(unsigned int iphi, unsigned int itheta){

   if( itheta == 0 ) return 0;
   if( itheta + 1 == n_direction_bins_theta ) return n_direction_bins - 1;

   return 1 + (itheta - 1) * n_direction_bins_phi + iphi;

}

unsigned int get_direction_index_at_pmt(unsigned int pmt_id, unsigned int vertex_index, unsigned int direction_index){

  //                                                     pmt id 1                        ...        pmt id p
  // [                      vertex 1                              vertex 2 ... vertex m] ... [vertex 1 ... vertex m]
  // [(dir 1 ... dir n) (dir 1 ... dir n) ... (dir 1 ... dir n)] ...

  return n_direction_bins * (pmt_id * n_test_vertices  + vertex_index) + direction_index ;

}

unsigned int get_direction_index_at_time(unsigned int time_bin, unsigned int vertex_index, unsigned int direction_index){

  //                                                     time 1                        ...        time p
  // [                      vertex 1                              vertex 2 ... vertex m] ... [vertex 1 ... vertex m]
  // [(dir 1 ... dir n) (dir 1 ... dir n) ... (dir 1 ... dir n)] ...

  return n_direction_bins* (time_bin * n_test_vertices  + vertex_index ) + direction_index ;

}


__device__ unsigned int device_get_distance_index(unsigned int pmt_id, unsigned int vertex_block){
  // block = (npmts) * (vertex index)

  return pmt_id - 1 + vertex_block;

}

__device__ unsigned int device_get_time_index(unsigned int hit_index, unsigned int vertex_block){
  // block = (n time bins) * (vertex index)

  return hit_index + vertex_block;

}

__device__ unsigned int device_get_direction_index_at_pmt(unsigned int pmt_id, unsigned int vertex_index, unsigned int direction_index){

  //                                                     pmt id 1                        ...        pmt id p
  // [                      vertex 1                              vertex 2 ... vertex m] ... [vertex 1 ... vertex m]
  // [(dir 1 ... dir n) (dir 1 ... dir n) ... (dir 1 ... dir n)] ...

  return constant_n_direction_bins * (pmt_id * constant_n_test_vertices  + vertex_index) + direction_index ;

}

__device__ unsigned int device_get_direction_index_at_angles(unsigned int iphi, unsigned int itheta){

   if( itheta == 0 ) return 0;
   if( itheta + 1 == constant_n_direction_bins_theta ) return constant_n_direction_bins - 1;

   return 1 + (itheta - 1) * constant_n_direction_bins_phi + iphi;

}

__device__ unsigned int device_get_direction_index_at_time(unsigned int time_bin, unsigned int vertex_index, unsigned int direction_index){

  //                                                     time 1                        ...        time p
  // [                      vertex 1                              vertex 2 ... vertex m] ... [vertex 1 ... vertex m]
  // [(dir 1 ... dir n) (dir 1 ... dir n) ... (dir 1 ... dir n)] ...

  return constant_n_direction_bins* (time_bin * constant_n_test_vertices  + vertex_index ) + direction_index ;

}

// Print device properties
void print_gpu_properties(){

  int devCount;
  cudaGetDeviceCount(&devCount);
  printf(" CUDA Device Query...\n");
  printf(" There are %d CUDA devices.\n", devCount);
  cudaDeviceProp devProp;
  for (int i = 0; i < devCount; ++i){
    // Get device properties
    printf(" CUDA Device #%d\n", i);
    cudaGetDeviceProperties(&devProp, i);
    printf("Major revision number:          %d\n",  devProp.major);
    printf("Minor revision number:          %d\n",  devProp.minor);
    printf("Name:                           %s\n",  devProp.name);
    printf("Total global memory:            %lu\n",  devProp.totalGlobalMem);
    printf("Total shared memory per block:  %lu\n",  devProp.sharedMemPerBlock);
    printf("Total registers per block:      %d\n",  devProp.regsPerBlock);
    printf("Warp size:                      %d\n",  devProp.warpSize);
    printf("Maximum memory pitch:           %lu\n",  devProp.memPitch);
    max_n_threads_per_block = devProp.maxThreadsPerBlock;
    printf("Maximum threads per block:      %d\n",  max_n_threads_per_block);
    for (int i = 0; i < 3; ++i)
      printf("Maximum dimension %d of block:   %d\n", i, devProp.maxThreadsDim[i]);
    max_n_blocks = devProp.maxGridSize[0];
    for (int i = 0; i < 3; ++i)
      printf("Maximum dimension %d of grid:    %d\n", i, devProp.maxGridSize[i]);
    printf("Clock rate:                     %d\n",  devProp.clockRate);
    printf("Total constant memory:          %lu\n",  devProp.totalConstMem);
    printf("Texture alignment:              %lu\n",  devProp.textureAlignment);
    printf("Concurrent copy and execution:  %s\n",  (devProp.deviceOverlap ? "Yes" : "No"));
    printf("Number of multiprocessors:      %d\n",  devProp.multiProcessorCount);
    printf("Kernel execution timeout:       %s\n",  (devProp.kernelExecTimeoutEnabled ?"Yes" : "No"));
    printf("Memory Clock Rate (KHz):        %d\n", devProp.memoryClockRate);
    printf("Memory Bus Width (bits):        %d\n", devProp.memoryBusWidth);
    printf("Peak Memory Bandwidth (GB/s):   %f\n", 2.0*devProp.memoryClockRate*(devProp.memoryBusWidth/8)/1.0e6);
    printf("Concurrent kernels:             %d\n",  devProp.concurrentKernels);
  }
  size_t available_memory, total_memory;
  cudaMemGetInfo(&available_memory, &total_memory);
  size_t stack_memory;
  cudaDeviceGetLimit(&stack_memory, cudaLimitStackSize);
  size_t fifo_memory;
  cudaDeviceGetLimit(&fifo_memory, cudaLimitPrintfFifoSize);
  size_t heap_memory;
  cudaDeviceGetLimit(&heap_memory, cudaLimitMallocHeapSize);
  printf(" memgetinfo: available_memory %f MB, total_memory %f MB, stack_memory %f MB, fifo_memory %f MB, heap_memory %f MB \n", (double)available_memory/1.e6, (double)total_memory/1.e6, (double)stack_memory/1.e6, (double)fifo_memory/1.e6, (double)heap_memory/1.e6);


  return;
}


__global__ void kernel_find_vertex_with_max_npmts_in_timebin(unsigned int * np, unsigned int * mnp, unsigned int * vmnp){


  // get unique id for each thread in each block == time bin
  unsigned int time_bin_index = threadIdx.x + blockDim.x*blockIdx.x;

  // skip if thread is assigned to nonexistent time bin
  if( time_bin_index >= constant_n_time_bins ) return;


  unsigned int number_of_pmts_in_time_bin = 0;
  unsigned int time_index;
  unsigned int max_number_of_pmts_in_time_bin=0;
  unsigned int vertex_with_max_n_pmts = 0;

  for(unsigned int iv=0;iv<constant_n_test_vertices;iv++) { // loop over test vertices
    // sum the number of hit PMTs in this time window and the next
    
    time_index = time_bin_index + constant_n_time_bins*iv;
    if( time_index >= constant_n_time_bins*constant_n_test_vertices - 1 ) continue;
    number_of_pmts_in_time_bin = np[time_index] + np[time_index+1];
    if( number_of_pmts_in_time_bin >= max_number_of_pmts_in_time_bin ){
      max_number_of_pmts_in_time_bin = number_of_pmts_in_time_bin;
      vertex_with_max_n_pmts = iv;
    }
  }

  mnp[time_bin_index] = max_number_of_pmts_in_time_bin;
  vmnp[time_bin_index] = vertex_with_max_n_pmts;

  return;

}

__global__ void kernel_find_vertex_with_max_npmts_in_timebin_and_directionbin(unsigned int * np, unsigned int * mnp, unsigned int * vmnp){


  // get unique id for each thread in each block == time bin
  unsigned int time_bin_index = threadIdx.x + blockDim.x*blockIdx.x;

  // skip if thread is assigned to nonexistent time bin
  if( time_bin_index >= constant_n_time_bins - 1 ) return;


  unsigned int number_of_pmts_in_time_bin = 0;
  unsigned int max_number_of_pmts_in_time_bin=0;
  unsigned int vertex_with_max_n_pmts = 0;
  unsigned int dir_index_1, dir_index_2;

  for(unsigned int iv=0;iv<constant_n_test_vertices;iv++) { // loop over test vertices
    // sum the number of hit PMTs in this time window
    
    for(unsigned int idir = 0; idir < constant_n_direction_bins; idir++){

      dir_index_1 = device_get_direction_index_at_time(time_bin_index, iv, idir);
      dir_index_2 = device_get_direction_index_at_time(time_bin_index + 1, iv, idir);

      number_of_pmts_in_time_bin = np[dir_index_1]
	+ np[dir_index_2];
      if( number_of_pmts_in_time_bin > max_number_of_pmts_in_time_bin ){
	max_number_of_pmts_in_time_bin = number_of_pmts_in_time_bin;
	vertex_with_max_n_pmts = iv;
      }
    }
  }
  
  mnp[time_bin_index] = max_number_of_pmts_in_time_bin;
  vmnp[time_bin_index] = vertex_with_max_n_pmts;

  return;

}

void free_event_memories(){

  checkCudaErrors(cudaUnbindTexture(tex_ids));
  checkCudaErrors(cudaUnbindTexture(tex_times));
  free(host_ids);
  free(host_times);
  checkCudaErrors(cudaFree(device_ids));
  checkCudaErrors(cudaFree(device_times));
  if( correct_mode == 1 ){
    checkCudaErrors(cudaFree(device_time_bin_of_hit));
  } else if( correct_mode == 2 ){
    checkCudaErrors(cudaFree(device_time_bin_of_hit));
    free(host_time_bin_of_hit);
    free(host_n_pmts_per_time_bin);
  } else if( correct_mode == 3 ){
    checkCudaErrors(cudaFree(device_time_bin_of_hit));
  } else if( correct_mode == 4 ){
    checkCudaErrors(cudaFree(device_time_bin_of_hit));
  } else if( correct_mode == 5 ){
    checkCudaErrors(cudaFree(device_time_bin_of_hit));
  } else if( correct_mode == 6 ){
    checkCudaErrors(cudaFree(device_time_bin_of_hit));
  } else if( correct_mode == 7 ){
    checkCudaErrors(cudaFree(device_time_bin_of_hit));
  }
  if( correct_mode != 9 ){
    checkCudaErrors(cudaFree(device_n_pmts_per_time_bin));
  }else{
    checkCudaErrors(cudaFree(device_n_pmts_per_time_bin_and_direction_bin));
  }
  free(host_max_number_of_pmts_in_time_bin);
  free(host_vertex_with_max_n_pmts);
  checkCudaErrors(cudaFree(device_max_number_of_pmts_in_time_bin));
  checkCudaErrors(cudaFree(device_vertex_with_max_n_pmts));

  return;
}


void free_event_memories_nhits(){

  checkCudaErrors(cudaUnbindTexture(tex_ids));
  checkCudaErrors(cudaUnbindTexture(tex_times));
  free(host_ids);
  free(host_times);
  checkCudaErrors(cudaFree(device_ids));
  checkCudaErrors(cudaFree(device_times));
  checkCudaErrors(cudaFree(device_n_pmts_nhits));
  //  checkCudaErrors(cudaFree(device_time_nhits));
  free(host_n_pmts_nhits);
  //  free(host_time_nhits);

  return;
}

void free_global_memories(){

  //unbind texture reference to free resource 
  checkCudaErrors(cudaUnbindTexture(tex_times_of_flight));

  if( correct_mode == 9 ){
    //    checkCudaErrors(cudaUnbindTexture(tex_directions_for_vertex_and_pmt));
    checkCudaErrors(cudaFree(device_directions_for_vertex_and_pmt));
    free(host_directions_for_vertex_and_pmt);
  }

  free(PMT_x);
  free(PMT_y);
  free(PMT_z);
  free(vertex_x);
  free(vertex_y);
  free(vertex_z);
  checkCudaErrors(cudaFree(device_times_of_flight));
  free(host_times_of_flight);

  return;
}

void copy_candidates_from_device_to_host(){

  checkCudaErrors(cudaMemcpy(host_max_number_of_pmts_in_time_bin,
			     device_max_number_of_pmts_in_time_bin,
			     n_time_bins*sizeof(unsigned int),
			     cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(host_vertex_with_max_n_pmts,
			     device_vertex_with_max_n_pmts,
			     n_time_bins*sizeof(unsigned int),
			     cudaMemcpyDeviceToHost));


}


void choose_candidates_above_threshold(){

  candidate_trigger_pair_vertex_time.clear();
  candidate_trigger_npmts_in_time_bin.clear();

  unsigned int the_threshold;

  for(unsigned int time_bin = 0; time_bin<n_time_bins - 1; time_bin++){ // loop over time bins
    // n_time_bins - 1 as we are checking the i and i+1 at the same time
    
    if(host_vertex_with_max_n_pmts[time_bin] < n_water_like_test_vertices )
      the_threshold = water_like_threshold_number_of_pmts;
    else
      the_threshold = wall_like_threshold_number_of_pmts;

    if(host_max_number_of_pmts_in_time_bin[time_bin] > the_threshold) {

      if( use_verbose )
	printf(" time %f vertex (%f, %f, %f) npmts %d \n", (time_bin + 2)*time_step_size - time_offset, vertex_x[host_vertex_with_max_n_pmts[time_bin]], vertex_y[host_vertex_with_max_n_pmts[time_bin]], vertex_z[host_vertex_with_max_n_pmts[time_bin]], host_max_number_of_pmts_in_time_bin[time_bin]);

      candidate_trigger_pair_vertex_time.push_back(std::make_pair(host_vertex_with_max_n_pmts[time_bin],time_bin+1));
      candidate_trigger_npmts_in_time_bin.push_back(host_max_number_of_pmts_in_time_bin[time_bin]);
    }

  }

  if( use_verbose )
    printf(" n candidates: %d \n", candidate_trigger_pair_vertex_time.size());
}

bool set_input_file_for_event(int n){

  int nchar = (ceil(log10(n+1))+1);
  char * num =  (char*)malloc(sizeof(char)*nchar);
  sprintf(num, "%d", n+1);
  event_file = event_file_base + num + event_file_suffix;

  bool file_exists = (access( event_file.c_str(), F_OK ) != -1);

  free(num);

  return file_exists;

}

void set_output_file(){

  int nchar = (ceil(log10(water_like_threshold_number_of_pmts))+1);
  char * num =  (char*)malloc(sizeof(char)*nchar);
  sprintf(num, "%d", water_like_threshold_number_of_pmts);
  output_file = output_file_base + num + event_file_suffix;

  free(num);

  return ;

}

void set_output_file_nhits(unsigned int threshold){

  int nchar = (ceil(log10(threshold))+1);
  char * num =  (char*)malloc(sizeof(char)*nchar);
  sprintf(num, "%d", threshold);
  output_file = output_file_base + num + event_file_suffix;

  free(num);

  return ;

}

void write_output_nhits(unsigned int * ntriggers){

  if( output_txt ){

    for(unsigned int u=nhits_threshold_min; u<=nhits_threshold_max; u++){
      set_output_file_nhits(u);
      FILE *of=fopen(output_file.c_str(), "a");
      
      //    int trigger = (ntriggers[u - nhits_threshold_min] > 0 ? 1 : 0);
      int trigger = ntriggers[u - nhits_threshold_min];
      fprintf(of, " %d \n", trigger);
      
      fclose(of);
    }

  }


}


void write_output(){

  if( output_txt ){
    FILE *of=fopen(output_file.c_str(), "a");

    // output 1 if there is a trigger, 0 otherwise
    //int trigger = (trigger_pair_vertex_time.size() > 0 ? 1 : 0);

    // output the n of triggers
    int trigger = trigger_pair_vertex_time.size();

    // output the n of water-like triggers
    /* int trigger = 0; */
    /* for(std::vector<std::pair<unsigned int,unsigned int> >::const_iterator itrigger=trigger_pair_vertex_time.begin(); itrigger != trigger_pair_vertex_time.end(); ++itrigger){ */
    /*   if( itrigger->first  < n_water_like_test_vertices ) */
    /* 	trigger ++; */
    /* } */

    fprintf(of, " %d \n", trigger);

    // output reconstructed vertices
    /* for(std::vector<std::pair<unsigned int,unsigned int> >::const_iterator itrigger=trigger_pair_vertex_time.begin(); itrigger != trigger_pair_vertex_time.end(); ++itrigger){ */
    /*   unsigned int triggertime = itrigger->second*time_step_size - time_offset; */
    /*   fprintf(of, " %d %f %f %f %d \n", n_events, vertex_x[itrigger->first], vertex_y[itrigger->first], vertex_z[itrigger->first], triggertime); */
    /* } */


    /* // output non-corrected and corrected times for best vertex */
    /* int max_n_pmts = 0; */
    /* unsigned int best_vertex; */
    /* for(std::vector<std::pair<unsigned int,unsigned int> >::const_iterator itrigger=trigger_pair_vertex_time.begin(); itrigger != trigger_pair_vertex_time.end(); ++itrigger){ */
    /*   unsigned int vertex_index = itrigger->first; */
    /*   unsigned int time_index = itrigger->second; */
    /*   unsigned int local_n_pmts = host_max_number_of_pmts_in_time_bin[itrigger->second]; */
    /*   if( local_n_pmts > max_n_pmts ){ */
    /* 	max_n_pmts = local_n_pmts; */
    /* 	best_vertex = vertex_index; */
    /*   } */
    /* } */
    /* unsigned int distance_index; */
    /* double tof; */
    /* double corrected_time; */
    
    /* for(unsigned int i=0; i<n_hits; i++){ */
      
    /*   distance_index = get_distance_index(host_ids[i], n_PMTs*best_vertex); */
    /*   tof = host_times_of_flight[distance_index]; */
    /*   corrected_time = host_times[i]-tof; */
      
    /*   fprintf(of, " %d %d %f \n", host_ids[i], host_times[i], corrected_time); */
    /*   //fprintf(of, " %d %f \n", host_ids[i], corrected_time); */
    /* } */

    
    fclose(of);
  }


}


void initialize_output(){

  if( output_txt )
    remove( output_file.c_str() );

}

void initialize_output_nhits(){

  if( output_txt )
    for(unsigned int u=nhits_threshold_min; u<=nhits_threshold_max; u++){
      set_output_file_nhits(u);
      remove( output_file.c_str() );
    }
}

float read_value_from_file(std::string paramname, std::string filename){

  FILE * pFile = fopen (filename.c_str(),"r");

  char name[256];
  float value;

  while( EOF != fscanf(pFile, "%s %e", name, &value) ){
    if( paramname.compare(name) == 0 ){
      fclose(pFile);
      return value;
    }
  }

  printf(" warning: could not find parameter %s in file %s \n", paramname.c_str(), filename.c_str());

  fclose(pFile);
  exit(0);

  return 0.;

}

void read_user_parameters(){

  std::string parameter_file = "parameters.txt";

  twopi = 2.*acos(-1.);
  speed_light_water = 29.9792/1.3330; // speed of light in water, cm/ns
  //double speed_light_water = 22.490023;

  cerenkov_angle_water = acos(1./1.3330);

  dark_rate = read_value_from_file("dark_rate", parameter_file); // Hz
  cylindrical_grid = (bool)read_value_from_file("cylindrical_grid", parameter_file);
  distance_between_vertices = read_value_from_file("distance_between_vertices", parameter_file); // cm
  wall_like_distance = read_value_from_file("wall_like_distance", parameter_file); // units of distance between vertices
  time_step_size = (unsigned int)(sqrt(3.)*distance_between_vertices/(4.*speed_light_water)); // ns
  int extra_threshold = (int)(dark_rate*n_PMTs*2.*time_step_size*1.e-9); // to account for dark current occupancy
  water_like_threshold_number_of_pmts = read_value_from_file("water_like_threshold_number_of_pmts", parameter_file) + extra_threshold;
  wall_like_threshold_number_of_pmts = read_value_from_file("wall_like_threshold_number_of_pmts", parameter_file) + extra_threshold;
  coalesce_time = read_value_from_file("coalesce_time", parameter_file); // ns
  trigger_gate_up = read_value_from_file("trigger_gate_up", parameter_file); // ns
  trigger_gate_down = read_value_from_file("trigger_gate_down", parameter_file); // ns
  max_n_hits_per_job = read_value_from_file("max_n_hits_per_job", parameter_file);
  output_txt = (bool)read_value_from_file("output_txt", parameter_file);
  correct_mode = read_value_from_file("correct_mode", parameter_file);
  number_of_kernel_blocks_3d.y = read_value_from_file("num_blocks_y", parameter_file);
  number_of_threads_per_block_3d.y = read_value_from_file("num_threads_per_block_y", parameter_file);
  number_of_threads_per_block_3d.x = read_value_from_file("num_threads_per_block_x", parameter_file);

  n_direction_bins_theta = read_value_from_file("n_direction_bins_theta", parameter_file);
  n_direction_bins_phi = 2*(n_direction_bins_theta - 1);
  n_direction_bins = n_direction_bins_phi*n_direction_bins_theta - 2*(n_direction_bins_phi - 1);


}


void read_user_parameters_nhits(){

  std::string parameter_file = "parameters.txt";

  twopi = 2.*acos(-1.);
  speed_light_water = 29.9792/1.3330; // speed of light in water, cm/ns
  //double speed_light_water = 22.490023;

  double dark_rate = read_value_from_file("dark_rate", parameter_file); // Hz
  distance_between_vertices = read_value_from_file("distance_between_vertices", parameter_file); // cm
  wall_like_distance = read_value_from_file("wall_like_distance", parameter_file); // units of distance between vertices
  time_step_size = read_value_from_file("nhits_step_size", parameter_file); // ns
  nhits_window = read_value_from_file("nhits_window", parameter_file); // ns
  int extra_threshold = (int)(dark_rate*n_PMTs*nhits_window*1.e-9); // to account for dark current occupancy
  extra_threshold = 0;
  water_like_threshold_number_of_pmts = read_value_from_file("water_like_threshold_number_of_pmts", parameter_file) + extra_threshold;
  wall_like_threshold_number_of_pmts = read_value_from_file("wall_like_threshold_number_of_pmts", parameter_file) + extra_threshold;
  coalesce_time = read_value_from_file("coalesce_time", parameter_file); // ns
  trigger_gate_up = read_value_from_file("trigger_gate_up", parameter_file); // ns
  trigger_gate_down = read_value_from_file("trigger_gate_down", parameter_file); // ns
  max_n_hits_per_job = read_value_from_file("max_n_hits_per_job", parameter_file);
  output_txt = (bool)read_value_from_file("output_txt", parameter_file);
  correct_mode = read_value_from_file("correct_mode", parameter_file);
  number_of_kernel_blocks_3d.y = read_value_from_file("num_blocks_y", parameter_file);
  number_of_threads_per_block_3d.y = read_value_from_file("num_threads_per_block_y", parameter_file);
  number_of_threads_per_block_3d.x = read_value_from_file("num_threads_per_block_x", parameter_file);
  nhits_threshold_min = read_value_from_file("nhits_threshold_min", parameter_file);
  nhits_threshold_max = read_value_from_file("nhits_threshold_max", parameter_file);

}


void check_cudamalloc_float(unsigned int size){

  unsigned int bytes_per_float = 4;
  size_t available_memory, total_memory;
  cudaMemGetInfo(&available_memory, &total_memory);
  if( size*bytes_per_float > available_memory*1000/1024 ){
    printf(" cannot allocate %d floats, or %d B, available %d B \n", 
	   size, size*bytes_per_float, available_memory*1000/1024);
  }

}

void check_cudamalloc_int(unsigned int size){

  unsigned int bytes_per_int = 4;
  size_t available_memory, total_memory;
  cudaMemGetInfo(&available_memory, &total_memory);
  if( size*bytes_per_int > available_memory*1000/1024 ){
    printf(" cannot allocate %d ints, or %d B, available %d B \n", 
	   size, size*bytes_per_int, available_memory*1000/1024);
  }

}

void check_cudamalloc_unsigned_int(unsigned int size){

  unsigned int bytes_per_unsigned_int = 4;
  size_t available_memory, total_memory;
  cudaMemGetInfo(&available_memory, &total_memory);
  if( size*bytes_per_unsigned_int > available_memory*1000/1024 ){
    printf(" cannot allocate %d unsigned_ints, or %d B, available %d B \n", 
	   size, size*bytes_per_unsigned_int, available_memory*1000/1024);
  }

}


void check_cudamalloc_bool(unsigned int size){

  unsigned int bytes_per_bool = 1;
  size_t available_memory, total_memory;
  cudaMemGetInfo(&available_memory, &total_memory);
  if( size*bytes_per_bool > available_memory*1000/1024 ){
    printf(" cannot allocate %d unsigned_ints, or %d B, available %d B \n", 
	   size, size*bytes_per_bool, available_memory*1000/1024);
  }

}

unsigned int find_greatest_divisor(unsigned int n, unsigned int max){

  if( n == 1 ){
    return 1;
  }

  if (n % 2 == 0){
    if( n <= 2*max ){
      return n / 2;
    } 
    else{
      float sqrtN = sqrt(n); // square root of n in float precision.
      unsigned int start = ceil(std::max((double)2,(double)n/(double)max));
      for(unsigned int d = start; d <= n; d += 1)
	if (n % d == 0)
	  return n/d;
      return 1;
    }
  }

  // Now, the least prime divisor of n is odd.
  // So, we increment the counter by 2 in the loop, by starting in 3.
  
  float sqrtN = sqrt(n); // square root of n in float precision.
  unsigned int start = ceil(std::max((double)3,(double)n/(double)max));
  for(unsigned int d = start; d <= n; d += 2)
    if (n % d == 0)
      return n/d;
  
  // If the loop has reached its end normally, 
  // it means that N is prime.

  return 1;

}


void setup_threads_for_histo(unsigned int n){

  number_of_kernel_blocks_3d.x = n/greatest_divisor;
  number_of_kernel_blocks_3d.y = 1;

  number_of_threads_per_block_3d.x = greatest_divisor;
  number_of_threads_per_block_3d.y = 1;

}

void setup_threads_for_histo(){

  number_of_kernel_blocks_3d.x = 1000;
  number_of_kernel_blocks_3d.y = 1;

  number_of_threads_per_block_3d.x = max_n_threads_per_block;
  number_of_threads_per_block_3d.y = 1;

}

void setup_threads_for_histo_iterated(bool last){

  number_of_kernel_blocks_3d.x = 1000;
  number_of_kernel_blocks_3d.y = 1;

  unsigned int size = n_time_bins*n_test_vertices;
  number_of_threads_per_block_3d.x = (last ? size - size/max_n_threads_per_block*max_n_threads_per_block : max_n_threads_per_block);
  number_of_threads_per_block_3d.y = 1;

}

void setup_threads_for_histo_per(unsigned int n){

  number_of_kernel_blocks_3d.x = n;

  print_parameters_2d();

  if( number_of_threads_per_block_3d.x * number_of_threads_per_block_3d.y > max_n_threads_per_block ){
    printf(" --------------------- warning: number_of_threads_per_block (x*y) = %d cannot exceed max value %d \n", number_of_threads_per_block_3d.x * number_of_threads_per_block_3d.y, max_n_threads_per_block );
  }

  if( number_of_kernel_blocks_3d.x > max_n_blocks ){
    printf(" warning: number_of_kernel_blocks x = %d cannot exceed max value %d \n", number_of_kernel_blocks_3d.x, max_n_blocks );
  }

  if( number_of_kernel_blocks_3d.y > max_n_blocks ){
    printf(" warning: number_of_kernel_blocks y = %d cannot exceed max value %d \n", number_of_kernel_blocks_3d.y, max_n_blocks );
  }

  if( std::numeric_limits<int>::max() / (number_of_kernel_blocks_3d.x*number_of_kernel_blocks_3d.y)  < number_of_threads_per_block_3d.x*number_of_threads_per_block_3d.y ){
    printf(" --------------------- warning: grid size cannot exceed max value %u \n", std::numeric_limits<int>::max() );
  }

}


#endif

