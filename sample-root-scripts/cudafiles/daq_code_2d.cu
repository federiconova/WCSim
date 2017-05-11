
//
// include files
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <helper_cuda.h>
#include <sys/time.h>
#include <library_daq.h>

// CUDA = Computer Device Unified Architecture

__global__ void kernel_correct_times_and_get_n_pmts_per_time_bin(unsigned int *ct);
__global__ void kernel_correct_times_and_get_n_pmts_per_time_bin_and_direction_bin(unsigned int *ct, bool * dirs);
__global__ void kernel_correct_times(unsigned int *ct);
__global__ void kernel_histo_one_thread_one_vertex( unsigned int *ct, unsigned int *histo );
__global__ void kernel_histo_stride( unsigned int *ct, unsigned int *histo);
__global__ void kernel_histo_iterated( unsigned int *ct, unsigned int *histo, unsigned int offset );
__device__ int get_time_bin();
__device__ int get_time_bin_for_vertex_and_hit(unsigned int vertex_index, unsigned int hit_index);
__global__ void kernel_histo_stride_2d( unsigned int *ct, unsigned int *histo);
__global__ void kernel_histo_per_vertex( unsigned int *ct, unsigned int *histo);
__global__ void kernel_histo_per_vertex_shared( unsigned int *ct, unsigned int *histo);
__global__ void kernel_correct_times_and_get_histo_per_vertex_shared(unsigned int *ct);

int gpu_daq_initialize();
int gpu_daq_execute();
int gpu_daq_finalize();


//
// main code
//

int main()
{

  bool good_init = (bool)gpu_daq_initialize();

  bool good_exec = (bool)gpu_daq_execute();

  bool good_fin = (bool)gpu_daq_finalize();



  return 1;
}



//
// kernel routine
// 

// __global__ identifier says it's a kernel function
__global__ void kernel_correct_times_and_get_n_pmts_per_time_bin(unsigned int *ct){

  int time_bin = get_time_bin();

  if( time_bin < 0 ) return;

  atomicAdd(&ct[time_bin],1);

  //  printf( " hit %d (nh %d) id %d t %d; vertex %d (nv %d) tof %f  %d \n", hit_index, constant_n_hits, ids[hit_index], t[hit_index], vertex_index, constant_n_test_vertices, tof, ct[time_index]);

  return;

}


// __global__ identifier says it's a kernel function
__global__ void kernel_correct_times_and_get_n_pmts_per_time_bin_and_direction_bin(unsigned int *ct, bool * dirs){

  // get unique id for each thread in each block
  unsigned int tid_x = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int tid_y = threadIdx.y + blockDim.y*blockIdx.y;

  // map the two 2D indices to a single linear, 1D index
  int tid = tid_y * gridDim.x * blockDim.x + tid_x;

  // tid runs from 0 to n_test_vertices * n_hits:
  //      vertex 0           vertex 1       ...     vertex m
  // (hit 0, ..., hit n; hit 0, ..., hit n; ...; hit 0, ..., hit n);

  unsigned int vertex_index = (int)(tid/constant_n_hits);
  unsigned int hit_index = tid - vertex_index * constant_n_hits;

  // skip if thread is assigned to nonexistent vertex
  if( vertex_index >= constant_n_test_vertices ) return ;

  // skip if thread is assigned to nonexistent hit
  if( hit_index >= constant_n_hits ) return ;

  int time_bin = get_time_bin_for_vertex_and_hit(vertex_index, hit_index) - constant_n_time_bins*vertex_index;

  if( time_bin < 0 ) return;

  for(unsigned int idir = 0; idir < constant_n_direction_bins; idir++){

    unsigned int dir_index = device_get_direction_index_at_pmt(
							       tex1Dfetch(tex_ids,hit_index), 
							       vertex_index, 
							       idir
							       );
    
    //    bool good_direction = (bool)tex1Dfetch(tex_directions_for_vertex_and_pmt, dir_index);

    bool good_direction = dirs[dir_index];


    if( good_direction ){
      atomicAdd(&ct[device_get_direction_index_at_time(time_bin, vertex_index, idir)],1);
    }
    
  }





  //  printf( " hit %d (nh %d) id %d t %d; vertex %d (nv %d) tof %f  %d \n", hit_index, constant_n_hits, ids[hit_index], t[hit_index], vertex_index, constant_n_test_vertices, tof, ct[time_index]);

  return;

}


__global__ void kernel_correct_times(unsigned int *ct){


  int time_bin = get_time_bin();

  if( time_bin < 0 ) return;

  // get unique id for each thread in each block
  unsigned int tid_x = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int tid_y = threadIdx.y + blockDim.y*blockIdx.y;

  // map the two 2D indices to a single linear, 1D index
  int tid = tid_y * gridDim.x * blockDim.x + tid_x;

  ct[tid] = time_bin;


  return;

}

__device__ int get_time_bin_for_vertex_and_hit(unsigned int vertex_index, unsigned int hit_index){

  // skip if thread is assigned to nonexistent vertex
  if( vertex_index >= constant_n_test_vertices ) return -1;

  // skip if thread is assigned to nonexistent hit
  if( hit_index >= constant_n_hits ) return -1;

  unsigned int vertex_block = constant_n_time_bins*vertex_index;

  unsigned int vertex_block2 = constant_n_PMTs*vertex_index;

  return device_get_time_index(
			       (tex1Dfetch(tex_times,hit_index)
				- tex1Dfetch(tex_times_of_flight,
					     device_get_distance_index(
								       tex1Dfetch(tex_ids,hit_index),
								       vertex_block2
								       )
					     )
				+ constant_time_offset)/constant_time_step_size
			       ,
			       vertex_block
			       );
  

}


__device__ int get_time_bin(){

  // get unique id for each thread in each block
  unsigned int tid_x = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int tid_y = threadIdx.y + blockDim.y*blockIdx.y;

  // map the two 2D indices to a single linear, 1D index
  int tid = tid_y * gridDim.x * blockDim.x + tid_x;

  // tid runs from 0 to n_test_vertices * n_hits:
  //      vertex 0           vertex 1       ...     vertex m
  // (hit 0, ..., hit n; hit 0, ..., hit n; ...; hit 0, ..., hit n);

  unsigned int vertex_index = (int)(tid/constant_n_hits);
  unsigned int hit_index = tid - vertex_index * constant_n_hits;

  return get_time_bin_for_vertex_and_hit(vertex_index, hit_index);


}






__global__ void kernel_histo_one_thread_one_vertex( unsigned int *ct, unsigned int *histo ){

  
  // get unique id for each thread in each block
  unsigned int tid_x = threadIdx.x + blockDim.x*blockIdx.x;

  unsigned int vertex_index = tid_x;
  unsigned int bin ;
  unsigned int max = constant_n_test_vertices*constant_n_hits;
  unsigned int size = vertex_index * constant_n_hits;

  for( unsigned int ihit=0; ihit<constant_n_hits; ihit++){
    bin = size + ihit;
    if( bin < max)
      atomicAdd(&histo[ct[bin]],1);
  }
  
}

__global__ void kernel_histo_stride( unsigned int *ct, unsigned int *histo){

  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  while( i < constant_n_hits*constant_n_test_vertices ){
    atomicAdd( &histo[ct[i]], 1);
    i += stride;
  }


}



__global__ void kernel_histo_iterated( unsigned int *ct, unsigned int *histo, unsigned int offset ){

  
  extern __shared__ unsigned int temp[];
  unsigned int index = threadIdx.x + offset;
  temp[index] = 0;
  __syncthreads();
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int size = blockDim.x * gridDim.x;
  unsigned int max = constant_n_hits*constant_n_test_vertices;
  while( i < max ){
    atomicAdd( &temp[ct[i]], 1);
    i += size;
  }
  __syncthreads();
  atomicAdd( &(histo[index]), temp[index] );


}


__global__ void kernel_histo_stride_2d( unsigned int *ct, unsigned int *histo){

  // get unique id for each thread in each block
  unsigned int tid_x = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int tid_y = threadIdx.y + blockDim.y*blockIdx.y;

  unsigned int size = blockDim.x * gridDim.x;
  unsigned int max = constant_n_hits*constant_n_test_vertices;

  // map the two 2D indices to a single linear, 1D index
  int tid = tid_y * size + tid_x;

  /*
  unsigned int vertex_index = (int)(tid/constant_n_time_bins);
  unsigned int time_index = tid - vertex_index * constant_n_time_bins;

  // skip if thread is assigned to nonexistent vertex
  if( vertex_index >= constant_n_test_vertices ) return;

  // skip if thread is assigned to nonexistent hit
  if( time_index >= constant_n_time_bins ) return;

  unsigned int vertex_block = constant_n_time_bins*vertex_index;

  unsigned int vertex_block2 = constant_n_PMTs*vertex_index;
  */

  unsigned int stride = blockDim.y * gridDim.y*size;

  while( tid < max ){
    atomicAdd( &histo[ct[tid]], 1);
    tid += stride;
  }


}


__global__ void kernel_histo_per_vertex( unsigned int *ct, unsigned int *histo){

  // get unique id for each thread in each block
  unsigned int tid_x = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int tid_y = threadIdx.y + blockDim.y*blockIdx.y;

  if( tid_x >= constant_n_test_vertices ) return;

  unsigned int vertex_offset = tid_x*constant_n_hits;
  unsigned int bin;
  unsigned int stride = blockDim.y*gridDim.y;
  unsigned int ihit = vertex_offset + tid_y;

  while( ihit<vertex_offset+constant_n_hits){

    bin = ct[ihit];
    //histo[bin]++;
    atomicAdd( &histo[bin], 1);
    ihit += stride;

  }
  __syncthreads();
}

__global__ void kernel_histo_per_vertex_shared( unsigned int *ct, unsigned int *histo){
  // get unique id for each thread in each block
  unsigned int tid_x = threadIdx.x + blockDim.x*blockIdx.x;
  unsigned int tid_y = threadIdx.y + blockDim.y*blockIdx.y;

  if( tid_x >= constant_n_test_vertices ) return;

  unsigned int vertex_offset = tid_x*constant_n_hits;
  unsigned int bin;
  unsigned int stride = blockDim.y*gridDim.y;
  unsigned int stride_block = blockDim.y;
  unsigned int ihit = vertex_offset + tid_y;
  unsigned int time_offset = tid_x*constant_n_time_bins;

  unsigned int local_ihit = threadIdx.y;
  extern __shared__ unsigned int temp[];
  while( local_ihit<constant_n_time_bins ){
    temp[local_ihit] = 0;
    local_ihit += stride_block;
  }

  __syncthreads();

  while( ihit<vertex_offset+constant_n_hits){

    bin = ct[ihit];
    atomicAdd(&temp[bin - time_offset],1);
    ihit += stride;

  }

  __syncthreads();

  local_ihit = threadIdx.y;
  while( local_ihit<constant_n_time_bins ){
    atomicAdd( &histo[local_ihit+time_offset], temp[local_ihit]);
    local_ihit += stride_block;
  }


}

__global__ void kernel_correct_times_and_get_histo_per_vertex_shared(unsigned int *ct){

  unsigned int vertex_index = blockIdx.x;
  if( vertex_index >= constant_n_test_vertices ) return;

  unsigned int local_ihit_initial = threadIdx.x + threadIdx.y*blockDim.x;
  unsigned int local_ihit = local_ihit_initial;
  unsigned int stride_block = blockDim.x*blockDim.y;
  unsigned int stride = stride_block*gridDim.y;
  unsigned int hit_index = local_ihit + stride_block*blockIdx.y;

  unsigned int bin;
  unsigned int time_offset = vertex_index*constant_n_time_bins;

  extern __shared__ unsigned int temp[];
  while( local_ihit<constant_n_time_bins ){
    temp[local_ihit] = 0;
    local_ihit += stride_block;
  }

  __syncthreads();

  while( hit_index<constant_n_hits){

    bin = get_time_bin_for_vertex_and_hit(vertex_index, hit_index);
    atomicAdd(&temp[bin - time_offset],1);
    hit_index += stride;

  }

  __syncthreads();

  local_ihit = local_ihit_initial;
  while( local_ihit<constant_n_time_bins ){
    atomicAdd( &ct[local_ihit+time_offset], temp[local_ihit]);
    local_ihit += stride_block;
  }


}


int gpu_daq_initialize(){

  int argc = 0;
  const char* n_argv[] = {};
  const char **argv = n_argv;

  /////////////////////
  // initialise card //
  /////////////////////
  findCudaDevice(argc, argv);


  // initialise CUDA timing
  use_timing = true;
  if( use_timing ){
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  }
  cudaEventCreate(&total_start);
  cudaEventCreate(&total_stop);
  elapsed_parameters = 0; elapsed_pmts = 0; elapsed_detector = 0; elapsed_vertices = 0;
  elapsed_threads = 0; elapsed_tof = 0; elapsed_directions = 0; elapsed_memory_tofs_dev = 0; elapsed_memory_directions_dev = 0; elapsed_memory_candidates_host = 0; elapsed_tofs_copy_dev = 0; elapsed_directions_copy_dev = 0;
  elapsed_input = 0; elapsed_memory_dev = 0; elapsed_copy_dev = 0; elapsed_kernel_correct_times_and_get_n_pmts_per_time_bin = 0; 
  elapsed_threads_candidates = 0; elapsed_candidates_memory_dev = 0; elapsed_candidates_kernel = 0;
  elapsed_candidates_copy_host = 0; choose_candidates = 0; elapsed_coalesce = 0; elapsed_gates = 0; elapsed_free = 0; elapsed_total = 0;
  elapsed_tofs_free = 0; elapsed_reset = 0;
  use_verbose = false;


  ////////////////////
  // inspect device //
  ////////////////////
  // set: max_n_threads_per_block, max_n_blocks
  print_gpu_properties();




  ////////////////
  // read PMTs  //
  ////////////////
  // set: n_PMTs, PMT_x, PMT_y, PMT_z
  if( use_timing )
    start_c_clock();
  event_file_base = "all_hits_";
  event_file_suffix = ".txt";
  detector_file = "detector.txt";
  pmts_file = "all_pmts.txt";
  output_file_base = "all_hits_emerald_threshold_";
  if( !read_the_pmts() ) return 0;
  if( use_timing )
    elapsed_pmts = stop_c_clock();


  ///////////////////////
  // define parameters //
  ///////////////////////
  if( use_timing )
    start_c_clock();
  read_user_parameters();
  set_output_file();
  printf(" --- user parameters \n");
  printf(" dark_rate %f \n", dark_rate);
  printf(" distance between test vertices = %f cm \n", distance_between_vertices);
  printf(" wall_like_distance %f \n", wall_like_distance);
  printf(" water_like_threshold_number_of_pmts = %d \n", water_like_threshold_number_of_pmts);
  printf(" wall_like_threshold_number_of_pmts %d \n", wall_like_threshold_number_of_pmts);
  printf(" coalesce_time = %f ns \n", coalesce_time);
  printf(" trigger_gate_up = %f ns \n", trigger_gate_up);
  printf(" trigger_gate_down = %f ns \n", trigger_gate_down);
  printf(" max_n_hits_per_job = %d \n", max_n_hits_per_job);
  printf(" output_txt %d \n", output_txt);
  printf(" correct_mode %d \n", correct_mode);
  printf(" num_blocks_y %d \n", number_of_kernel_blocks_3d.y);
  printf(" num_threads_per_block_x %d \n", number_of_threads_per_block_3d.x);
  printf(" num_threads_per_block_y %d \n", number_of_threads_per_block_3d.y);
  printf(" cylindrical_grid %d \n", cylindrical_grid);
  printf(" time step size = %d ns \n", time_step_size);
  if( correct_mode == 9 ){
    printf(" n_direction_bins_theta %d, n_direction_bins_phi %d, n_direction_bins %d \n",
	   n_direction_bins_theta, n_direction_bins_phi, n_direction_bins);
  }
  if( use_timing )
    elapsed_parameters = stop_c_clock();




  /////////////////////
  // read detector ////
  /////////////////////
  // set: detector_height, detector_radius, pmt_radius
  if( use_timing )
    start_c_clock();
  if( !read_the_detector() ) return 0;
  if( use_timing )
    elapsed_detector = stop_c_clock();




  ////////////////////////
  // make test vertices //
  ////////////////////////
  // set: n_test_vertices, n_water_like_test_vertices, vertex_x, vertex_y, vertex_z
  // use: detector_height, detector_radius
  if( use_timing )
    start_c_clock();
  make_test_vertices();
  if( use_timing )
    elapsed_vertices = stop_c_clock();



  //////////////////////////////
  // table of times_of_flight //
  //////////////////////////////
  // set: host_times_of_flight, time_offset
  // use: n_test_vertices, vertex_x, vertex_y, vertex_z, n_PMTs, PMT_x, PMT_y, PMT_z
  // malloc: host_times_of_flight
  if( use_timing )
    start_c_clock();
  make_table_of_tofs();
  if( use_timing )
    elapsed_tof = stop_c_clock();



  if( correct_mode == 9 ){
    //////////////////////////////
    // table of directions //
    //////////////////////////////
    // set: host_directions_phi, host_directions_cos_theta
    // use: n_test_vertices, vertex_x, vertex_y, vertex_z, n_PMTs, PMT_x, PMT_y, PMT_z
    // malloc: host_directions_phi, host_directions_cos_theta
    if( use_timing )
      start_c_clock();
    make_table_of_directions();
    if( use_timing )
      elapsed_directions = stop_c_clock();
  }


  ////////////////////////////////////
  // allocate tofs memory on device //
  ////////////////////////////////////
  // use: n_test_vertices, n_PMTs
  // cudamalloc: device_times_of_flight
  if( use_timing )
    start_cuda_clock();
  allocate_tofs_memory_on_device();
  if( use_timing )
    elapsed_memory_tofs_dev = stop_cuda_clock();


  if( correct_mode == 9 ){
    ////////////////////////////////////
    // allocate direction memory on device //
    ////////////////////////////////////
    // use: n_test_vertices, n_PMTs
    // cudamalloc: device_directions_phi, device_directions_cos_theta
    if( use_timing )
      start_cuda_clock();
    allocate_directions_memory_on_device();
    if( use_timing )
      elapsed_memory_directions_dev = stop_cuda_clock();
  }


  ////////////////////////////////
  // fill tofs memory on device //
  ////////////////////////////////
  // use: n_test_vertices, n_water_like_test_vertices, n_PMTs
  // memcpy: device_times_of_flight, constant_time_step_size, constant_n_test_vertices, constant_n_water_like_test_vertices, constant_n_PMTs
  // texture: tex_times_of_flight
  if( use_timing )
    start_cuda_clock();
  fill_tofs_memory_on_device();
  if( use_timing )
    elapsed_tofs_copy_dev = stop_cuda_clock();


  if( correct_mode == 9 ){
    ////////////////////////////////
    // fill directions memory on device //
    ////////////////////////////////
    // use: n_test_vertices, n_water_like_test_vertices, n_PMTs
    // memcpy: device_directions_phi, device_directions_cos_theta, constant_time_step_size, constant_n_test_vertices, constant_n_water_like_test_vertices, constant_n_PMTs
    // texture: tex_directions_phi, tex_directions_cos_theta
    if( use_timing )
      start_cuda_clock();
    fill_directions_memory_on_device();
    if( use_timing )
      elapsed_directions_copy_dev = stop_cuda_clock();
  }

  ///////////////////////
  // initialize output //
  ///////////////////////
  initialize_output();


  return 1;

}

int gpu_daq_execute(){

  start_total_cuda_clock();

  n_events = 0;

  while( set_input_file_for_event(n_events) ){

    printf(" ------ analyzing event %d \n", n_events+1);

    ////////////////
    // read input //
    ////////////////
    // set: n_hits, host_ids, host_times, time_offset, n_time_bins
    // use: time_offset, n_test_vertices
    // memcpy: constant_n_time_bins, constant_n_hits
    if( use_timing )
      start_c_clock();
    if( !read_the_input() ){
      if( use_timing )
	elapsed_input += stop_c_clock();
      write_output();
      n_events ++;
      continue;
    }
    if( use_timing )
      elapsed_input += stop_c_clock();
  


    ////////////////////////////////////////
    // allocate candidates memory on host //
    ////////////////////////////////////////
    // use: n_time_bins
    // malloc: host_max_number_of_pmts_in_time_bin, host_vertex_with_max_n_pmts
    if( use_timing )
      start_cuda_clock();
    allocate_candidates_memory_on_host();
    if( use_timing )
      elapsed_memory_candidates_host += stop_cuda_clock();


    if( correct_mode != 8 ){
      ////////////////////////////////////////////////
      // set number of blocks and threads per block //
      ////////////////////////////////////////////////
      // set: number_of_kernel_blocks, number_of_threads_per_block
      // use: n_test_vertices, n_hits
      if( use_timing )
	start_c_clock();
      if( !setup_threads_for_tof_2d(n_test_vertices, n_hits) ) return 0;
      if( use_timing )
	elapsed_threads += stop_c_clock();
    }


    ///////////////////////////////////////
    // allocate correct memory on device //
    ///////////////////////////////////////
    // use: n_test_vertices, n_hits, n_time_bins
    // cudamalloc: device_ids, device_times, device_n_pmts_per_time_bin
    if( use_timing )
      start_cuda_clock();
    allocate_correct_memory_on_device();
    if( use_timing )
      elapsed_memory_dev += stop_cuda_clock();


    //////////////////////////////////////
    // copy input into device variables //
    //////////////////////////////////////
    // use: n_hits
    // memcpy: device_ids, device_times, constant_time_offset
    // texture: tex_ids, tex_times
    if( use_timing )
      start_cuda_clock();
    fill_correct_memory_on_device();
    if( use_timing )
      elapsed_copy_dev += stop_cuda_clock();



    ////////////////////
    // execute kernel //
    ////////////////////
    if( use_timing )
      start_cuda_clock();
    if( correct_mode == 0 ){
      printf(" --- execute kernel to correct times and get n pmts per time bin \n");
      kernel_correct_times_and_get_n_pmts_per_time_bin<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_n_pmts_per_time_bin);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");
    }else if( correct_mode == 1 ){
      printf(" --- execute kernel to correct times \n");
      kernel_correct_times<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");

      setup_threads_for_histo(n_test_vertices);
      printf(" --- execute kernel to get n pmts per time bin \n");
      kernel_histo_one_thread_one_vertex<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit, device_n_pmts_per_time_bin);
      cudaThreadSynchronize();
      getLastCudaError("kernel_histo_one_thread_one_vertex execution failed\n");
    }else if( correct_mode == 2 ){
      printf(" --- execute kernel to correct times \n");
      kernel_correct_times<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");

      checkCudaErrors(cudaMemcpy(host_time_bin_of_hit,
				 device_time_bin_of_hit,
				 n_hits*n_test_vertices*sizeof(unsigned int),
				 cudaMemcpyDeviceToHost));

      for( unsigned int u=0; u<n_hits*n_test_vertices; u++){
	unsigned int bin = host_time_bin_of_hit[u];
	if( bin < n_time_bins*n_test_vertices )
	  host_n_pmts_per_time_bin[ bin ] ++;
      }

      checkCudaErrors(cudaMemcpy(device_n_pmts_per_time_bin,
				 host_n_pmts_per_time_bin,
				 n_time_bins*n_test_vertices*sizeof(unsigned int),
				 cudaMemcpyHostToDevice));
      

    }else if( correct_mode == 3 ){
      printf(" --- execute kernel to correct times \n");
      kernel_correct_times<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");
      
      setup_threads_for_histo();
      printf(" --- execute kernel to get n pmts per time bin \n");
      kernel_histo_stride<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit, device_n_pmts_per_time_bin);
      cudaThreadSynchronize();
      getLastCudaError("kernel_histo_one_thread_one_vertex execution failed\n");
    }else if( correct_mode == 4 ){
      printf(" --- execute kernel to correct times \n");
      kernel_correct_times<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");

      unsigned int njobs = n_time_bins*n_test_vertices/max_n_threads_per_block + 1;
      printf(" executing %d njobs to get n pmts per time bin \n", njobs); 
      for( unsigned int iter=0; iter<njobs; iter++){

	setup_threads_for_histo_iterated((bool)(iter + 1 == njobs));

	kernel_histo_iterated<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d,n_time_bins*n_test_vertices*sizeof(unsigned int) >>>(device_time_bin_of_hit, device_n_pmts_per_time_bin, iter*max_n_threads_per_block);
	cudaThreadSynchronize();
	getLastCudaError("kernel_histo execution failed\n");
      }

    }else if( correct_mode == 5 ){
      printf(" --- execute kernel to correct times \n");
      kernel_correct_times<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");

      if( !setup_threads_for_tof_2d(n_test_vertices, n_time_bins) ) return 0;

      printf(" executing kernel to get n pmts per time bin \n"); 

      kernel_histo_stride_2d<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d >>>(device_time_bin_of_hit, device_n_pmts_per_time_bin);
      cudaThreadSynchronize();
      getLastCudaError("kernel_histo execution failed\n");

    }else if( correct_mode == 6 ){
      printf(" --- execute kernel to correct times \n");
      kernel_correct_times<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");
      
      setup_threads_for_histo_per(n_test_vertices);
      printf(" --- execute kernel to get n pmts per time bin \n");
      kernel_histo_per_vertex<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit, device_n_pmts_per_time_bin);
      cudaThreadSynchronize();
      getLastCudaError("kernel_histo_one_thread_one_vertex execution failed\n");
    }else if( correct_mode == 7 ){
      printf(" --- execute kernel to correct times \n");
      kernel_correct_times<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_time_bin_of_hit);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");
      
      setup_threads_for_histo_per(n_test_vertices);
      printf(" --- execute kernel to get n pmts per time bin \n");
      kernel_histo_per_vertex_shared<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d,n_time_bins*sizeof(unsigned int)>>>(device_time_bin_of_hit, device_n_pmts_per_time_bin);
      cudaThreadSynchronize();
      getLastCudaError("kernel_histo_one_thread_one_vertex execution failed\n");
    }else if( correct_mode == 8 ){
      setup_threads_for_histo_per(n_test_vertices);
      printf(" --- execute kernel to correct times and get n pmts per time bin \n");
      kernel_correct_times_and_get_histo_per_vertex_shared<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d,n_time_bins*sizeof(unsigned int)>>>(device_n_pmts_per_time_bin);
      cudaThreadSynchronize();
      getLastCudaError("kernel_correct_times_and_get_histo_per_vertex_shared execution failed\n");
    }else if( correct_mode == 9 ){
      printf(" --- execute kernel to correct times and get n pmts per time bin \n");
      kernel_correct_times_and_get_n_pmts_per_time_bin_and_direction_bin<<<number_of_kernel_blocks_3d,number_of_threads_per_block_3d>>>(device_n_pmts_per_time_bin_and_direction_bin, device_directions_for_vertex_and_pmt);
      cudaThreadSynchronize();
      getLastCudaError("correct_kernel execution failed\n");
    }
    if( use_timing )
      elapsed_kernel_correct_times_and_get_n_pmts_per_time_bin += stop_cuda_clock();



    //////////////////////////////////
    // setup threads for candidates //
    //////////////////////////////////
    // set: number_of_kernel_blocks, number_of_threads_per_block
    // use: n_time_bins
    if( use_timing )
      start_c_clock();
    if( !setup_threads_to_find_candidates() ) return 0;
    if( use_timing )
      elapsed_threads_candidates += stop_c_clock();



    //////////////////////////////////////////
    // allocate candidates memory on device //
    //////////////////////////////////////////
    // use: n_time_bins
    // cudamalloc: device_max_number_of_pmts_in_time_bin, device_vertex_with_max_n_pmts
    if( use_timing )
      start_cuda_clock();
    allocate_candidates_memory_on_device();
    if( use_timing )
      elapsed_candidates_memory_dev += stop_cuda_clock();



    /////////////////////////////////////
    // find candidates above threshold //
    /////////////////////////////////////
    if( use_timing )
      start_cuda_clock();
    if( use_verbose )
      printf(" --- execute candidates kernel \n");
    if( correct_mode != 9 ){
      kernel_find_vertex_with_max_npmts_in_timebin<<<number_of_kernel_blocks,number_of_threads_per_block>>>(device_n_pmts_per_time_bin, device_max_number_of_pmts_in_time_bin, device_vertex_with_max_n_pmts);
    }else{
      kernel_find_vertex_with_max_npmts_in_timebin_and_directionbin<<<number_of_kernel_blocks,number_of_threads_per_block>>>(device_n_pmts_per_time_bin_and_direction_bin, device_max_number_of_pmts_in_time_bin, device_vertex_with_max_n_pmts);
    }
    getLastCudaError("candidates_kernel execution failed\n");
    if( use_timing )
      elapsed_candidates_kernel += stop_cuda_clock();



    /////////////////////////////////////////
    // copy candidates from device to host //
    /////////////////////////////////////////
    // use: n_time_bins
    // memcpy: host_max_number_of_pmts_in_time_bin, host_vertex_with_max_n_pmts
    if( use_timing )
      start_cuda_clock();
    if( use_verbose )
      printf(" --- copy candidates from device to host \n");
    copy_candidates_from_device_to_host();
    if( use_timing )
      elapsed_candidates_copy_host += stop_cuda_clock();



    ///////////////////////////////////////
    // choose candidates above threshold //
    ///////////////////////////////////////
    if( use_timing )
      start_cuda_clock();
    if( use_verbose )
      printf(" --- choose candidates above threshold \n");
    choose_candidates_above_threshold();
    if( use_timing )
      choose_candidates = stop_cuda_clock();



    ///////////////////////
    // coalesce triggers //
    ///////////////////////
    if( use_timing )
      start_cuda_clock();
    coalesce_triggers();
    if( use_timing )
      elapsed_coalesce += stop_cuda_clock();




    //////////////////////////////////
    // separate triggers into gates //
    //////////////////////////////////
    if( use_timing )
      start_cuda_clock();
    separate_triggers_into_gates();
    if( use_timing )
      elapsed_gates += stop_cuda_clock();



    //////////////////
    // write output //
    //////////////////
    if( use_timing )
      start_cuda_clock();
    write_output();
    if( use_timing )
      elapsed_write_output += stop_cuda_clock();

    /////////////////////////////
    // deallocate event memory //
    /////////////////////////////
    if( use_timing )
      start_cuda_clock();
    if( use_verbose )
      printf(" --- deallocate memory \n");
    free_event_memories();
    if( use_timing )
      elapsed_free += stop_cuda_clock();

    n_events ++;

  }

  elapsed_total += stop_total_cuda_clock();


  printf(" ------ analyzed %d events \n", n_events);

  ///////////////////////
  // normalize timings //
  ///////////////////////
  if( use_timing ){
    elapsed_input /= n_events;
    elapsed_memory_candidates_host /= n_events;
    elapsed_threads /= n_events;
    elapsed_memory_dev /= n_events;
    elapsed_copy_dev /= n_events;
    elapsed_kernel_correct_times_and_get_n_pmts_per_time_bin /= n_events;
    elapsed_threads_candidates /= n_events;
    elapsed_candidates_memory_dev /= n_events;
    elapsed_candidates_kernel /= n_events;
    elapsed_candidates_copy_host /= n_events;
    elapsed_coalesce /= n_events;
    elapsed_gates /= n_events;
    elapsed_write_output /= n_events;
    elapsed_free /= n_events;
  }
  elapsed_total /= n_events;

  return 1;
}

int gpu_daq_finalize(){


  //////////////////////////////
  // deallocate global memory //
  //////////////////////////////
  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- deallocate tofs memory \n");
  free_global_memories();
  if( use_timing )
    elapsed_tofs_free = stop_cuda_clock();



  //////////////////
  // reset device //
  //////////////////
  // -- needed to flush the buffer which holds printf from each thread
  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- reset device \n");
  //  cudaDeviceReset();
  if( use_timing )
    elapsed_reset = stop_cuda_clock();



  //////////////////
  // print timing //
  //////////////////
  if( use_timing ){
    printf(" user parameters time : %f ms \n", elapsed_parameters);
    printf(" read pmts execution time : %f ms \n", elapsed_pmts);
    printf(" read detector execution time : %f ms \n", elapsed_detector);
    printf(" make test vertices execution time : %f ms \n", elapsed_vertices);
    printf(" setup threads candidates execution time : %f ms \n", elapsed_threads_candidates);
    printf(" make table of tofs execution time : %f ms \n", elapsed_tof);
    printf(" make table of directions execution time : %f ms \n", elapsed_directions);
    printf(" allocate tofs memory on device execution time : %f ms \n", elapsed_memory_tofs_dev);
    printf(" allocate directions memory on device execution time : %f ms \n", elapsed_memory_directions_dev);
    printf(" fill tofs memory on device execution time : %f ms \n", elapsed_tofs_copy_dev);
    printf(" fill directions memory on device execution time : %f ms \n", elapsed_directions_copy_dev);
    printf(" deallocate tofs memory execution time : %f ms \n", elapsed_tofs_free);
    printf(" device reset execution time : %f ms \n", elapsed_reset);
    printf(" read input execution time : %f ms (%f) \n", elapsed_input, elapsed_input/elapsed_total);
    printf(" allocate candidates memory on host execution time : %f ms (%f) \n", elapsed_memory_candidates_host, elapsed_memory_candidates_host/elapsed_total);
    printf(" setup threads execution time : %f ms (%f) \n", elapsed_threads, elapsed_threads/elapsed_total);
    printf(" allocate memory on device execution time : %f ms (%f) \n", elapsed_memory_dev, elapsed_memory_dev/elapsed_total);
    printf(" fill memory on device execution time : %f ms (%f) \n", elapsed_copy_dev, elapsed_copy_dev/elapsed_total);
    printf(" correct kernel execution time : %f ms (%f) \n", elapsed_kernel_correct_times_and_get_n_pmts_per_time_bin, elapsed_kernel_correct_times_and_get_n_pmts_per_time_bin/elapsed_total);
    printf(" allocate candidates memory on device execution time : %f ms (%f) \n", elapsed_candidates_memory_dev, elapsed_candidates_memory_dev/elapsed_total);
    printf(" copy candidates to host execution time : %f ms (%f) \n", elapsed_candidates_copy_host, elapsed_candidates_copy_host/elapsed_total);
    printf(" choose candidates execution time : %f ms (%f) \n", choose_candidates, choose_candidates/elapsed_total);
    printf(" candidates kernel execution time : %f ms (%f) \n", elapsed_candidates_kernel, elapsed_candidates_kernel/elapsed_total);
    printf(" coalesce triggers execution time : %f ms (%f) \n", elapsed_coalesce, elapsed_coalesce/elapsed_total);
    printf(" separate triggers into gates execution time : %f ms (%f) \n", elapsed_gates, elapsed_gates/elapsed_total);
    printf(" write output execution time : %f ms (%f) \n", elapsed_write_output, elapsed_write_output/elapsed_total);
    printf(" deallocate memory execution time : %f ms (%f) \n", elapsed_free, elapsed_free/elapsed_total);
  }
  printf(" total execution time : %f ms \n", elapsed_total);

  return 1;

}

