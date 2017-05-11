
INC	:= -I$(CUDA_HOME)/include -I.
LIB	:= -L$(CUDA_HOME)/lib64 -lcudart

NVCCFLAGS	:= -lineinfo -arch=sm_20 --ptxas-options=-v --use_fast_math
# NVCCFLAGS	:= -lineinfo -arch=sm_50 --ptxas-options=-v --use_fast_math

all:	daq_code_2d

daq_code_2d:	daq_code_2d.cu library_daq.h Makefile
	nvcc daq_code_2d.cu -o daq_code_2d $(INC) $(NVCCFLAGS) $(LIB)

clean:
	rm -f daq_code_2d
