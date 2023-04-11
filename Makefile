#37 for K80 GPUs and 75 for Turing architecture GPUs
ARCH = 37 
OBJS = configurations ellipsoid energy functions initial read_parameters relaxations kernel

all: lc_cuda

lc_cuda: $(OBJS).o
	nvcc $(OBJS).o -O3 -w --fmad=false -gencode=arch=compute_$(ARCH),code=sm_$(ARCH)

configurations.o: configurations.cpp definitions.cuh
	nvcc -c configurations.cpp

ellipsoid.o: ellipsoid.cpp geometry.hpp
	nvcc -c ellipsoid.cpp

energy.o: energy.cpp definitions.cuh
	nvcc -c energy.cpp

functions.o: functions.cu definitions.cuh
	nvcc -c functions.cu

initial.o: initial.cpp initial.hpp
	nvcc -c initial.cpp

read_parameters.o: read_parameters.cpp read_parameters.hpp
	nvcc -c read_parameters.cpp

relaxations.o: relaxations.cu definitions.cuh
	nvcc -c relaxation.cu

kernel.o: kernel.cu definitions.cuh
	nvcc -c kernel.cu

clean:
	rm -f lc_cuda $(OBJS).o

#nvcc *.cpp *.cu -O2 -w --fmad=false -gencode=arch=compute_37,code=sm_37 -o lc_cuda.x