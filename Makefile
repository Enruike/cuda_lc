# Auto-detect GPU compute capability via nvidia-smi, fallback to Ada (89).
DETECTED_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | grep -Eo '[0-9]+\.[0-9]+' | head -n 1 | tr -d '.')
ARCH ?= $(if $(DETECTED_ARCH),$(DETECTED_ARCH),89)
OBJS = configurations.o energy.o functions.o initial.o read_parameters.o solver_bulk.o solver_surface.o main.o geometry.o
HEADERS = globals.cuh solver_kernels.cuh solver_common.cuh geometry.hpp initial.hpp read_parameters.hpp
NO_ARCH_WARNING = -Wno-deprecated-gpu-targets

all: lc_cuda.x

lc_cuda.x: $(OBJS)
	nvcc -O3 -w --fmad=true -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -o lc_cuda.x $(OBJS)

configurations.o: configurations.cu $(HEADERS)
	nvcc -O3 -w --fmad=true -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -c configurations.cu

#ellipsoid.o: ellipsoid.cpp geometry.hpp
#	nvcc -c ellipsoid.cpp

energy.o: energy.cpp energy.hpp
	nvcc -O3 -w --fmad=true -x cu -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -c energy.cpp

#nanochannel.o: nanochannel.cpp geometry.hpp
#	nvcc -c nanochannel.cpp

geometry.o: geometry.cpp geometry.hpp
	nvcc -O3 -w --fmad=true -c geometry.cpp

functions.o: functions.cpp functions.hpp
	nvcc -O3 -w --fmad=true -c functions.cpp

initial.o: initial.cpp initial.hpp $(HEADERS)
	nvcc -O3 -w --fmad=true -c initial.cpp

read_parameters.o: read_parameters.cpp read_parameters.hpp $(HEADERS)
	nvcc -O3 -w --fmad=true -c read_parameters.cpp

solver_bulk.o: solver_bulk.cu globals.cuh solver_kernels.cuh solver_common.cuh
	nvcc -O3 -w --fmad=true -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -c solver_bulk.cu

solver_surface.o: solver_surface.cu globals.cuh solver_kernels.cuh solver_common.cuh
	nvcc -O3 -w --fmad=true -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -c solver_surface.cu

main.o: main.cpp globals.cuh solver_kernels.cuh
	nvcc -O3 -w --fmad=true -x cu -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -c main.cpp

clean:
	rm -f lc_cuda.x $(OBJS)

#nvcc *.cpp *.cu -O2 -w --fmad=false -gencode=arch=compute_37,code=sm_37 -o lc_cuda.x
