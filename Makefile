# Auto-detect GPU compute capability via nvidia-smi, fallback to Ada (89).
DETECTED_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | grep -Eo '[0-9]+\.[0-9]+' | head -n 1 | tr -d '.')
ARCH ?= $(if $(DETECTED_ARCH),$(DETECTED_ARCH),89)
FMAD ?= 1
DEBUG_ENERGY_COMPARE ?= 0
DEBUG_CHANNEL_SURF_TRACE ?= 0

COMMON_NVCC_FLAGS = -O3 -w -DDEBUG_ENERGY_COMPARE=$(DEBUG_ENERGY_COMPARE) -DDEBUG_CHANNEL_SURF_TRACE=$(DEBUG_CHANNEL_SURF_TRACE)
FMAD_FLAG = --fmad=$(if $(filter 1,$(FMAD)),true,false)
GPU_ARCH_FLAGS = -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING)

OBJS = configurations.o energy.o functions.o initial.o read_parameters.o solver_bulk.o solver_surface.o main.o geometry.o
HEADERS = globals.cuh solver_kernels.cuh solver_common.cuh geometry.hpp initial.hpp read_parameters.hpp
NO_ARCH_WARNING = -Wno-deprecated-gpu-targets

all: lc_cuda.x

lc_cuda.x: $(OBJS)
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) $(GPU_ARCH_FLAGS) -o lc_cuda.x $(OBJS)

configurations.o: configurations.cu $(HEADERS)
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) $(GPU_ARCH_FLAGS) -c configurations.cu

#ellipsoid.o: ellipsoid.cpp geometry.hpp
#	nvcc -c ellipsoid.cpp

energy.o: energy.cpp energy.hpp
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) -x cu $(GPU_ARCH_FLAGS) -c energy.cpp

#nanochannel.o: nanochannel.cpp geometry.hpp
#	nvcc -c nanochannel.cpp

geometry.o: geometry.cpp geometry.hpp
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) -c geometry.cpp

functions.o: functions.cpp functions.hpp
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) -c functions.cpp

initial.o: initial.cpp initial.hpp $(HEADERS)
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) -c initial.cpp

read_parameters.o: read_parameters.cpp read_parameters.hpp $(HEADERS)
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) -c read_parameters.cpp

solver_bulk.o: solver_bulk.cu globals.cuh solver_kernels.cuh solver_common.cuh
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) $(GPU_ARCH_FLAGS) -c solver_bulk.cu

solver_surface.o: solver_surface.cu globals.cuh solver_kernels.cuh solver_common.cuh
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) $(GPU_ARCH_FLAGS) -c solver_surface.cu

main.o: main.cpp globals.cuh solver_kernels.cuh
	nvcc $(COMMON_NVCC_FLAGS) $(FMAD_FLAG) -x cu $(GPU_ARCH_FLAGS) -c main.cpp

clean:
	rm -f lc_cuda.x $(OBJS)

#nvcc *.cpp *.cu -O2 -w --fmad=false -gencode=arch=compute_37,code=sm_37 -o lc_cuda.x
