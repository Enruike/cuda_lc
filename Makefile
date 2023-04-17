#37 for K80 GPUs and 75 for Turing architecture GPUs
ARCH = 75
OBJS = configurations.o ellipsoid.o energy.o functions.o initial.o read_parameters.o relaxations.o kernel.o
FILES = configurations.cu ellipsoid.cpp energy.cpp functions.cu initial.cpp kernel.cu read_parameters.cpp relaxations.cu
HEADERS = definitions.cuh geometry.hpp initial.hpp read_parameters.hpp
NO_ARCH_WARNING = -Wno-deprecated-gpu-targets

all: lc_cuda.x

lc_cuda.x: $(FILES) $(HEADERS) #$(OBJS)
	nvcc -O2 -w -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -o lc_cuda.x $(FILES)

#configurations.o: configurations.cu definitions.cuh
#	nvcc -c configurations.cu

#ellipsoid.o: ellipsoid.cpp geometry.hpp
#	nvcc -c ellipsoid.cpp

# energy.o: energy.cpp definitions.cuh
# 	nvcc -c energy.cpp

# functions.o: functions.cu definitions.cuh
# 	nvcc -c functions.cu

# initial.o: initial.cpp initial.hpp
# 	nvcc -c initial.cpp

# read_parameters.o: read_parameters.cpp read_parameters.hpp
# 	nvcc -c read_parameters.cpp

# relaxations.o: relaxations.cu definitions.cuh
# 	nvcc -c relaxations.cu

# kernel.o: kernel.cu definitions.cuh
# 	nvcc -c kernel.cu

clean:
	rm -f lc_cuda.x $(OBJS)

#nvcc *.cpp *.cu -O2 -w --fmad=false -gencode=arch=compute_37,code=sm_37 -o lc_cuda.x