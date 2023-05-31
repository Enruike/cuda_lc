#37 for K80 GPUs and 75 for Turing architecture GPUs
ARCH = 75
OBJS = configurations.o energy.o functions.o initial.o read_parameters.o relaxations.o kernel.o geometry.o
OBJ = energy.o geometry.o functions.o
FILES = configurations.cu initial.cpp kernel.cu read_parameters.cpp relaxations.cu
HEADERS = definitions.cuh geometry.hpp initial.hpp read_parameters.hpp
NO_ARCH_WARNING = -Wno-deprecated-gpu-targets

all: lc_cuda.x

lc_cuda.x: $(FILES) $(HEADERS) $(OBJ)
	nvcc -O3 -w --fmad=false -gencode=arch=compute_$(ARCH),code=sm_$(ARCH) $(NO_ARCH_WARNING) -o lc_cuda.x $(FILES) $(OBJ)

#configurations.o: configurations.cu definitions.cuh
#	nvcc -c configurations.cu

#ellipsoid.o: ellipsoid.cpp geometry.hpp
#	nvcc -c ellipsoid.cpp

energy.o: energy.cpp energy.hpp
	nvcc -c energy.cpp

#nanochannel.o: nanochannel.cpp geometry.hpp
#	nvcc -c nanochannel.cpp

geometry.o: geometry.cpp geometry.hpp
	nvcc -c geometry.cpp

functions.o: functions.cpp functions.hpp
	nvcc -c functions.cpp

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
