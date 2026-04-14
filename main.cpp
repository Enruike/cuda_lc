#include "globals.cuh"
#include<time.h>
//#include <bits/stdc++.h>
cudaError_t cudaStatus;

int total_points = 0;
double dt = 0.0;
double dtime = 0.0;
int cycle = 0;
bool flag;
double S, S2, S0;
double* d_Qold = nullptr;
double* d_Qnew = nullptr;
unsigned char* d_bulktype = nullptr;
signed char* d_signal = nullptr;
double* d_Qo = nullptr;
int* d_neighbor = nullptr;
unsigned int* d_Qtensor_index = nullptr;
unsigned int* d_Nvector_index = nullptr;
unsigned char* d_Qtensor_signal = nullptr;
unsigned char* d_Nvector_signal = nullptr;
double* d_nu = nullptr;

static inline unsigned int ceil_div(unsigned int n, unsigned int d) {
	return (n + d - 1) / d;
}

static bool report_cuda_malloc(void** ptr, size_t bytes, const char* name) {
	cudaStatus = cudaMalloc(ptr, bytes);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed for %s (%zu bytes, %.2f MiB): %s\n",
			name, bytes, bytes / 1024.0 / 1024.0, cudaGetErrorString(cudaStatus));
		return false;
	}
	return true;
}

static bool report_cuda_memcpy(void* dst, const void* src, size_t bytes, cudaMemcpyKind kind, const char* name) {
	cudaStatus = cudaMemcpy(dst, src, bytes, kind);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed for %s (%zu bytes, %.2f MiB): %s\n",
			name, bytes, bytes / 1024.0 / 1024.0, cudaGetErrorString(cudaStatus));
		return false;
	}
	return true;
}

static void print_cuda_environment() {
	int device = 0;
	int driver_version = 0;
	int runtime_version = 0;
	cudaDeviceProp prop;
	size_t free_mem = 0;
	size_t total_mem = 0;

	cudaStatus = cudaGetDevice(&device);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaGetDevice failed: %s\n", cudaGetErrorString(cudaStatus));
		return;
	}
	if (cudaGetDeviceProperties(&prop, device) != cudaSuccess) {
		fprintf(stderr, "cudaGetDeviceProperties failed.\n");
		return;
	}
	cudaDriverGetVersion(&driver_version);
	cudaRuntimeGetVersion(&runtime_version);
	cudaMemGetInfo(&free_mem, &total_mem);

	printf("\n~ CUDA Device ~\n");
	printf("Device %d: %s\n", device, prop.name);
	printf("Compute capability: %d.%d\n", prop.major, prop.minor);
	printf("Driver/runtime: %d / %d\n", driver_version, runtime_version);
	printf("Global memory: %.2f MiB total, %.2f MiB free\n\n",
		total_mem / 1024.0 / 1024.0, free_mem / 1024.0 / 1024.0);
}

__device__ double d_trqq(double Qin[6]){
        double ans = 0.;
        ans = Qin[0] * Qin[0] + Qin[3] * Qin[3] + Qin[5] * Qin[5]\
                + 2 * (Qin[1] * Qin[1] + Qin[2] * Qin[2] + Qin[4] * Qin[4]);
        return ans;
}

/* __global__ void test_symbol(void){
	printf("L1_dev is %lf", L1_dev);
}
 */
__global__ void d_checktrace(double* d_Qold, unsigned int droplet){

	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;

	if(indx < droplet){

		double tr = 0;
		double third =  1.0 / 3.0;
	
		tr = (d_Qold[indx * 6 + 0] + d_Qold[indx * 6 + 3] + d_Qold[indx * 6 + 5]) * third;

		if(tr > 1e-5) {
			
			printf("Correcting trace %lf for node %d!\n %f %f %f %f %f %f\n", tr, indx, d_Qold[indx * 6 + 0], d_Qold[indx * 6 + 1],\
				d_Qold[indx * 6 + 2], d_Qold[indx * 6 + 3], d_Qold[indx * 6 + 4], d_Qold[indx * 6 + 5]);

			d_Qold[indx * 6 + 0] -= tr;
			d_Qold[indx * 6 + 3] -= tr;
			d_Qold[indx * 6 + 5] -= tr;
			//printf("Non-tracelss.\n");			
		}

		double Qin[6] = { 0. };

		for(int i = 0; i < 6; i++){
			Qin[i] = d_Qold[indx * 6 + i];
		}

		if(d_trqq(Qin) > 1.){
	//              for(n = 0; n < 6; n ++){
	//                      Q[n] /= 1.3;
	//              }
			printf("Order parameter exceed 1. For node #%d\nQ info %f %f %f %f %f %f\nTrQQ:%lf\n", indx, d_Qold[indx * 6 + 0], d_Qold[indx * 6 + 1],\
				d_Qold[indx * 6 + 2], d_Qold[indx * 6 + 3], d_Qold[indx * 6 + 4], d_Qold[indx * 6 + 5], d_trqq(Qin));
		}
	}
}

	int main() {

	flag = true;
	double time_taken;
	time_t start, end;

	time(&start);

	//Lectura de parámetros.
	//Si los parámetros son verdaderos continuará.
	if (!read_param()) {
		printf("No file param.in found!\n");
		//EXIT_SUCCESS();
		return 0;
		exit(1);
	}

	//continua con la funcion initial ubicada en initialization.cpp
	else {
		old_en = 1.;
		//S and U are the values for the inner LC. Also irx, iry and irz are the corresponding radii.
		S = 0.25 * (1.0 + 3.0 * sqrt(1.0 - 8.0 / (3.0 * U)));
		//S2 and U2 are the values for outer LC. Rx, Ry, Rz would be the radii of the whole droplet.
		S2 = 0.25 * (1.0 + 3.0 * sqrt(1.0 - 8.0 / (3.0 * U2)));
		dt = tmin;
		dtime = (tmax - tmin) / increment;
		total_points = Nx * Ny * Nz;

		printf("Value for S1 is: %lf\n", S);
		printf("Value for S2 is: %lf\n", S2);

		if (!initial()) {
			printf("Geometry couldn't be initialized!\n");
			//EXIT_SUCCESS();
			return 0;
		}

			else {
				printf("Geometry successfully initialized!\n");
			}
			print_cuda_environment();
			//Freeing vectors used in geometry.
		

		unsigned int* h_Nvector_index;
		unsigned char* h_Nvector_signal;

		h_Nvector_index = (unsigned int*)malloc((surf + nsurf) * sizeof(unsigned int));
		h_Nvector_signal = (unsigned char*)malloc((surf + nsurf) * sizeof(unsigned char));

		unsigned int nb = 0;

		for (int i = 0; i < droplet; i++) {
			if ((signal[i] >= 2 && signal[i] <= 8) || (signal[i] == 12 || signal[i] == 13) || (signal[i] >= 20 && signal[i] <= 23)) {
				h_Nvector_index[nb] = i;				//We can find the Qtensor index of the Surface Vector point (nu_p or nu).
				h_Nvector_signal[nb] = signal[i];		//Type of point.
				nb++;
			}
		}

		if (nb != surf + nsurf) {
			printf("Error in transfer index and types for surface Qtensors!\n");
			printf("Count is %d, surf is %d & nsurf is %d!\n", nb, surf, nsurf);
			exit(1);
		}

		nb = 0;

		for (int i = 0; i < droplet; i++) {
			if ((signal[i] >= 2 && signal[i] <= 8) || (signal[i] == 12 || signal[i] == 13) || (signal[i] >= 20 && signal[i] <= 23)) {
				if (Qold[h_Nvector_index[nb]] != Qold[i]) {
					printf("Error in transfer from Qtensor to Surface Index Vector!\n");
					exit(1);
				}
				if (h_Nvector_signal[nb] != signal[i]) {
					printf("Error in transfer from Signal Vector to Surface Index Type Vector!\n");
					exit(1);
				}
				nb++;
			}
		}		

			unsigned int* h_Qtensor_index;
			unsigned char* h_Qtensor_signal;
			h_Qtensor_index = (unsigned int*)malloc(bulk * sizeof(unsigned int));
			h_Qtensor_signal = (unsigned char*)malloc(bulk * sizeof(unsigned char));

			nb = 0;

			for (int i = 0; i < droplet; i++) {
				if (signal[i] == 0 || signal[i] == 1) {
					h_Qtensor_index[nb] = i;
					h_Qtensor_signal[nb] = signal[i];
					nb++;
				}
			}

		if (nb != bulk) {
			printf("Error in transfer index and types for bulk Qtensors!\n");
			printf("Count is %d and bulk is %d!\n", nb, bulk);
			exit(1);
		}

			nb = 0;
			for(int i = 0; i < bulk;  i++) {
				if (signal[i] == 0 || signal[i] == 1) {
					if (Qold[h_Qtensor_index[nb]] != Qold[i]) {
						printf("Error in transfer from Qtensor to new Tensor Index Vector!\n");
						exit(1);
					}
					if (h_Qtensor_signal[nb] != signal[i]) {
						printf("Error in transfer from Signal Vector to new Tensor Index Type Vector!\n");
						exit(1);
					}
					nb++;
				}
			}
		nb = 0;

		//Remove signal vector in order to implement h_index_nu signal.

			// ********************************************** ----------------------- ****************************//
			//Allocating arrays in device memory. We add the error checker.
			if (!report_cuda_malloc((void**)&d_Qold, sizeof(double) * droplet * 6, "d_Qold")) {
				return 0;
			}
			if (!report_cuda_malloc((void**)&d_Qnew, sizeof(double) * droplet * 6, "d_Qnew")) {
				return 0;
			}
			if (!report_cuda_malloc((void**)&d_bulktype, sizeof(unsigned char) * droplet, "d_bulktype")) {
				return 0;
			}
			if (infinite == 0 && degenerate == 0){
				//d_Q0 allocation
				if (!report_cuda_malloc((void**)&d_Qo, sizeof(double) * (surf + nsurf) * 6, "d_Qo")) {
					return 0;
				}
			}
		// else{
		// 	cudaStatus = cudaMalloc((void**)&d_Qo, 1 * surf * 6);
		// 	if (cudaStatus != cudaSuccess) {
		// 		fprintf(stderr, "cudaMalloc failed!");
		// 		return 0;
		// 	}
		// }

			//*************we now change de sizes of new signal vectors.********************//

				if (!report_cuda_malloc((void**)&d_Qtensor_signal, sizeof(unsigned char) * bulk, "d_Qtensor_signal")) {
					return 0;
				}

				if (!report_cuda_malloc((void**)&d_Nvector_signal, sizeof(unsigned char) * (surf + nsurf), "d_Nvector_signal")) {
					return 0;
				}

			if (!report_cuda_malloc((void**)&d_Qtensor_index, sizeof(unsigned int) * bulk, "d_Qtensor_index")) {
				return 0;
			}

			if (!report_cuda_malloc((void**)&d_Nvector_index, sizeof(unsigned int) * (surf + nsurf), "d_Nvector_index")) {
				return 0;
			}
			//****************************************************************************//

			//Neighboor must have sign for -1 value.
			if (!report_cuda_malloc((void**)&d_neighbor, sizeof(int) * droplet * 6, "d_neighbor")) {
				return 0;
			}
			if (!report_cuda_malloc((void**)&d_nu, sizeof(double) * (surf + nsurf) * 3, "d_nu")) {
				return 0;
			}

		//We need h_bulktype to calculate energy in CPU. Don't free it.
			if (!report_cuda_memcpy(d_bulktype, h_bulktype, droplet * sizeof(unsigned char), cudaMemcpyHostToDevice, "d_bulktype")) {
				return 0;
			}
			if (infinite == 0 && degenerate == 0){
				if (!report_cuda_memcpy(d_Qo, Qo, (surf + nsurf) * 6 * sizeof(double), cudaMemcpyHostToDevice, "d_Qo")) {
					return 0;
				}
			}

			//Copy from host to device
			if (!report_cuda_memcpy(d_Qold, Qold, droplet * 6 * sizeof(double), cudaMemcpyHostToDevice, "d_Qold")) {
				return 0;
			}
			if (!report_cuda_memcpy(d_Qnew, Qold, droplet * 6 * sizeof(double), cudaMemcpyHostToDevice, "d_Qnew")) {
				return 0;
			}

		////****************************///
		//New vectors signal
			if (!report_cuda_memcpy(d_Nvector_signal, h_Nvector_signal, (surf + nsurf) * sizeof(unsigned char), cudaMemcpyHostToDevice, "d_Nvector_signal")) {
				return 0;
			}

				if (!report_cuda_memcpy(d_Qtensor_signal, h_Qtensor_signal, bulk * sizeof(unsigned char), cudaMemcpyHostToDevice, "d_Qtensor_signal")) {
					return 0;
				}

				//New vectors for index
				if (!report_cuda_memcpy(d_Nvector_index, h_Nvector_index, (surf + nsurf) * sizeof(unsigned int), cudaMemcpyHostToDevice, "d_Nvector_index")) {
					return 0;
				}
			if (!report_cuda_memcpy(d_Qtensor_index, h_Qtensor_index, bulk * sizeof(unsigned int), cudaMemcpyHostToDevice, "d_Qtensor_index")) {
				return 0;
			}

			//Free host-side compacted index buffers once they are on device.
			free(h_Qtensor_index);
			free(h_Qtensor_signal);
			free(h_Nvector_index);
			free(h_Nvector_signal);

		////****************************///

			if (!report_cuda_memcpy(d_neighbor, neighbor, droplet * 6 * sizeof(signed int), cudaMemcpyHostToDevice, "d_neighbor")) {
				return 0;
			}
			if (!report_cuda_memcpy(d_nu, nu, (surf + nsurf) * 3 * sizeof(double), cudaMemcpyHostToDevice, "d_nu")) {
				return 0;
			}

		/*printQ << <1, 10 >> > (d_Nvector_signal);
		cudaDeviceSynchronize();*/

			unsigned int threads_per_block = 256;
			unsigned int blk_thrds = 128;
		//size for surface
		unsigned int surfBlocks = ceil_div(surf + nsurf, threads_per_block);

		//size for bulk
		unsigned int bulkBlocks = ceil_div(bulk, blk_thrds);

		printf("The number of Bulk Blocks is %d\n", bulkBlocks);

		printf("The number of Surf Blocks is %d\n", surfBlocks);

		unsigned int dropletBlocks = ceil_div(droplet, threads_per_block);

		printf("The number of Droplet Blocks is %d\n\n", dropletBlocks);

		//__device__ double devThird;
		//double third = 1. / 3.;
		//cudaMemcpyToSymbol(devThird, &third, sizeof(double));
		//__host__ __device__ __constant__ double d_idx;
		//cudaMemcpyToSymbol(d_idx, &idx, sizeof(double));

		/* Copy symbols */
	/* 	cudaMemcpyToSymbol(chiral_dev, &chiral, sizeof(int));
		
		printf("host L1 value %lf", L1);

		cudaStatus = cudaMemcpyToSymbol(L1_dev, &L1, sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "failed Memcpy Symbol!");
			return 0;
		}

		test_symbol<<<1,16>>>();

		cudaMemcpyToSymbol(L2_dev, &L2, sizeof(float));
		cudaMemcpyToSymbol(L3_dev, &L3, sizeof(float));
		cudaMemcpyToSymbol(L4_dev, &L4, sizeof(float)); */

		//Progress bar
		const char *shade = "\u2592";
    	const char *shade2 = "\u2588";

		if(DoubleU && geo == 4){
			S0 = S2;
		}
		else{
			S0 = S;
		}
		
		if(stopat != 0){

			printf("Total Progress\n");
			printf("[");

			for(int i = 0; i < 50; i++){
				printf(" ");
			}
			
			printf("] 0.00%\n");
		}
	
		while (flag) {

			printf("\t\t ~Computing Energy~ \n");
			free_energy();
			
			if(fabs(dE) < accuracy || (stopat != 0 && cycle == stopat)){
				printf("Stopping condition reached; cycle: %d.\n", cycle);
				flag = false;
				break;
			}

			if(trace_checker > 0 && (cycle % trace_checker) == 0){
				printf("\t\t ~Checking Trace~ \n");
				printf("\033[1;31m");
				d_checktrace<<<dropletBlocks, threads_per_block>>>(d_Qold, droplet);
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "d_checktrace launch failed: %s\n", cudaGetErrorString(cudaStatus));
					return 0;
				}
				cudaDeviceSynchronize();
				printf("\033[0m");
			}

			if((cycle % save_every) == 0){
				printf("\x1b[32m");
				printf("\t\t ~Saving Data~ \n");
				printf("\033[0m");
				cudaStatus = cudaMemcpy(Qold, d_Qold, droplet * 6 * sizeof(double), cudaMemcpyDeviceToHost);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy failed while saving output!\n");
					return 0;
				}
				output();
			}

			printf("\t\t ~Relaxing~ \n");

			printf("\033[1;33m");
			//printf("\t");
			for(int i = 0; i < 50; i++){
				//std::cout << "\x2592";
				printf(shade);
			}

			printf("\r");
			//printf("\t");
			for (int i = 0; i < check_every; i++) {

				if((cycle%(check_every/50)==0)){
					printf(shade2);
				}

				/* test_symbol<<<1,32>>>();
				cudaDeviceSynchronize(); */

				if (!report_cuda_memcpy(d_Qnew, d_Qold, droplet * 6 * sizeof(double), cudaMemcpyDeviceToDevice, "d_Qnew <- d_Qold (bulk)")) {
					return 0;
				}

					if (L2 != 0.) {
						if (chiral == 1) {
							relax_bulk<true, true><<<bulkBlocks, blk_thrds>>>(d_Qold, d_Qnew, d_bulktype, d_neighbor, d_Qtensor_index, d_Qtensor_signal,
								U, U2, chiral, qch, L1, L2, bulk, idx, idy, idz, iddx, iddy, iddz, dt);
						}
						else {
							relax_bulk<true, false><<<bulkBlocks, blk_thrds>>>(d_Qold, d_Qnew, d_bulktype, d_neighbor, d_Qtensor_index, d_Qtensor_signal,
								U, U2, chiral, qch, L1, L2, bulk, idx, idy, idz, iddx, iddy, iddz, dt);
						}
					}
					else if (chiral == 1) {
						relax_bulk<false, true><<<bulkBlocks, blk_thrds>>>(d_Qold, d_Qnew, d_bulktype, d_neighbor, d_Qtensor_index, d_Qtensor_signal,
							U, U2, chiral, qch, L1, L2, bulk, idx, idy, idz, iddx, iddy, iddz, dt);
					}
					else {
						relax_bulk<false, false><<<bulkBlocks, blk_thrds>>>(d_Qold, d_Qnew, d_bulktype, d_neighbor, d_Qtensor_index, d_Qtensor_signal,
							U, U2, chiral, qch, L1, L2, bulk, idx, idy, idz, iddx, iddy, iddz, dt);
					}
				double* tmp = d_Qold;
				d_Qold = d_Qnew;
				d_Qnew = tmp;

				if (!report_cuda_memcpy(d_Qnew, d_Qold, droplet * 6 * sizeof(double), cudaMemcpyDeviceToDevice, "d_Qnew <- d_Qold (surf)")) {
					return 0;
				}
				if (chiral == 1) {
					relax_surf<true><<<surfBlocks, threads_per_block>>>(d_Qold, d_Qnew, d_neighbor, d_Nvector_index, d_Nvector_signal, d_Qo, chiral, qch, L1, L2, L3, L4,
						tiltAngle, (surf + nsurf), degenerate, infinite, W, Wp, d_nu, idx, idy, idz, dt, S0);
				}
				else {
					relax_surf<false><<<surfBlocks, threads_per_block>>>(d_Qold, d_Qnew, d_neighbor, d_Nvector_index, d_Nvector_signal, d_Qo, chiral, qch, L1, L2, L3, L4,
						tiltAngle, (surf + nsurf), degenerate, infinite, W, Wp, d_nu, idx, idy, idz, dt, S0);
				}
				tmp = d_Qold;
				d_Qold = d_Qnew;
				d_Qnew = tmp;
 
				if(dt < tmax){
					dt += dtime;
					if(dt >= tmax){
						dt = tmax;
					}
				}
				cycle++;
				//

			}
			printf("\033[0m\n");
			printf("\n");

			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Relaxation kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				return 0;
			}

			cudaStatus = cudaDeviceSynchronize();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Relaxation kernel execution failed: %s\n", cudaGetErrorString(cudaStatus));
				return 0;
			}

			if(stopat != 0){

				double percentage = (double)cycle / (double)stopat * 100.;

				printf("Total Progress\n");
				printf("[");
				for(int i = 0; i < rint(percentage / 2); i++){
					printf("#");
				}
				for(int i = 0; i < rint(100 / 2) - rint(percentage / 2); i++){
					printf(" ");
				}
				
				printf("] %.2lf%\n", percentage);
			}

			//system("clear");
			//printf("\n\t\t\t ~Done~ \n\n");

		}
	
		

		cudaStatus = cudaMemcpy(Qold, d_Qold, droplet * 6 * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed before final output!\n");
			return 0;
		}

		output();
		
		time(&end);

		time_taken = (double)(end - start);

		FILE* energy;

		energy = fopen("energy.out", "a");
		if(time_taken < 60){
			fprintf(energy, "\nTime used:	%lf min.\n", time_taken);
			printf("\nTime used:	%lf min.\n", time_taken);
		}
		else{
			fprintf(energy, "\nTime used:	%lf h.\n", time_taken / 60. / 60.);
			printf("\nTime used:	%lf h.\n", time_taken / 60. / 60.);
		}
		fclose(energy);	

	//free device variables
		cudaFree(d_Qold);
		cudaFree(d_Qnew);
		if(infinite == 0 && degenerate == 0){
			cudaFree(d_Qo);
		}
		
		cudaFree(d_bulktype);
			cudaFree(d_neighbor);
			cudaFree(d_nu);
			cudaFree(d_Qtensor_index);
			cudaFree(d_Qtensor_signal);
			cudaFree(d_Nvector_index);
			cudaFree(d_Nvector_signal);

		//free host variables
		free(signal);
		if(infinite == 0 && degenerate == 0){
			free(Qo);
		}
		free(nu);
		free(Qold);
		free(neighbor);
		free(h_bulktype);
		free(drop);
		free(boundary);
		return true;
	}

}
