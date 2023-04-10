#include "definitions.cuh"
#include<time.h>
//#include <bits/stdc++.h>
cudaError_t cudaStatus;

int total_points = 0;
double S = 0.0, S2 = 0.0;
double dE;
double dt = 0.0;
double dtime = 0.0;
int cycle = 0;
bool flag;
double en_tot;
double old_en;
double en_ldg;
double en_surf[2];
double en_el[5];

__device__ double d_trqq(double Qin[6]){
        double ans = 0.;
        ans = Qin[0] * Qin[0] + Qin[3] * Qin[3] + Qin[5] * Qin[5]\
                + 2 * (Qin[1] * Qin[1] + Qin[2] * Qin[2] + Qin[4] * Qin[4]);
        return ans;
}

__global__ void d_checktrace(double* d_Qold, unsigned int droplet){

	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;

	if(indx < droplet){

		double tr = 0;
		double third =  1.0 / 3.0;
	
		tr = (d_Qold[indx * 6 + 0] + d_Qold[indx * 6 + 3] + d_Qold[indx * 6 + 5]) * third;

		if(tr > 1e-5) {
			printf("\033[1;31m");
			printf("Correcting trace for node %d!\033[0m\n %f %f %f %f %f %f\n", indx, d_Qold[indx * 6 + 0], d_Qold[indx * 6 + 1],\
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
                printf("Order parameter exceed 1.\n");
            }
	}
}

int main() {

	//device variable callings
	double* d_Qold;
	unsigned char* d_Nvector_signal;
	//unsigned char* d_Qtensor_signal;
	unsigned int* d_Nvector_index;
	unsigned int* d_Qtensor_index;

	unsigned char* d_bulktype;
	double* d_Qnew;
	int* d_neighbor;
	double* d_nu;
	double* d_Qo;

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
	}

	//continua con la funcion initial ubicada en initialization.cpp
	else {


		dE = 1.0;
		//S and U are the values for the inner LC. Also irx, iry and irz are the corresponding radii.
		S = 0.25 * (1.0 + 3.0 * sqrt(1.0 - 8.0 / (3.0 * U)));
		//S2 and U2 are the values for outer LC. Rx, Ry, Rz would be the radii of the whole droplet.
		S2 = 0.25 * (1.0 + 3.0 * sqrt(1.0 - 8.0 / (3.0 * U2)));
		dt = tmin;
		dtime = (tmax - tmin) / increment;
		total_points = Nx * Ny * Nz;

		printf("Value for S1 is: %lf\n", S);
		printf("Value for S2 is: %lf\n\n", S2);

		if (!initial()) {
			printf("Geometry couldn't be initialized!\n");
			//EXIT_SUCCESS();
			return 0;
		}

		else {
			printf("Geometry successfully initialized!\n");
		}
		//Freeing vectors used in geometry.
		

		unsigned int* h_Nvector_index;
		unsigned char* h_Nvector_signal;

		h_Nvector_index = (unsigned int*)malloc(surf * sizeof(unsigned int));
		h_Nvector_signal = (unsigned char*)malloc(surf * sizeof(unsigned char));

		unsigned int nb = 0;

		for (int i = 0; i < droplet; i++) {
			if ((signal[i] >= 2 && signal[i] <= 8) || (signal[i] == 12 || signal[i] == 13)) {
				h_Nvector_index[nb] = i;						//We can find the Qtensor index of the Surface Vector point (nu_p or nu).
				h_Nvector_signal[nb] = signal[i];		//Type of point.
				nb++;
			}
		}

		if (nb != surf) {
			printf("Error in transfer index and types for surface Qtensors!\n");
			exit(1);
		}

		nb = 0;

		for (int i = 0; i < droplet; i++) {
			if ((signal[i] >= 2 && signal[i] <= 8) || (signal[i] == 12 || signal[i] == 13)) {
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
		//unsigned char* h_Qtensor_signal;

		unsigned int bulk = droplet - surf;

		h_Qtensor_index = (unsigned int*)malloc(bulk * sizeof(unsigned int));
		//h_Qtensor_signal = (unsigned char*)malloc(bulk * sizeof(unsigned char));

		nb = 0;

		for (int i = 0; i < droplet; i++) {
			if (signal[i] == 0 || signal[i] == 1) {
				h_Qtensor_index[nb] = i;
				//h_Qtensor_signal[nb] = signal[i];
				nb++;
			}
		}

		if (nb != bulk) {
			printf("Error in transfer index and types for bulk Qtensors!\n");
			exit(1);
		}

		nb = 0;
		for(int i = 0; i < bulk;  i++) {
			if (signal[i] == 0 || signal[i] == 1) {
				if (Qold[h_Qtensor_index[nb]] != Qold[i]) {
					printf("Error in transfer from Qtensor to new Tensor Index Vector!\n");
					exit(1);
				}
				/*if (h_Qtensor_signal[nb] != signal[i]) {
					printf("Error in transfer from Signal Vector to new Tensor Index Type Vector!\n");
					exit(1);
				}*/
				nb++;
			}
		}
		nb = 0;

		//Remove signal vector in order to implement h_index_nu signal.

		// ********************************************** ----------------------- ****************************//
		//Allocating arrays in device memory. We add the error checker.
		cudaStatus = cudaMalloc((void**)&d_Qold, sizeof(double) * droplet * 6);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}
		cudaStatus = cudaMalloc((void**)&d_bulktype, sizeof(unsigned char) * droplet);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}
		//d_Q0 allocation
		cudaStatus = cudaMalloc((void**)&d_Qo, sizeof(double) * surf * 6);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}

		//*************we now change de sizes of new signal vectors.********************//

	/*	cudaStatus = cudaMalloc((unsigned char**)&d_Qtensor_signal, sizeof(unsigned char) * bulk);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}*/

		cudaStatus = cudaMalloc((void**)&d_Nvector_signal, sizeof(unsigned char) * surf);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}

		cudaStatus = cudaMalloc((void**)&d_Qtensor_index, sizeof(unsigned int) * bulk);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}

		cudaStatus = cudaMalloc((void**)&d_Nvector_index, sizeof(unsigned int) * surf);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}
		//****************************************************************************//

		//Neighboor must have sign for -1 value.
		cudaStatus = cudaMalloc((void**)&d_neighbor, sizeof(int) * droplet * 6);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}
		cudaStatus = cudaMalloc((void**)&d_nu, sizeof(double) * surf * 3);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return 0;
		}

		//We need h_bulktype to calculate energy in CPU. Don't free it.
		cudaStatus = cudaMemcpy(d_bulktype, h_bulktype, droplet * sizeof(unsigned char), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}

		cudaStatus = cudaMemcpy(d_Qo, Qo, surf * 6 * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}

		//Copy from host to device
		cudaStatus = cudaMemcpy(d_Qold, Qold, droplet * 6 * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}

		////****************************///
		//New vectors signal
		cudaStatus = cudaMemcpy(d_Nvector_signal, h_Nvector_signal, surf * sizeof(unsigned char), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}

		//New vectors for index
		cudaStatus = cudaMemcpy(d_Nvector_index, h_Nvector_index, surf * sizeof(unsigned int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}
		cudaStatus = cudaMemcpy(d_Qtensor_index, h_Qtensor_index, bulk * sizeof(unsigned int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}

		//Freeing h_Qtensor_index, h_Nvector_index, h_Nvector_signal
		cudaFree(h_Qtensor_index);
		cudaFree(h_Nvector_index);
		cudaFree(h_Nvector_signal);

		////****************************///

		cudaStatus = cudaMemcpy(d_neighbor, neighbor, droplet * 6 * sizeof(signed int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}
		cudaStatus = cudaMemcpy(d_nu, nu, surf * 3 * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			return 0;
		}

		/*printQ << <1, 10 >> > (d_Nvector_signal);
		cudaDeviceSynchronize();*/


		unsigned int threads_per_block = 512;
		//size for surface
		unsigned int surfBlocks = rint(surf / threads_per_block) + 1;		

		//size for bulk
		unsigned int bulkBlocks = rint(bulk / threads_per_block) + 1;

		unsigned int dropletBlocks = rint(droplet / threads_per_block) + 1;

		printf("The number of Bulk Blocks is %d\n", bulkBlocks);

		printf("The number of Surf Blocks is %d\n", surfBlocks);

		printf("The number of Droplet Blocks is %d\n\n", dropletBlocks);

		//__device__ double devThird;
		//double third = 1. / 3.;
		//cudaMemcpyToSymbol(devThird, &third, sizeof(double));
		//__host__ __device__ __constant__ double d_idx;
		//cudaMemcpyToSymbol(d_idx, &idx, sizeof(double));

		//Progress bar
		const char *shade = "\u2592";
    	const char *shade2 = "\u2588";

		while (flag) {

			printf("\t\t ~Computing Energy~ \n");
 			free_energy();
			
			if(fabs(dE) < accuracy || (stopat != 0 && cycle == stopat)){
				printf("Error in the trace of q; cycle : %d.\n", cycle);
				flag = false;
				break;
			}

			printf("\t\t ~Checking Trace~ \n");
			d_checktrace<<<dropletBlocks, threads_per_block>>>(d_Qold, droplet);
			cudaDeviceSynchronize();

			//if((cycle % check_every) == 0){

				
				// for(int i = 0; i < droplet; i++){
				// 	//				checktr(&q[i * 6]);
				// 	if(!checktr(&Qold[i * 6])){
				// 		//flag = false;
				// 		printf("Trace corrected!\n");
				// 		printf("%d\n", i);

				// 	}
				// }				
				//if(!flag){
				//	printf("Error in the trace of q; cycle : %d.\n", cycle);
				//}
			//}

			if((cycle % save_every) == 0){
				printf("\t\t ~Saving Data~ \n");
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

				relax_bulk<<<bulkBlocks, threads_per_block>>>(d_Qold, d_bulktype, d_neighbor, d_Qtensor_index, chiral, U, U2, qch, L1, bulk, idx, idy, idz,
					iddx, iddy, iddz, dt);
				cudaDeviceSynchronize();

				relax_surf<<<surfBlocks, threads_per_block>>>(d_Qold, d_neighbor, d_Nvector_index, d_Nvector_signal, d_Qo, chiral, qch, L1, surf, degenerate,
					infinite, W, Wp, d_nu, idx, idy, idz, dt);

				cudaDeviceSynchronize();

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

			if(stopat != 0){

				double percentage = (double)cycle / (double)stopat * 100.;

				printf("Total Progress\n");
				printf("[");
				for(int i = 0; i < rint(percentage / 2); i++){
					printf("#");
				}
				for(int i = 0; i < rint(100 / 2) - 1 - percentage; i++){
					printf(" ");
				}
				
				printf("] %.2lf%\n", percentage);
			}

			//system("clear");
			//printf("\n\t\t\t ~Done~ \n\n");

			// flag=false;
			//printf("      <<=== ~Copying Q-Tensor back to Host Memory~ ===>>>\n");
			cudaStatus = cudaMemcpy(Qold, d_Qold, droplet * 6 * sizeof(double), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy failed!\n");
				return 0;
				break;
			}
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
		cudaFree(d_Qo);
		cudaFree(d_bulktype);
		cudaFree(d_neighbor);
		cudaFree(d_nu);
		cudaFree(d_Qtensor_index);
		cudaFree(d_Nvector_index);
		cudaFree(d_Nvector_signal);

		//free host variables
		free(signal);
		free(Qo);
		free(nu);
		free(Qold);
		free(neighbor);
		free(h_bulktype);
		free(drop);
		free(boundary);
		return true;
	}

}
