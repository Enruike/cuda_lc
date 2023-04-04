#include "definitions.cuh"

//cuda variables. 
__device__ double fabs(double x);
__device__ __constant__ double devThird = 1. / 3.;
__device__ double delta[6] = { 1., 0., 0., 1., 0., 1. };
//__device__ __constant__ double d_idx, d_idy, d_idz, d_iddx, d_iddy, d_iddz;

__global__ void relax_bulk(double* d_Qold, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, int chiral,
	double U, double U2, double qch, int L1, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt) 
	{
	
	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;

	if (indx < bulk) {

		//double third = 1.0 / 3.0;
		double Qin[6];
		double QQ[6];
		double Qldg[6];
		double trQQ;
		int xm, xp, ym, yp, zm, zp;
		double dQ[3][6] = { {0} };
		double ddQ[6][6] = { {0} };
		double Qelas[6] = { 0 };
		double Qch[6] = { 0 };

		// if(indx < 10){
		// 	printf("idx %lf idy %lf idz %lf\n", idx, idy, idz);
		// }

		// if(indx < 10){
		// 	printf("iddx %lf iddy %lf iddz %lf\n", iddx, iddy, iddz);
		// }

		Qin[0] = d_Qold[d_Qtensor_index[indx] * 6 + 0];
		Qin[1] = d_Qold[d_Qtensor_index[indx] * 6 + 1];
		Qin[2] = d_Qold[d_Qtensor_index[indx] * 6 + 2];
		Qin[3] = d_Qold[d_Qtensor_index[indx] * 6 + 3];
		Qin[4] = d_Qold[d_Qtensor_index[indx] * 6 + 4];
		Qin[5] = d_Qold[d_Qtensor_index[indx] * 6 + 5];

		QQ[0] = Qin[0] * Qin[0] + Qin[1] * Qin[1] + Qin[2] * Qin[2];
		QQ[1] = Qin[0] * Qin[1] + Qin[1] * Qin[3] + Qin[2] * Qin[4];
		QQ[2] = Qin[0] * Qin[2] + Qin[1] * Qin[4] + Qin[2] * Qin[5];
		QQ[3] = Qin[1] * Qin[1] + Qin[3] * Qin[3] + Qin[4] * Qin[4];
		QQ[4] = Qin[1] * Qin[2] + Qin[3] * Qin[4] + Qin[4] * Qin[5];
		QQ[5] = Qin[2] * Qin[2] + Qin[4] * Qin[4] + Qin[5] * Qin[5];

		trQQ = Qin[0] * Qin[0] + Qin[3] * Qin[3] + Qin[5] * Qin[5]\
			+ 2 * (Qin[1] * Qin[1] + Qin[2] * Qin[2] + Qin[4] * Qin[4]);

		if (d_bulktype[d_Qtensor_index[indx]] == 1) {
			for (int i = 0; i < 6; i++) {
				Qldg[i] = (1 - U * devThird) * Qin[i] - U * (QQ[i] - trQQ * (Qin[i] + delta[i] * devThird));
			}
		}

		else if (d_bulktype[d_Qtensor_index[indx]] == 2) {
			for (int i = 0; i < 6; i++) {
				Qldg[i] = (1 - U2 * devThird) * Qin[i] - U2 * (QQ[i] - trQQ * (Qin[i] + delta[i] * devThird));
			}
		}

		xm = d_neighbor[d_Qtensor_index[indx] * 6 + 0];
		xp = d_neighbor[d_Qtensor_index[indx] * 6 + 1];
		ym = d_neighbor[d_Qtensor_index[indx] * 6 + 2];
		yp = d_neighbor[d_Qtensor_index[indx] * 6 + 3];
		zm = d_neighbor[d_Qtensor_index[indx] * 6 + 4];
		zp = d_neighbor[d_Qtensor_index[indx] * 6 + 5];

		//if(indx < 20){
		//	printf("I'm thread %d and my neighbors are %d %d %d %d %d %d\n", indx, xm, xp, ym, yp, zm, zp);
		//}

		for (int i = 0; i < 6; i++) {
			ddQ[0][i] = (d_Qold[xp * 6 + i] + d_Qold[xm * 6 + i] - 2. * Qin[i]) * iddx;
			ddQ[3][i] = (d_Qold[yp * 6 + i] + d_Qold[ym * 6 + i] - 2. * Qin[i]) * iddy;
			ddQ[5][i] = (d_Qold[zp * 6 + i] + d_Qold[zm * 6 + i] - 2. * Qin[i]) * iddz;
			Qelas[i] = ddQ[0][i] + ddQ[3][i] + ddQ[5][i];
		}

		if (chiral == 1) {
			for (int n = 0; n < 6; n++) {
				dQ[0][n] = (d_Qold[xp * 6 + n] - d_Qold[xm * 6 + n]) * 0.5 * idx;
				dQ[1][n] = (d_Qold[yp * 6 + n] - d_Qold[ym * 6 + n]) * 0.5 * idy;
				dQ[2][n] = (d_Qold[zp * 6 + n] - d_Qold[zm * 6 + n]) * 0.5 * idz;
			}
		}

		if (chiral == 1) {
			Qch[0] = 2. * (dQ[1][2] - dQ[2][1]);
			Qch[3] = 2. * (dQ[2][1] - dQ[0][4]);
			Qch[5] = 2. * (dQ[0][4] - dQ[1][2]);
			Qch[1] = dQ[1][4] - dQ[2][3] + dQ[2][0] - dQ[0][2];
			Qch[2] = dQ[1][5] - dQ[2][4] + dQ[0][1] - dQ[1][0];
			Qch[4] = dQ[2][2] - dQ[0][5] + dQ[0][3] - dQ[1][1];
		}

		__syncthreads();
		for (int n = 0; n < 6; n++) {
			//d_Qnew[indx * 6 + n] = Qin[n] + dt * (-Qldg[n] + L1 * Qelas[n] + (L2 + L4) * Qelas2[n] + L3 * Qelas3[n] - 2 * chiral * qch * L1 * Qch[n]);
			d_Qold[d_Qtensor_index[indx] * 6 + n] = Qin[n] + dt * (-Qldg[n] + L1 * Qelas[n] - 2 * chiral * qch * L1 * Qch[n]);
		}
		__syncthreads();
		// for (int i = 0; i < 6; i++) {
		//  	d_Qold[d_Qtensor_index[indx] * 6 + i] = d_Qnew[d_Qtensor_index[indx] * 6 + i];
		// }

		// __syncthreads();

	}

}

__global__ void relax_surf(double* d_Qold, signed int* d_neighbor, unsigned int* d_Nvector_index, unsigned char* d_Nvector_signal, double* d_Qo, 
	int chiral, double qch, int L1, unsigned int surf, int degenerate, int infinite, double W, double Wp, double* d_nu, double d_idx, double d_idy, double d_idz, double dt) {

	//__device__ double fabs(double x);
	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;

	if (indx < surf) {

		bool degen;
		bool inf;
		double Wstr;
		double loc_nu[3];
		double Qin[6];
		int xm, xp, ym, yp, zm, zp;
		double dQ[3][6] = { { 0. } };
		double Qelas[6] = { 0. };
		double Qch[6] = { 0 };

		// if(indx == 0){
		// 	printf("Surf is %d\n", surf);
		// }

		//for geo boundary
		
		if (d_Nvector_signal[indx] == 2 || d_Nvector_signal[indx] == 3) {
			degen = degenerate;
			inf = infinite;
			Wstr = W;
		}

		//for nanoparticle boundary
		//4 is for NInf
		//6 is for degenerate
		//8 is for Infinite homeotropic
		else if (d_Nvector_signal[indx] == 4 || d_Nvector_signal[indx] == 5) {
			degen = 0;
			inf = 0;
			Wstr = Wp;
		}
		else if (d_Nvector_signal[indx] == 6 || d_Nvector_signal[indx] == 7) {
			degen = 1;
			inf = 0;
			Wstr = Wp;
		}
		// else{
		// 	printf("Problems in surface node!\n");
		// 	return;
		// }

		if (inf == 0) {

			loc_nu[0] = d_nu[indx * 3 + 0];
			loc_nu[1] = d_nu[indx * 3 + 1];
			loc_nu[2] = d_nu[indx * 3 + 2];

			for (int i = 0; i < 6; i++) {
				Qin[i] = d_Qold[d_Nvector_index[indx] * 6 + i];
			}

			xm = d_neighbor[d_Nvector_index[indx] * 6 + 0];
			xp = d_neighbor[d_Nvector_index[indx] * 6 + 1];
			ym = d_neighbor[d_Nvector_index[indx] * 6 + 2];
			yp = d_neighbor[d_Nvector_index[indx] * 6 + 3];
			zm = d_neighbor[d_Nvector_index[indx] * 6 + 4];
			zp = d_neighbor[d_Nvector_index[indx] * 6 + 5];

			if ((d_Nvector_signal[indx] % 2) == 0) {
				for (int n = 0; n < 6; n++) {
					dQ[0][n] = (-d_Qold[xp * 6 + n] + 4 * d_Qold[xm * 6 + n] - 3 * Qin[n]) * 0.5 * d_idx;
					dQ[1][n] = (-d_Qold[yp * 6 + n] + 4 * d_Qold[ym * 6 + n] - 3 * Qin[n]) * 0.5 * d_idy;
					dQ[2][n] = (-d_Qold[zp * 6 + n] + 4 * d_Qold[zm * 6 + n] - 3 * Qin[n]) * 0.5 * d_idz;
					Qelas[n] = dQ[0][n] * fabs(loc_nu[0]) + dQ[1][n] * fabs(loc_nu[1]) + dQ[2][n] * fabs(loc_nu[2]);
				}
			}

			else if ((d_Nvector_signal[indx] % 2) == 1) {
				for (int n = 0; n < 6; n++) {
					if (xm == -1) {
						dQ[0][n] = 0;
					}
					else if (xp == -1) {
						dQ[0][n] = (d_Qold[xm * 6 + n] - Qin[n]) * d_idx;
					}
					else {
						dQ[0][n] = (-d_Qold[xp * 6 + n] + 4 * d_Qold[xm * 6 + n] - 3 * Qin[n]) * 0.5 * d_idx;
					}

					if (ym == -1) {
						dQ[1][n] = 0;
					}
					else if (yp  == -1) {
						dQ[1][n] = (d_Qold[ym * 6 + n] - Qin[n]) * d_idy;
					}
					else {
						dQ[1][n] = (-d_Qold[yp * 6 + n] + 4 * d_Qold[ym * 6 + n] - 3 * Qin[n]) * 0.5 * d_idy;
					}

					if (zm == -1) {
						dQ[2][n] = 0;
					}
					else if (zp == -1) {
						dQ[2][n] = (d_Qold[zm * 6 + n] - Qin[n]) * d_idz;
					}
					else {
						dQ[2][n] = (-d_Qold[zp * 6 + n] + 4 * d_Qold[zm * 6 + n] - 3 * Qin[n]) * 0.5 * d_idz;
					}

					Qelas[n] = dQ[0][n] * fabs(loc_nu[0]) + dQ[1][n] * fabs(loc_nu[1]) + dQ[2][n] * fabs(loc_nu[2]);
				}
			}

			if (chiral == 1) {
				Qch[0] = loc_nu[2] * Qin[1] - loc_nu[1] * Qin[2];
				Qch[3] = loc_nu[0] * Qin[4] - loc_nu[2] * Qin[1];
				Qch[5] = loc_nu[1] * Qin[2] - loc_nu[0] * Qin[4];
				Qch[1] = 0.5 * (loc_nu[2] * Qin[3] - loc_nu[1] * Qin[4] + loc_nu[0] * Qin[2] - loc_nu[2] * Qin[0]);
				Qch[2] = 0.5 * (loc_nu[2] * Qin[4] - loc_nu[1] * Qin[5] + loc_nu[1] * Qin[0] - loc_nu[0] * Qin[1]);
				Qch[4] = 0.5 * (loc_nu[0] * Qin[5] - loc_nu[2] * Qin[2] + loc_nu[1] * Qin[1] - loc_nu[0] * Qin[3]);
			}

		}

		__syncthreads();

		for (int n = 0; n < 6; n++) {
			d_Qold[d_Nvector_index[indx] * 6 + n] = Qin[n] + dt * (L1 * Qelas[n] + chiral * 2 * qch * Qch[n] - Wstr * (Qin[n] - d_Qo[indx * 6 + n]));
		}
		__syncthreads();
		
	}
}
