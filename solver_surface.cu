#include "solver_kernels.cuh"
#include "solver_common.cuh"

//variable degenerate it's only relevant for the system surface. For nanoparticle surface, we will use degen and signal types.
template <bool UseChiral>
__global__ void relax_surf(const double* d_Qin_old, double* d_Qout, signed int* d_neighbor, unsigned int* d_Nvector_index, unsigned char* d_Nvector_signal, double* d_Qo, 
	int chiral, double qch, double L1, double L2, double L3, double L4, double tiltAngle, unsigned int surf, int degenerate, int infinite, double Wstr, double Wp, double* d_nu, double d_idx, 
	double d_idy, double d_idz, double dt, double S){

	//__device__ double fabs(double x);
	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;
	
	if (indx < surf) {

		unsigned int signal = d_Nvector_signal[indx];
		unsigned int nv_indx = d_Nvector_index[indx];

		//8 is for Infinite homeotropic
		if(signal == 8){
			return;
		}

		double loc_nu[3] = { 0. };
		double Qin[6] = { 0. };
		double Qdiff[6] = { 0. };
		int xm, xp, ym, yp, zm, zp;
		double dQ[3][6] = { { 0. } };
		double Qelas[6] = { 0. };
		double Qelas2[6] = { 0. };
		double Qch[6] = { 0. };

		

		// if(indx == 0){
		// 	printf("Surf is %d\n", surf);
		// }

		//for geo boundary
		
		//for nanoparticle boundary
		//4 is for NInf
		//6 is for degenerate
		
		if (signal == 4 || signal == 5) {
			degenerate = 0;
			infinite = 0;
			Wstr = Wp;
		}
		else if (signal == 6 || signal == 7) {
			degenerate = 1;
			infinite = 0;
			Wstr = Wp;
		}
		else if(signal == 20 || signal == 21){
			degenerate = 1;
			infinite = 0;
			Wstr = Wp;
		}
		else if(signal == 22 || signal == 23){
			degenerate = 2;
			infinite = 0;
			Wstr = Wp;
		}
		// else{
		// 	printf("Problems in surface node!\n");
		// 	return;
		// }

		if (infinite == 0) {

			loc_nu[0] = d_nu[indx * 3 + 0];
			loc_nu[1] = d_nu[indx * 3 + 1];
			loc_nu[2] = d_nu[indx * 3 + 2];
			Qin[0] = d_Qin_old[nv_indx * 6 + 0];
			Qin[1] = d_Qin_old[nv_indx * 6 + 1];
			Qin[2] = d_Qin_old[nv_indx * 6 + 2];
			Qin[3] = d_Qin_old[nv_indx * 6 + 3];
			Qin[4] = d_Qin_old[nv_indx * 6 + 4];
			Qin[5] = d_Qin_old[nv_indx * 6 + 5];

			xm = d_neighbor[nv_indx * 6 + 0];
			xp = d_neighbor[nv_indx * 6 + 1];
			ym = d_neighbor[nv_indx * 6 + 2];
			yp = d_neighbor[nv_indx * 6 + 3];
			zm = d_neighbor[nv_indx * 6 + 4];
			zp = d_neighbor[nv_indx * 6 + 5];

			if ((signal % 2) == 0) {
				dQ[0][0] = (-d_Qin_old[xp * 6 + 0] + 4 * d_Qin_old[xm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idx;
				dQ[1][0] = (-d_Qin_old[yp * 6 + 0] + 4 * d_Qin_old[ym * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idy;
				dQ[2][0] = (-d_Qin_old[zp * 6 + 0] + 4 * d_Qin_old[zm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idz;
				Qelas[0] = dQ[0][0] * fabs(loc_nu[0]) + dQ[1][0] * fabs(loc_nu[1]) + dQ[2][0] * fabs(loc_nu[2]);

				dQ[0][1] = (-d_Qin_old[xp * 6 + 1] + 4 * d_Qin_old[xm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idx;
				dQ[1][1] = (-d_Qin_old[yp * 6 + 1] + 4 * d_Qin_old[ym * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idy;
				dQ[2][1] = (-d_Qin_old[zp * 6 + 1] + 4 * d_Qin_old[zm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idz;
				Qelas[1] = dQ[0][1] * fabs(loc_nu[0]) + dQ[1][1] * fabs(loc_nu[1]) + dQ[2][1] * fabs(loc_nu[2]);

				dQ[0][2] = (-d_Qin_old[xp * 6 + 2] + 4 * d_Qin_old[xm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idx;
				dQ[1][2] = (-d_Qin_old[yp * 6 + 2] + 4 * d_Qin_old[ym * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idy;
				dQ[2][2] = (-d_Qin_old[zp * 6 + 2] + 4 * d_Qin_old[zm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idz;
				Qelas[2] = dQ[0][2] * fabs(loc_nu[0]) + dQ[1][2] * fabs(loc_nu[1]) + dQ[2][2] * fabs(loc_nu[2]);

				dQ[0][3] = (-d_Qin_old[xp * 6 + 3] + 4 * d_Qin_old[xm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idx;
				dQ[1][3] = (-d_Qin_old[yp * 6 + 3] + 4 * d_Qin_old[ym * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idy;
				dQ[2][3] = (-d_Qin_old[zp * 6 + 3] + 4 * d_Qin_old[zm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idz;
				Qelas[3] = dQ[0][3] * fabs(loc_nu[0]) + dQ[1][3] * fabs(loc_nu[1]) + dQ[2][3] * fabs(loc_nu[2]);

				dQ[0][4] = (-d_Qin_old[xp * 6 + 4] + 4 * d_Qin_old[xm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idx;
				dQ[1][4] = (-d_Qin_old[yp * 6 + 4] + 4 * d_Qin_old[ym * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idy;
				dQ[2][4] = (-d_Qin_old[zp * 6 + 4] + 4 * d_Qin_old[zm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idz;
				Qelas[4] = dQ[0][4] * fabs(loc_nu[0]) + dQ[1][4] * fabs(loc_nu[1]) + dQ[2][4] * fabs(loc_nu[2]);

				dQ[0][5] = (-d_Qin_old[xp * 6 + 5] + 4 * d_Qin_old[xm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idx;
				dQ[1][5] = (-d_Qin_old[yp * 6 + 5] + 4 * d_Qin_old[ym * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idy;
				dQ[2][5] = (-d_Qin_old[zp * 6 + 5] + 4 * d_Qin_old[zm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idz;
				Qelas[5] = dQ[0][5] * fabs(loc_nu[0]) + dQ[1][5] * fabs(loc_nu[1]) + dQ[2][5] * fabs(loc_nu[2]);
			}
			else if ((signal % 2) == 1) {
				if (xm == -1) {
					dQ[0][0] = 0;
					dQ[0][1] = 0;
					dQ[0][2] = 0;
					dQ[0][3] = 0;
					dQ[0][4] = 0;
					dQ[0][5] = 0;
				}
				else if (xp == -1) {
					dQ[0][0] = (d_Qin_old[xm * 6 + 0] - Qin[0]) * d_idx;
					dQ[0][1] = (d_Qin_old[xm * 6 + 1] - Qin[1]) * d_idx;
					dQ[0][2] = (d_Qin_old[xm * 6 + 2] - Qin[2]) * d_idx;
					dQ[0][3] = (d_Qin_old[xm * 6 + 3] - Qin[3]) * d_idx;
					dQ[0][4] = (d_Qin_old[xm * 6 + 4] - Qin[4]) * d_idx;
					dQ[0][5] = (d_Qin_old[xm * 6 + 5] - Qin[5]) * d_idx;
				}
				else {
					dQ[0][0] = (-d_Qin_old[xp * 6 + 0] + 4 * d_Qin_old[xm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idx;
					dQ[0][1] = (-d_Qin_old[xp * 6 + 1] + 4 * d_Qin_old[xm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idx;
					dQ[0][2] = (-d_Qin_old[xp * 6 + 2] + 4 * d_Qin_old[xm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idx;
					dQ[0][3] = (-d_Qin_old[xp * 6 + 3] + 4 * d_Qin_old[xm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idx;
					dQ[0][4] = (-d_Qin_old[xp * 6 + 4] + 4 * d_Qin_old[xm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idx;
					dQ[0][5] = (-d_Qin_old[xp * 6 + 5] + 4 * d_Qin_old[xm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idx;
				}

				if (ym == -1) {
					dQ[1][0] = 0;
					dQ[1][1] = 0;
					dQ[1][2] = 0;
					dQ[1][3] = 0;
					dQ[1][4] = 0;
					dQ[1][5] = 0;
				}
				else if (yp == -1) {
					dQ[1][0] = (d_Qin_old[ym * 6 + 0] - Qin[0]) * d_idy;
					dQ[1][1] = (d_Qin_old[ym * 6 + 1] - Qin[1]) * d_idy;
					dQ[1][2] = (d_Qin_old[ym * 6 + 2] - Qin[2]) * d_idy;
					dQ[1][3] = (d_Qin_old[ym * 6 + 3] - Qin[3]) * d_idy;
					dQ[1][4] = (d_Qin_old[ym * 6 + 4] - Qin[4]) * d_idy;
					dQ[1][5] = (d_Qin_old[ym * 6 + 5] - Qin[5]) * d_idy;
				}
				else {
					dQ[1][0] = (-d_Qin_old[yp * 6 + 0] + 4 * d_Qin_old[ym * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idy;
					dQ[1][1] = (-d_Qin_old[yp * 6 + 1] + 4 * d_Qin_old[ym * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idy;
					dQ[1][2] = (-d_Qin_old[yp * 6 + 2] + 4 * d_Qin_old[ym * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idy;
					dQ[1][3] = (-d_Qin_old[yp * 6 + 3] + 4 * d_Qin_old[ym * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idy;
					dQ[1][4] = (-d_Qin_old[yp * 6 + 4] + 4 * d_Qin_old[ym * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idy;
					dQ[1][5] = (-d_Qin_old[yp * 6 + 5] + 4 * d_Qin_old[ym * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idy;
				}

				if (zm == -1) {
					dQ[2][0] = 0;
					dQ[2][1] = 0;
					dQ[2][2] = 0;
					dQ[2][3] = 0;
					dQ[2][4] = 0;
					dQ[2][5] = 0;
				}
				else if (zp == -1) {
					dQ[2][0] = (d_Qin_old[zm * 6 + 0] - Qin[0]) * d_idz;
					dQ[2][1] = (d_Qin_old[zm * 6 + 1] - Qin[1]) * d_idz;
					dQ[2][2] = (d_Qin_old[zm * 6 + 2] - Qin[2]) * d_idz;
					dQ[2][3] = (d_Qin_old[zm * 6 + 3] - Qin[3]) * d_idz;
					dQ[2][4] = (d_Qin_old[zm * 6 + 4] - Qin[4]) * d_idz;
					dQ[2][5] = (d_Qin_old[zm * 6 + 5] - Qin[5]) * d_idz;
				}
				else {
					dQ[2][0] = (-d_Qin_old[zp * 6 + 0] + 4 * d_Qin_old[zm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idz;
					dQ[2][1] = (-d_Qin_old[zp * 6 + 1] + 4 * d_Qin_old[zm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idz;
					dQ[2][2] = (-d_Qin_old[zp * 6 + 2] + 4 * d_Qin_old[zm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idz;
					dQ[2][3] = (-d_Qin_old[zp * 6 + 3] + 4 * d_Qin_old[zm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idz;
					dQ[2][4] = (-d_Qin_old[zp * 6 + 4] + 4 * d_Qin_old[zm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idz;
					dQ[2][5] = (-d_Qin_old[zp * 6 + 5] + 4 * d_Qin_old[zm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idz;
				}

				Qelas[0] = dQ[0][0] * fabs(loc_nu[0]) + dQ[1][0] * fabs(loc_nu[1]) + dQ[2][0] * fabs(loc_nu[2]);
				Qelas[1] = dQ[0][1] * fabs(loc_nu[0]) + dQ[1][1] * fabs(loc_nu[1]) + dQ[2][1] * fabs(loc_nu[2]);
				Qelas[2] = dQ[0][2] * fabs(loc_nu[0]) + dQ[1][2] * fabs(loc_nu[1]) + dQ[2][2] * fabs(loc_nu[2]);
				Qelas[3] = dQ[0][3] * fabs(loc_nu[0]) + dQ[1][3] * fabs(loc_nu[1]) + dQ[2][3] * fabs(loc_nu[2]);
				Qelas[4] = dQ[0][4] * fabs(loc_nu[0]) + dQ[1][4] * fabs(loc_nu[1]) + dQ[2][4] * fabs(loc_nu[2]);
				Qelas[5] = dQ[0][5] * fabs(loc_nu[0]) + dQ[1][5] * fabs(loc_nu[1]) + dQ[2][5] * fabs(loc_nu[2]);
			}

			if (L2 != 0 || L3 != 0 || L4 != 0) {
				for (int j = 0; j < 3; j++) {
					if (loc_nu[j] < 0.) {
						for (int n = 0; n < 6; n++) {
							dQ[j][n] = -dQ[j][n];
						}
					}
				}
			}

			if(L2 != 0){
				const double temp0 = dQ[0][0] + dQ[1][1] + dQ[2][2];
				const double temp1 = dQ[0][1] + dQ[1][3] + dQ[2][4];
				const double temp2 = dQ[0][2] + dQ[1][4] + dQ[2][5];
				const double trace = (loc_nu[0] * temp0 + loc_nu[1] * temp1 + loc_nu[2] * temp2) * devThird;
				Qelas2[0] = loc_nu[0] * temp0 - trace;
				Qelas2[3] = loc_nu[1] * temp1 - trace;
				Qelas2[5] = loc_nu[2] * temp2 - trace;
				Qelas2[1] = 0.5 * (loc_nu[0] * temp1 + loc_nu[1] * temp0);
				Qelas2[2] = 0.5 * (loc_nu[0] * temp2 + loc_nu[2] * temp0);
				Qelas2[4] = 0.5 * (loc_nu[2] * temp1 + loc_nu[1] * temp2);
			}

			if (UseChiral) {
				Qch[0] = loc_nu[2] * Qin[1] - loc_nu[1] * Qin[2];
				Qch[3] = loc_nu[0] * Qin[4] - loc_nu[2] * Qin[1];
				Qch[5] = loc_nu[1] * Qin[2] - loc_nu[0] * Qin[4];
				Qch[1] = 0.5 * (loc_nu[2] * Qin[3] - loc_nu[1] * Qin[4] + loc_nu[0] * Qin[2] - loc_nu[2] * Qin[0]);
				Qch[2] = 0.5 * (loc_nu[2] * Qin[4] - loc_nu[1] * Qin[5] + loc_nu[1] * Qin[0] - loc_nu[0] * Qin[1]);
				Qch[4] = 0.5 * (loc_nu[0] * Qin[5] - loc_nu[2] * Qin[2] + loc_nu[1] * Qin[1] - loc_nu[0] * Qin[3]);
			}

		}

		if(degenerate == 1 || (signal == 12 || signal == 13)){
			relax_degen(Qin, loc_nu, Qdiff, S);
			for (int n = 0; n < 6; n++) {
				d_Qout[nv_indx * 6 + n] = Qin[n] + dt * (L1 * Qelas[n] + L2 * Qelas2[n] + (UseChiral ? 2. * qch * Qch[n] : 0.) - 2. * Wstr * Qdiff[n]);
			}
		}
		else if(degenerate == 2){
			relax_conic(Qin, loc_nu, Qdiff, S, tiltAngle);
			for (int n = 0; n < 6; n++) {
				d_Qout[nv_indx * 6 + n] = Qin[n] + dt * (L1 * Qelas[n] + L2 * Qelas2[n] + (UseChiral ? 2. * qch * Qch[n] : 0.) - 2. * Wstr * Qdiff[n]);
			}
		}
		else if(degenerate == 0 && infinite == 0){
			for (int n = 0; n < 6; n++) {
				d_Qout[nv_indx * 6 + n] = Qin[n] + dt * (L1 * Qelas[n] + L2 * Qelas2[n] + (UseChiral ? 2. * qch * Qch[n] : 0.) - Wstr * (Qin[n] - d_Qo[indx * 6 + n]));
			}
		}
		
	}
}

template __global__ void relax_surf<false>(const double* d_Qin_old, double* d_Qout, signed int* d_neighbor, unsigned int* d_Nvector_index, unsigned char* d_Nvector_signal, double* d_Qo,
	int chiral, double qch, double L1, double L2, double L3, double L4, double tiltAngle, unsigned int surf, int degenerate, int infinite, double Wstr, double Wp, double* d_nu, double d_idx,
	double d_idy, double d_idz, double dt, double S);
template __global__ void relax_surf<true>(const double* d_Qin_old, double* d_Qout, signed int* d_neighbor, unsigned int* d_Nvector_index, unsigned char* d_Nvector_signal, double* d_Qo,
	int chiral, double qch, double L1, double L2, double L3, double L4, double tiltAngle, unsigned int surf, int degenerate, int infinite, double Wstr, double Wp, double* d_nu, double d_idx,
	double d_idy, double d_idz, double dt, double S);

__device__ void relax_degen(double Qin[6], double loc_nu[3], double Qdiff[6], double S){
	double Qtemp[3][3] = { {0.} };
	double ptemp[3][3] = { {0.} };
	double proj[3][3] = { {0.} };
	double nuQnu = 0.;

	Qtemp[0][0] = Qin[0] + devThird * S;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + devThird * S;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + devThird * S;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) ptemp[i][j] = 1. - loc_nu[i] * loc_nu[j];
			else ptemp[i][j] = -loc_nu[i] * loc_nu[j];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int l = 0; l < 3; l++) {
				for (int m = 0; m < 3; m++) {
					proj[i][j] += ptemp[i][l] * Qtemp[l][m] * ptemp[m][j];
				}
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			nuQnu += loc_nu[i] * Qtemp[i][j] * loc_nu[j];
		}
	}
	nuQnu *= devThird;

	Qdiff[0] = Qtemp[0][0] - proj[0][0] - nuQnu;
	Qdiff[1] = Qtemp[0][1] - proj[0][1];
	Qdiff[2] = Qtemp[0][2] - proj[0][2];
	Qdiff[3] = Qtemp[1][1] - proj[1][1] - nuQnu;
	Qdiff[4] = Qtemp[1][2] - proj[1][2];
	Qdiff[5] = Qtemp[2][2] - proj[2][2] - nuQnu;
}

__device__ void relax_conic(double Qin[6], double loc_nu[3], double Qdiff[6], double S, double tiltAngle){
	double Qtemp[3][3] = { {0.} };
	double ptemp[3][3] = { {0.} };
	double proj[3][3] = { {0.} };

	Qtemp[0][0] = Qin[0] + devThird * S;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + devThird * S;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + devThird * S;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			ptemp[i][j] = loc_nu[i] * loc_nu[j];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int l = 0; l < 3; l++) {
				for (int m = 0; m < 3; m++) {
					proj[i][j] += ptemp[i][l] * Qtemp[l][m] * ptemp[m][j];
				}
			}
		}
	}

	const double cos_tilt = cos(tiltAngle / 180.0 * devPI);
	const double cosTiltAngleSq = cos_tilt * cos_tilt;

	Qdiff[0] = proj[0][0] - cosTiltAngleSq * S * ptemp[0][0];
	Qdiff[1] = proj[0][1] - cosTiltAngleSq * S * ptemp[0][1];
	Qdiff[2] = proj[0][2] - cosTiltAngleSq * S * ptemp[0][2];
	Qdiff[3] = proj[1][1] - cosTiltAngleSq * S * ptemp[1][1];
	Qdiff[4] = proj[1][2] - cosTiltAngleSq * S * ptemp[1][2];
	Qdiff[5] = proj[2][2] - cosTiltAngleSq * S * ptemp[2][2];

	const double trace = trace_f(Qdiff);
	Qdiff[0] -= trace;
	Qdiff[3] -= trace;
	Qdiff[5] -= trace;

}
