#include "relaxations.cuh"

//cuda variables. 
__device__ double fabs(double x);
__device__ double cos(double x);
__device__ double pow(double x, double y);
__device__ __constant__ double devThird = 1. / 3.;
__device__ __constant__ double devPI = 3.14159265358979323846;
__device__ double delta[6] = { 1., 0., 0., 1., 0., 1. };

//__device__ __constant__ double d_idx, d_idy, d_idz, d_iddx, d_iddy, d_iddz;

__device__ double trQQ_f(double Q[6]){

	return Q[0] * Q[0] + Q[3] * Q[3] + Q[5] * Q[5]\
			+ 2. * (Q[1] * Q[1] + Q[2] * Q[2] + Q[4] * Q[4]);

}

__device__ double trace_f(double Q[6]){

	return devThird * (Q[0] + Q[3] + Q[5]);

}

__global__ void relax_bulk(double* d_Qold, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, unsigned char* d_Qtensor_signal,
	double U, double U2, int chiral, double qch, double L1, double L2, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt) 
	{
	
	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;


	if (indx < bulk) {
		
		double Qin[6];
		double QQ[6];
		double Qldg[6];
		double trQQ;
		int xm, xp, ym, yp, zm, zp;
		double dQ[3][6] = { {0} };
		double ddQ[6][6] = { {0} };
		double Qelas[6] = { 0 };
		double Qelas2[6] = { 0 };
		double Qch[6] = { 0 };

		unsigned int Q_indx = d_Qtensor_index[indx];
		unsigned int Q_signal = d_Qtensor_signal[indx];

		if(d_bulktype[Q_indx] == 3){
			return;
		}

		// if(indx < 10){
		// 	printf("idx %lf idy %lf idz %lf\n", idx, idy, idz);
		// }

		// if(indx < 10){
		// 	printf("iddx %lf iddy %lf iddz %lf\n", iddx, iddy, iddz);
		// }

		//This could be converted to local memory
		Qin[0] = d_Qold[Q_indx * 6 + 0];
		Qin[1] = d_Qold[Q_indx * 6 + 1];
		Qin[2] = d_Qold[Q_indx * 6 + 2];
		Qin[3] = d_Qold[Q_indx * 6 + 3];
		Qin[4] = d_Qold[Q_indx * 6 + 4];
		Qin[5] = d_Qold[Q_indx * 6 + 5];

		//This could be registers.
		QQ[0] = Qin[0] * Qin[0] + Qin[1] * Qin[1] + Qin[2] * Qin[2];
		QQ[1] = Qin[0] * Qin[1] + Qin[1] * Qin[3] + Qin[2] * Qin[4];
		QQ[2] = Qin[0] * Qin[2] + Qin[1] * Qin[4] + Qin[2] * Qin[5];
		QQ[3] = Qin[1] * Qin[1] + Qin[3] * Qin[3] + Qin[4] * Qin[4];
		QQ[4] = Qin[1] * Qin[2] + Qin[3] * Qin[4] + Qin[4] * Qin[5];
		QQ[5] = Qin[2] * Qin[2] + Qin[4] * Qin[4] + Qin[5] * Qin[5];

		//Also a register.
		trQQ = trQQ_f(Qin);

		//We also need to change the thread access global memory and bring it to a lower level.
		if (d_bulktype[Q_indx] == 1) {
		
			Qldg[0] = (1. - U * devThird) * Qin[0] - U * (QQ[0] - trQQ * (Qin[0] + delta[0] * devThird));
			Qldg[1] = (1. - U * devThird) * Qin[1] - U * (QQ[1] - trQQ * (Qin[1] + delta[1] * devThird));
			Qldg[2] = (1. - U * devThird) * Qin[2] - U * (QQ[2] - trQQ * (Qin[2] + delta[2] * devThird));
			Qldg[3] = (1. - U * devThird) * Qin[3] - U * (QQ[3] - trQQ * (Qin[3] + delta[3] * devThird));
			Qldg[4] = (1. - U * devThird) * Qin[4] - U * (QQ[4] - trQQ * (Qin[4] + delta[4] * devThird));
			Qldg[5] = (1. - U * devThird) * Qin[5] - U * (QQ[5] - trQQ * (Qin[5] + delta[5] * devThird));

		}

		else if (d_bulktype[Q_indx] == 2) {

			Qldg[0] = (1. - U2 * devThird) * Qin[0] - U2 * (QQ[0] - trQQ * (Qin[0] + delta[0] * devThird));
			Qldg[1] = (1. - U2 * devThird) * Qin[1] - U2 * (QQ[1] - trQQ * (Qin[1] + delta[1] * devThird));
			Qldg[2] = (1. - U2 * devThird) * Qin[2] - U2 * (QQ[2] - trQQ * (Qin[2] + delta[2] * devThird));
			Qldg[3] = (1. - U2 * devThird) * Qin[3] - U2 * (QQ[3] - trQQ * (Qin[3] + delta[3] * devThird));
			Qldg[4] = (1. - U2 * devThird) * Qin[4] - U2 * (QQ[4] - trQQ * (Qin[4] + delta[4] * devThird));
			Qldg[5] = (1. - U2 * devThird) * Qin[5] - U2 * (QQ[5] - trQQ * (Qin[5] + delta[5] * devThird));

		}

		xm = d_neighbor[Q_indx * 6 + 0];
		xp = d_neighbor[Q_indx * 6 + 1];
		ym = d_neighbor[Q_indx * 6 + 2];
		yp = d_neighbor[Q_indx * 6 + 3];
		zm = d_neighbor[Q_indx * 6 + 4];
		zp = d_neighbor[Q_indx * 6 + 5];

		//if(indx < 20){
		//	printf("I'm thread %d and my neighbors are %d %d %d %d %d %d\n", indx, xm, xp, ym, yp, zm, zp);
		//}

		ddQ[0][0] = (d_Qold[xp * 6 + 0] + d_Qold[xm * 6 + 0] - 2. * Qin[0]) * iddx;
		ddQ[3][0] = (d_Qold[yp * 6 + 0] + d_Qold[ym * 6 + 0] - 2. * Qin[0]) * iddy;
		ddQ[5][0] = (d_Qold[zp * 6 + 0] + d_Qold[zm * 6 + 0] - 2. * Qin[0]) * iddz;
		Qelas[0] = ddQ[0][0] + ddQ[3][0] + ddQ[5][0];

		ddQ[0][1] = (d_Qold[xp * 6 + 1] + d_Qold[xm * 6 + 1] - 2. * Qin[1]) * iddx;
		ddQ[3][1] = (d_Qold[yp * 6 + 1] + d_Qold[ym * 6 + 1] - 2. * Qin[1]) * iddy;
		ddQ[5][1] = (d_Qold[zp * 6 + 1] + d_Qold[zm * 6 + 1] - 2. * Qin[1]) * iddz;
		Qelas[1] = ddQ[0][1] + ddQ[3][1] + ddQ[5][1];

		ddQ[0][2] = (d_Qold[xp * 6 + 2] + d_Qold[xm * 6 + 2] - 2. * Qin[2]) * iddx;
		ddQ[3][2] = (d_Qold[yp * 6 + 2] + d_Qold[ym * 6 + 2] - 2. * Qin[2]) * iddy;
		ddQ[5][2] = (d_Qold[zp * 6 + 2] + d_Qold[zm * 6 + 2] - 2. * Qin[2]) * iddz;
		Qelas[2] = ddQ[0][2] + ddQ[3][2] + ddQ[5][2];

		ddQ[0][3] = (d_Qold[xp * 6 + 3] + d_Qold[xm * 6 + 3] - 2. * Qin[3]) * iddx;
		ddQ[3][3] = (d_Qold[yp * 6 + 3] + d_Qold[ym * 6 + 3] - 2. * Qin[3]) * iddy;
		ddQ[5][3] = (d_Qold[zp * 6 + 3] + d_Qold[zm * 6 + 3] - 2. * Qin[3]) * iddz;
		Qelas[3] = ddQ[0][3] + ddQ[3][3] + ddQ[5][3];

		ddQ[0][4] = (d_Qold[xp * 6 + 4] + d_Qold[xm * 6 + 4] - 2. * Qin[4]) * iddx;
		ddQ[3][4] = (d_Qold[yp * 6 + 4] + d_Qold[ym * 6 + 4] - 2. * Qin[4]) * iddy;
		ddQ[5][4] = (d_Qold[zp * 6 + 4] + d_Qold[zm * 6 + 4] - 2. * Qin[4]) * iddz;
		Qelas[4] = ddQ[0][4] + ddQ[3][4] + ddQ[5][4];

		ddQ[0][5] = (d_Qold[xp * 6 + 5] + d_Qold[xm * 6 + 5] - 2. * Qin[5]) * iddx;
		ddQ[3][5] = (d_Qold[yp * 6 + 5] + d_Qold[ym * 6 + 5] - 2. * Qin[5]) * iddy;
		ddQ[5][5] = (d_Qold[zp * 6 + 5] + d_Qold[zm * 6 + 5] - 2. * Qin[5]) * iddz;
		Qelas[5] = ddQ[0][5] + ddQ[3][5] + ddQ[5][5];

		//printf("L1:%lf L2:%lf L3:%lf L4:%lf", L1_dev, L2_dev, L3_dev, L4_dev);

		if (chiral == 1) {

			dQ[0][0] = (d_Qold[xp * 6 + 0] - d_Qold[xm * 6 + 0]) * 0.5 * idx;
			dQ[1][0] = (d_Qold[yp * 6 + 0] - d_Qold[ym * 6 + 0]) * 0.5 * idy;
			dQ[2][0] = (d_Qold[zp * 6 + 0] - d_Qold[zm * 6 + 0]) * 0.5 * idz;

			dQ[0][1] = (d_Qold[xp * 6 + 1] - d_Qold[xm * 6 + 1]) * 0.5 * idx;
			dQ[1][1] = (d_Qold[yp * 6 + 1] - d_Qold[ym * 6 + 1]) * 0.5 * idy;
			dQ[2][1] = (d_Qold[zp * 6 + 1] - d_Qold[zm * 6 + 1]) * 0.5 * idz;

			dQ[0][2] = (d_Qold[xp * 6 + 2] - d_Qold[xm * 6 + 2]) * 0.5 * idx;
			dQ[1][2] = (d_Qold[yp * 6 + 2] - d_Qold[ym * 6 + 2]) * 0.5 * idy;
			dQ[2][2] = (d_Qold[zp * 6 + 2] - d_Qold[zm * 6 + 2]) * 0.5 * idz;

			dQ[0][3] = (d_Qold[xp * 6 + 3] - d_Qold[xm * 6 + 3]) * 0.5 * idx;
			dQ[1][3] = (d_Qold[yp * 6 + 3] - d_Qold[ym * 6 + 3]) * 0.5 * idy;
			dQ[2][3] = (d_Qold[zp * 6 + 3] - d_Qold[zm * 6 + 3]) * 0.5 * idz;

			dQ[0][4] = (d_Qold[xp * 6 + 4] - d_Qold[xm * 6 + 4]) * 0.5 * idx;
			dQ[1][4] = (d_Qold[yp * 6 + 4] - d_Qold[ym * 6 + 4]) * 0.5 * idy;
			dQ[2][4] = (d_Qold[zp * 6 + 4] - d_Qold[zm * 6 + 4]) * 0.5 * idz;

			dQ[0][5] = (d_Qold[xp * 6 + 5] - d_Qold[xm * 6 + 5]) * 0.5 * idx;
			dQ[1][5] = (d_Qold[yp * 6 + 5] - d_Qold[ym * 6 + 5]) * 0.5 * idy;
			dQ[2][5] = (d_Qold[zp * 6 + 5] - d_Qold[zm * 6 + 5]) * 0.5 * idz;
			
		}

		//if((L2 + L4) != 0 || L3 != 0){
		if(L2 != 0){
			if(Q_signal == 0 && d_bulktype[Q_indx] != 13){
				
				//neighbor xpyp is the yp neighbor of xp: neighbor[xp * 6 + 3]; same definition for other points
				ddQ[1][0] = (d_Qold[d_neighbor[xp * 6 + 3] * 6 + 0] + d_Qold[d_neighbor[xm * 6 + 2] * 6 + 0] - d_Qold[d_neighbor[xm * 6 + 3] * 6 + 0] - d_Qold[d_neighbor[xp * 6 + 2] * 6 + 0]) * idx * idy * 0.25;
				ddQ[2][0] = (d_Qold[d_neighbor[xp * 6 + 5] * 6 + 0] + d_Qold[d_neighbor[xm * 6 + 4] * 6 + 0] - d_Qold[d_neighbor[xm * 6 + 5] * 6 + 0] - d_Qold[d_neighbor[xp * 6 + 4] * 6 + 0]) * idx * idz * 0.25;
				ddQ[4][0] = (d_Qold[d_neighbor[yp * 6 + 5] * 6 + 0] + d_Qold[d_neighbor[ym * 6 + 4] * 6 + 0] - d_Qold[d_neighbor[ym * 6 + 5] * 6 + 0] - d_Qold[d_neighbor[yp * 6 + 4] * 6 + 0]) * idy * idz * 0.25;

				ddQ[1][1] = (d_Qold[d_neighbor[xp * 6 + 3] * 6 + 1] + d_Qold[d_neighbor[xm * 6 + 2] * 6 + 1] - d_Qold[d_neighbor[xm * 6 + 3] * 6 + 1] - d_Qold[d_neighbor[xp * 6 + 2] * 6 + 1]) * idx * idy * 0.25;
				ddQ[2][1] = (d_Qold[d_neighbor[xp * 6 + 5] * 6 + 1] + d_Qold[d_neighbor[xm * 6 + 4] * 6 + 1] - d_Qold[d_neighbor[xm * 6 + 5] * 6 + 1] - d_Qold[d_neighbor[xp * 6 + 4] * 6 + 1]) * idx * idz * 0.25;
				ddQ[4][1] = (d_Qold[d_neighbor[yp * 6 + 5] * 6 + 1] + d_Qold[d_neighbor[ym * 6 + 4] * 6 + 1] - d_Qold[d_neighbor[ym * 6 + 5] * 6 + 1] - d_Qold[d_neighbor[yp * 6 + 4] * 6 + 1]) * idy * idz * 0.25;
				
				ddQ[1][2] = (d_Qold[d_neighbor[xp * 6 + 3] * 6 + 2] + d_Qold[d_neighbor[xm * 6 + 2] * 6 + 2] - d_Qold[d_neighbor[xm * 6 + 3] * 6 + 2] - d_Qold[d_neighbor[xp * 6 + 2] * 6 + 2]) * idx * idy * 0.25;
				ddQ[2][2] = (d_Qold[d_neighbor[xp * 6 + 5] * 6 + 2] + d_Qold[d_neighbor[xm * 6 + 4] * 6 + 2] - d_Qold[d_neighbor[xm * 6 + 5] * 6 + 2] - d_Qold[d_neighbor[xp * 6 + 4] * 6 + 2]) * idx * idz * 0.25;
				ddQ[4][2] = (d_Qold[d_neighbor[yp * 6 + 5] * 6 + 2] + d_Qold[d_neighbor[ym * 6 + 4] * 6 + 2] - d_Qold[d_neighbor[ym * 6 + 5] * 6 + 2] - d_Qold[d_neighbor[yp * 6 + 4] * 6 + 2]) * idy * idz * 0.25;

				ddQ[1][3] = (d_Qold[d_neighbor[xp * 6 + 3] * 6 + 3] + d_Qold[d_neighbor[xm * 6 + 2] * 6 + 3] - d_Qold[d_neighbor[xm * 6 + 3] * 6 + 3] - d_Qold[d_neighbor[xp * 6 + 2] * 6 + 3]) * idx * idy * 0.25;
				ddQ[2][3] = (d_Qold[d_neighbor[xp * 6 + 5] * 6 + 3] + d_Qold[d_neighbor[xm * 6 + 4] * 6 + 3] - d_Qold[d_neighbor[xm * 6 + 5] * 6 + 3] - d_Qold[d_neighbor[xp * 6 + 4] * 6 + 3]) * idx * idz * 0.25;
				ddQ[4][3] = (d_Qold[d_neighbor[yp * 6 + 5] * 6 + 3] + d_Qold[d_neighbor[ym * 6 + 4] * 6 + 3] - d_Qold[d_neighbor[ym * 6 + 5] * 6 + 3] - d_Qold[d_neighbor[yp * 6 + 4] * 6 + 3]) * idy * idz * 0.25;

				ddQ[1][4] = (d_Qold[d_neighbor[xp * 6 + 3] * 6 + 4] + d_Qold[d_neighbor[xm * 6 + 2] * 6 + 4] - d_Qold[d_neighbor[xm * 6 + 3] * 6 + 4] - d_Qold[d_neighbor[xp * 6 + 2] * 6 + 4]) * idx * idy * 0.25;
				ddQ[2][4] = (d_Qold[d_neighbor[xp * 6 + 5] * 6 + 4] + d_Qold[d_neighbor[xm * 6 + 4] * 6 + 4] - d_Qold[d_neighbor[xm * 6 + 5] * 6 + 4] - d_Qold[d_neighbor[xp * 6 + 4] * 6 + 4]) * idx * idz * 0.25;
				ddQ[4][4] = (d_Qold[d_neighbor[yp * 6 + 5] * 6 + 4] + d_Qold[d_neighbor[ym * 6 + 4] * 6 + 4] - d_Qold[d_neighbor[ym * 6 + 5] * 6 + 4] - d_Qold[d_neighbor[yp * 6 + 4] * 6 + 4]) * idy * idz * 0.25;

				ddQ[1][5] = (d_Qold[d_neighbor[xp * 6 + 3] * 6 + 5] + d_Qold[d_neighbor[xm * 6 + 2] * 6 + 5] - d_Qold[d_neighbor[xm * 6 + 3] * 6 + 5] - d_Qold[d_neighbor[xp * 6 + 2] * 6 + 5]) * idx * idy * 0.25;
				ddQ[2][5] = (d_Qold[d_neighbor[xp * 6 + 5] * 6 + 5] + d_Qold[d_neighbor[xm * 6 + 4] * 6 + 5] - d_Qold[d_neighbor[xm * 6 + 5] * 6 + 5] - d_Qold[d_neighbor[xp * 6 + 4] * 6 + 5]) * idx * idz * 0.25;
				ddQ[4][5] = (d_Qold[d_neighbor[yp * 6 + 5] * 6 + 5] + d_Qold[d_neighbor[ym * 6 + 4] * 6 + 5] - d_Qold[d_neighbor[ym * 6 + 5] * 6 + 5] - d_Qold[d_neighbor[yp * 6 + 4] * 6 + 5]) * idy * idz * 0.25;

			}
			else{
				
				ddQ[1][0] = 0;
				ddQ[2][0] = 0;
				ddQ[4][0] = 0;

				ddQ[1][1] = 0;
				ddQ[2][1] = 0;
				ddQ[4][1] = 0;

				ddQ[1][2] = 0;
				ddQ[2][2] = 0;
				ddQ[4][2] = 0;

				ddQ[1][3] = 0;
				ddQ[2][3] = 0;
				ddQ[4][3] = 0;

				ddQ[1][4] = 0;
				ddQ[2][4] = 0;
				ddQ[4][4] = 0;

				ddQ[1][5] = 0;
				ddQ[2][5] = 0;
				ddQ[4][5] = 0;
				
			}

			Qelas2[0] = ddQ[0][0] + ddQ[1][1] + ddQ[2][2];
			Qelas2[1] = 0.5 * (ddQ[1][0] + ddQ[0][1] + ddQ[1][3] + ddQ[3][1] + ddQ[2][4] + ddQ[4][2]);
			Qelas2[2] = 0.5 * (ddQ[2][0] + ddQ[0][2] + ddQ[1][4] + ddQ[4][1] + ddQ[2][5] + ddQ[5][2]);
			Qelas2[3] = ddQ[1][1] + ddQ[3][3] + ddQ[4][4];
			Qelas2[4] = 0.5 * (ddQ[1][2] + ddQ[2][1] + ddQ[5][4] + ddQ[4][5] + ddQ[3][4] + ddQ[4][3]);
			Qelas2[5] = ddQ[2][2] + ddQ[4][4] + ddQ[5][5];
			Qelas2[0] -= trace_f(Qelas2);
			Qelas2[3] -= trace_f(Qelas2);
			Qelas2[5] -= trace_f(Qelas2);

		}
		/*if((L2 + L4) != 0){
			Qelas2[0] = ddQ[0][0] + ddQ[1][1] + ddQ[2][2];
			Qelas2[1] = 0.5 * (ddQ[1][0] + ddQ[0][1] + ddQ[1][3] + ddQ[3][1] + ddQ[2][4] + ddQ[4][2]);
			Qelas2[2] = 0.5 * (ddQ[2][0] + ddQ[0][2] + ddQ[1][4] + ddQ[4][1] + ddQ[2][5] + ddQ[5][2]);
			Qelas2[3] = ddQ[1][1] + ddQ[3][3] + ddQ[4][4];
			Qelas2[4] = 0.5 * (ddQ[1][2] + ddQ[2][1] + ddQ[5][4] + ddQ[4][5] + ddQ[3][4] + ddQ[4][3]);
			Qelas2[5] = ddQ[2][2] + ddQ[4][4] + ddQ[5][5];
			//trace = (Qelas2[0] + Qelas2[3] + Qelas2[5]) * third;
			//Qelas2[0] -= trace;
			//Qelas2[3] -= trace;
			//Qelas2[5] -= trace;

			Qelas2[0] -= trace_f(Qelas2);
			Qelas2[3] -= trace_f(Qelas2);
			Qelas2[5] -= trace_f(Qelas2);

		}
		/*if(L3 != 0){

			Qelas3[0] = - 0.5 * trqq(dQ[0]);
			Qelas3[1] = - 0.5 * (q_mult(dQ[0], dQ[1]));
			Qelas3[2] = - 0.5 * (q_mult(dQ[0], dQ[2]));
			Qelas3[3] = - 0.5 * trqq(dQ[1]);
			Qelas3[4] = - 0.5 * (q_mult(dQ[1], dQ[2]));
			Qelas3[5] = - 0.5 * trqq(dQ[2]);

			for (int n = 0; n < 6; n++){
				Qelas3[n] += Qin[0] * ddQ[0][n] + Qin[3] * ddQ[3][n] + Qin[5] * ddQ[5][n] + 2 * (Qin[1] * ddQ[1][n] + Qin[2] * ddQ[2][n] + Qin[4] * ddQ[4][n]);
				Qelas3[n] += dQ[0][n] * (dQ[0][0] + dQ[1][1] + dQ[2][2]) + dQ[1][n] * (dQ[0][1] + dQ[1][3] + dQ[2][4]) + dQ[2][n] * (dQ[0][2] + dQ[1][4] + dQ[2][5]);
			}		

			trace = (Qelas3[0] + Qelas3[3] + Qelas3[5]) * third;
			Qelas3[0] -= trace;
			Qelas3[3] -= trace;
			Qelas3[5] -= trace;
		} */

		if(chiral == 1) {
			Qch[0] = 2. * (dQ[1][2] - dQ[2][1]);
			Qch[3] = 2. * (dQ[2][1] - dQ[0][4]);
			Qch[5] = 2. * (dQ[0][4] - dQ[1][2]);
			Qch[1] = dQ[1][4] - dQ[2][3] + dQ[2][0] - dQ[0][2];
			Qch[2] = dQ[1][5] - dQ[2][4] + dQ[0][1] - dQ[1][0];
			Qch[4] = dQ[2][2] - dQ[0][5] + dQ[0][3] - dQ[1][1];
		} 

		//printf("Done %d %d", indx, Q_signal);
		__syncthreads();

		//d_Qold[Q_indx * 6 + 0] = Qin[0] + dt * (-Qldg[0] + L1 * Qelas[0] + (L2 + L4) * Qelas2[0] - 2. * (double)chiral * qch * L1 * Qch[0]);

		d_Qold[Q_indx * 6 + 0] = Qin[0] + dt * (-Qldg[0] + L1 * Qelas[0] + L2 * Qelas2[0] - 2. * (double)chiral * qch * L1 * Qch[0]);
		d_Qold[Q_indx * 6 + 1] = Qin[1] + dt * (-Qldg[1] + L1 * Qelas[1] + L2 * Qelas2[1] - 2. * (double)chiral * qch * L1 * Qch[1]);
		d_Qold[Q_indx * 6 + 2] = Qin[2] + dt * (-Qldg[2] + L1 * Qelas[2] + L2 * Qelas2[2] - 2. * (double)chiral * qch * L1 * Qch[2]);
		d_Qold[Q_indx * 6 + 3] = Qin[3] + dt * (-Qldg[3] + L1 * Qelas[3] + L2 * Qelas2[3] - 2. * (double)chiral * qch * L1 * Qch[3]);
		d_Qold[Q_indx * 6 + 4] = Qin[4] + dt * (-Qldg[4] + L1 * Qelas[4] + L2 * Qelas2[4] - 2. * (double)chiral * qch * L1 * Qch[4]);
		d_Qold[Q_indx * 6 + 5] = Qin[5] + dt * (-Qldg[5] + L1 * Qelas[5] + L2 * Qelas2[5] - 2. * (double)chiral * qch * L1 * Qch[5]);
		
	}
 

}

//variable degenerate it's only relevant for the system surface. For nanoparticle surface, we will use degen and signal types.
__global__ void relax_surf(double* d_Qold, signed int* d_neighbor, unsigned int* d_Nvector_index, unsigned char* d_Nvector_signal, double* d_Qo, 
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

		unsigned char degen = 0;
		double loc_nu[3] = { 0. };
		double Qin[6] = { 0. };
		double Qdiff[6] = { 0. };
		int xm, xp, ym, yp, zm, zp;

		double temp[3] = { 0. };
		double dQ[3][6] = { { 0. } };
		double Qelas[6] = { 0. };
		double Qelas2[6] = { 0. };
		double Qch[6] = { 0 };
		double trace = 0.;

		

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

			Qin[0] = d_Qold[nv_indx * 6 + 0];
			Qin[1] = d_Qold[nv_indx * 6 + 1];
			Qin[2] = d_Qold[nv_indx * 6 + 2];
			Qin[3] = d_Qold[nv_indx * 6 + 3];
			Qin[4] = d_Qold[nv_indx * 6 + 4];
			Qin[5] = d_Qold[nv_indx * 6 + 5];

			xm = d_neighbor[nv_indx * 6 + 0];
			xp = d_neighbor[nv_indx * 6 + 1];
			ym = d_neighbor[nv_indx * 6 + 2];
			yp = d_neighbor[nv_indx * 6 + 3];
			zm = d_neighbor[nv_indx * 6 + 4];
			zp = d_neighbor[nv_indx * 6 + 5];

			if ((signal % 2) == 0) {

				dQ[0][0] = (-d_Qold[xp * 6 + 0] + 4 * d_Qold[xm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idx;
				dQ[1][0] = (-d_Qold[yp * 6 + 0] + 4 * d_Qold[ym * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idy;
				dQ[2][0] = (-d_Qold[zp * 6 + 0] + 4 * d_Qold[zm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idz;
				Qelas[0] = dQ[0][0] * fabs(loc_nu[0]) + dQ[1][0] * fabs(loc_nu[1]) + dQ[2][0] * fabs(loc_nu[2]);

				dQ[0][1] = (-d_Qold[xp * 6 + 1] + 4 * d_Qold[xm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idx;
				dQ[1][1] = (-d_Qold[yp * 6 + 1] + 4 * d_Qold[ym * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idy;
				dQ[2][1] = (-d_Qold[zp * 6 + 1] + 4 * d_Qold[zm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idz;
				Qelas[1] = dQ[0][1] * fabs(loc_nu[0]) + dQ[1][1] * fabs(loc_nu[1]) + dQ[2][1] * fabs(loc_nu[2]);

				dQ[0][2] = (-d_Qold[xp * 6 + 2] + 4 * d_Qold[xm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idx;
				dQ[1][2] = (-d_Qold[yp * 6 + 2] + 4 * d_Qold[ym * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idy;
				dQ[2][2] = (-d_Qold[zp * 6 + 2] + 4 * d_Qold[zm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idz;
				Qelas[2] = dQ[0][2] * fabs(loc_nu[0]) + dQ[1][2] * fabs(loc_nu[1]) + dQ[2][2] * fabs(loc_nu[2]);

				dQ[0][3] = (-d_Qold[xp * 6 + 3] + 4 * d_Qold[xm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idx;
				dQ[1][3] = (-d_Qold[yp * 6 + 3] + 4 * d_Qold[ym * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idy;
				dQ[2][3] = (-d_Qold[zp * 6 + 3] + 4 * d_Qold[zm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idz;
				Qelas[3] = dQ[0][3] * fabs(loc_nu[0]) + dQ[1][3] * fabs(loc_nu[1]) + dQ[2][3] * fabs(loc_nu[2]);

				dQ[0][4] = (-d_Qold[xp * 6 + 4] + 4 * d_Qold[xm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idx;
				dQ[1][4] = (-d_Qold[yp * 6 + 4] + 4 * d_Qold[ym * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idy;
				dQ[2][4] = (-d_Qold[zp * 6 + 4] + 4 * d_Qold[zm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idz;
				Qelas[4] = dQ[0][4] * fabs(loc_nu[0]) + dQ[1][4] * fabs(loc_nu[1]) + dQ[2][4] * fabs(loc_nu[2]);

				dQ[0][5] = (-d_Qold[xp * 6 + 5] + 4 * d_Qold[xm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idx;
				dQ[1][5] = (-d_Qold[yp * 6 + 5] + 4 * d_Qold[ym * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idy;
				dQ[2][5] = (-d_Qold[zp * 6 + 5] + 4 * d_Qold[zm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idz;
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

					dQ[0][0] = (d_Qold[xm * 6 + 0] - Qin[0]) * d_idx;
					dQ[0][1] = (d_Qold[xm * 6 + 1] - Qin[1]) * d_idx;
					dQ[0][2] = (d_Qold[xm * 6 + 2] - Qin[2]) * d_idx;
					dQ[0][3] = (d_Qold[xm * 6 + 3] - Qin[3]) * d_idx;
					dQ[0][4] = (d_Qold[xm * 6 + 4] - Qin[4]) * d_idx;
					dQ[0][5] = (d_Qold[xm * 6 + 5] - Qin[5]) * d_idx;

				}
				else {

					dQ[0][0] = (-d_Qold[xp * 6 + 0] + 4 * d_Qold[xm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idx;
					dQ[0][1] = (-d_Qold[xp * 6 + 1] + 4 * d_Qold[xm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idx;
					dQ[0][2] = (-d_Qold[xp * 6 + 2] + 4 * d_Qold[xm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idx;
					dQ[0][3] = (-d_Qold[xp * 6 + 3] + 4 * d_Qold[xm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idx;
					dQ[0][4] = (-d_Qold[xp * 6 + 4] + 4 * d_Qold[xm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idx;
					dQ[0][5] = (-d_Qold[xp * 6 + 5] + 4 * d_Qold[xm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idx;

				}

				if (ym == -1) {

					dQ[1][0] = 0;
					dQ[1][1] = 0;
					dQ[1][2] = 0;
					dQ[1][3] = 0;
					dQ[1][4] = 0;
					dQ[1][5] = 0;

				}
				else if (yp  == -1) {

					dQ[1][0] = (d_Qold[ym * 6 + 0] - Qin[0]) * d_idy;
					dQ[1][1] = (d_Qold[ym * 6 + 1] - Qin[1]) * d_idy;
					dQ[1][2] = (d_Qold[ym * 6 + 2] - Qin[2]) * d_idy;
					dQ[1][3] = (d_Qold[ym * 6 + 3] - Qin[3]) * d_idy;
					dQ[1][4] = (d_Qold[ym * 6 + 4] - Qin[4]) * d_idy;
					dQ[1][5] = (d_Qold[ym * 6 + 5] - Qin[5]) * d_idy;

				}
				else {

					dQ[1][0] = (-d_Qold[yp * 6 + 0] + 4 * d_Qold[ym * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idy;
					dQ[1][1] = (-d_Qold[yp * 6 + 1] + 4 * d_Qold[ym * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idy;
					dQ[1][2] = (-d_Qold[yp * 6 + 2] + 4 * d_Qold[ym * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idy;
					dQ[1][3] = (-d_Qold[yp * 6 + 3] + 4 * d_Qold[ym * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idy;
					dQ[1][4] = (-d_Qold[yp * 6 + 4] + 4 * d_Qold[ym * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idy;
					dQ[1][5] = (-d_Qold[yp * 6 + 5] + 4 * d_Qold[ym * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idy;

				}

				if (zm == -1){
					
					dQ[2][0] = 0;
					dQ[2][1] = 0;
					dQ[2][2] = 0;
					dQ[2][3] = 0;
					dQ[2][4] = 0;
					dQ[2][5] = 0;

				}
				else if (zp == -1) {

					dQ[2][0] = (d_Qold[zm * 6 + 0] - Qin[0]) * d_idz;
					dQ[2][1] = (d_Qold[zm * 6 + 1] - Qin[1]) * d_idz;
					dQ[2][2] = (d_Qold[zm * 6 + 2] - Qin[2]) * d_idz;
					dQ[2][3] = (d_Qold[zm * 6 + 3] - Qin[3]) * d_idz;
					dQ[2][4] = (d_Qold[zm * 6 + 4] - Qin[4]) * d_idz;
					dQ[2][5] = (d_Qold[zm * 6 + 5] - Qin[5]) * d_idz;

				}
				else {

					dQ[2][0] = (-d_Qold[zp * 6 + 0] + 4 * d_Qold[zm * 6 + 0] - 3 * Qin[0]) * 0.5 * d_idz;
					dQ[2][1] = (-d_Qold[zp * 6 + 1] + 4 * d_Qold[zm * 6 + 1] - 3 * Qin[1]) * 0.5 * d_idz;
					dQ[2][2] = (-d_Qold[zp * 6 + 2] + 4 * d_Qold[zm * 6 + 2] - 3 * Qin[2]) * 0.5 * d_idz;
					dQ[2][3] = (-d_Qold[zp * 6 + 3] + 4 * d_Qold[zm * 6 + 3] - 3 * Qin[3]) * 0.5 * d_idz;
					dQ[2][4] = (-d_Qold[zp * 6 + 4] + 4 * d_Qold[zm * 6 + 4] - 3 * Qin[4]) * 0.5 * d_idz;
					dQ[2][5] = (-d_Qold[zp * 6 + 5] + 4 * d_Qold[zm * 6 + 5] - 3 * Qin[5]) * 0.5 * d_idz;

				}

				Qelas[0] = dQ[0][0] * fabs(loc_nu[0]) + dQ[1][0] * fabs(loc_nu[1]) + dQ[2][0] * fabs(loc_nu[2]);
				Qelas[1] = dQ[0][1] * fabs(loc_nu[0]) + dQ[1][1] * fabs(loc_nu[1]) + dQ[2][1] * fabs(loc_nu[2]);
				Qelas[2] = dQ[0][2] * fabs(loc_nu[0]) + dQ[1][2] * fabs(loc_nu[1]) + dQ[2][2] * fabs(loc_nu[2]);
				Qelas[3] = dQ[0][3] * fabs(loc_nu[0]) + dQ[1][3] * fabs(loc_nu[1]) + dQ[2][3] * fabs(loc_nu[2]);
				Qelas[4] = dQ[0][4] * fabs(loc_nu[0]) + dQ[1][4] * fabs(loc_nu[1]) + dQ[2][4] * fabs(loc_nu[2]);
				Qelas[5] = dQ[0][5] * fabs(loc_nu[0]) + dQ[1][5] * fabs(loc_nu[1]) + dQ[2][5] * fabs(loc_nu[2]);
				
			}

			//if(L2 != 0 || L3 != 0 || L4 != 0){
			if(L2 != 0){

				if(loc_nu[0] < 0.){
					dQ[0][0] = -dQ[0][0];
					dQ[0][1] = -dQ[0][1];
					dQ[0][2] = -dQ[0][2];
				}
				if(loc_nu[1] < 0.){
					dQ[1][0] = -dQ[1][0];
					dQ[1][1] = -dQ[1][1];
					dQ[1][2] = -dQ[1][2];
				}
				if(loc_nu[2] < 0.){
					dQ[2][0] = -dQ[2][0];
					dQ[2][1] = -dQ[2][1];
					dQ[2][2] = -dQ[2][2];
				}

			}

			if(L2 != 0){
				
				temp[0] = dQ[0][0] + dQ[1][1] + dQ[2][2];
				temp[1] = dQ[0][1] + dQ[1][3] + dQ[2][4];
				temp[2] = dQ[0][2] + dQ[1][4] + dQ[2][5];
				trace = (loc_nu[0] * temp[0] + loc_nu[1] * temp[1] + loc_nu[2] * temp[2]) * devThird;
				Qelas2[0] = loc_nu[0] * temp[0] - trace;
				Qelas2[3] = loc_nu[1] * temp[1] - trace;
				Qelas2[5] = loc_nu[2] * temp[2] - trace;
				Qelas2[1] = 0.5 * (loc_nu[0] * temp[1] + loc_nu[1] * temp[0]);
				Qelas2[2] = 0.5 * (loc_nu[0] * temp[2] + loc_nu[2] * temp[0]);
				Qelas2[4] = 0.5 * (loc_nu[2] * temp[1] + loc_nu[1] * temp[2]);
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
		
		if(degenerate == 1 || (signal == 12 || signal == 13)){

			relax_degen(Qin, loc_nu, Qdiff, S);

			d_Qold[d_Nvector_index[indx] * 6 + 0] = Qin[0] + dt * (L1 * Qelas[0] + L2 * Qelas2[0] + (double)chiral * 2. * qch * Qch[0] - 2. * Wstr * Qdiff[0]);
			d_Qold[d_Nvector_index[indx] * 6 + 1] = Qin[1] + dt * (L1 * Qelas[1] + L2 * Qelas2[1] + (double)chiral * 2. * qch * Qch[1] - 2. * Wstr * Qdiff[1]);
			d_Qold[d_Nvector_index[indx] * 6 + 2] = Qin[2] + dt * (L1 * Qelas[2] + L2 * Qelas2[2] + (double)chiral * 2. * qch * Qch[2] - 2. * Wstr * Qdiff[2]);
			d_Qold[d_Nvector_index[indx] * 6 + 3] = Qin[3] + dt * (L1 * Qelas[3] + L2 * Qelas2[3] + (double)chiral * 2. * qch * Qch[3] - 2. * Wstr * Qdiff[3]);
			d_Qold[d_Nvector_index[indx] * 6 + 4] = Qin[4] + dt * (L1 * Qelas[4] + L2 * Qelas2[4] + (double)chiral * 2. * qch * Qch[4] - 2. * Wstr * Qdiff[4]);
			d_Qold[d_Nvector_index[indx] * 6 + 5] = Qin[5] + dt * (L1 * Qelas[5] + L2 * Qelas2[5] + (double)chiral * 2. * qch * Qch[5] - 2. * Wstr * Qdiff[5]);
			
		}

		else if(degenerate == 2){

			relax_conic(Qin, loc_nu, Qdiff, S, tiltAngle);

			d_Qold[d_Nvector_index[indx] * 6 + 0] = Qin[0] + dt * (L1 * Qelas[0] + L2 * Qelas2[0] + (double)chiral * 2. * qch * Qch[0] - 2. * Wstr * Qdiff[0]);
			d_Qold[d_Nvector_index[indx] * 6 + 1] = Qin[1] + dt * (L1 * Qelas[1] + L2 * Qelas2[1] + (double)chiral * 2. * qch * Qch[1] - 2. * Wstr * Qdiff[1]);
			d_Qold[d_Nvector_index[indx] * 6 + 2] = Qin[2] + dt * (L1 * Qelas[2] + L2 * Qelas2[2] + (double)chiral * 2. * qch * Qch[2] - 2. * Wstr * Qdiff[2]);
			d_Qold[d_Nvector_index[indx] * 6 + 3] = Qin[3] + dt * (L1 * Qelas[3] + L2 * Qelas2[3] + (double)chiral * 2. * qch * Qch[3] - 2. * Wstr * Qdiff[3]);
			d_Qold[d_Nvector_index[indx] * 6 + 4] = Qin[4] + dt * (L1 * Qelas[4] + L2 * Qelas2[4] + (double)chiral * 2. * qch * Qch[4] - 2. * Wstr * Qdiff[4]);
			d_Qold[d_Nvector_index[indx] * 6 + 5] = Qin[5] + dt * (L1 * Qelas[5] + L2 * Qelas2[5] + (double)chiral * 2. * qch * Qch[5] - 2. * Wstr * Qdiff[5]);

		}

		else if(degenerate == 0 && infinite == 0){

			d_Qold[d_Nvector_index[indx] * 6 + 0] = Qin[0] + dt * (L1 * Qelas[0] + L2 * Qelas2[0] + chiral * 2 * qch * Qch[0] - Wstr * (Qin[0] - d_Qo[indx * 6 + 0]));
			d_Qold[d_Nvector_index[indx] * 6 + 1] = Qin[1] + dt * (L1 * Qelas[1] + L2 * Qelas2[1] + chiral * 2 * qch * Qch[1] - Wstr * (Qin[1] - d_Qo[indx * 6 + 1]));
			d_Qold[d_Nvector_index[indx] * 6 + 2] = Qin[2] + dt * (L1 * Qelas[2] + L2 * Qelas2[2] + chiral * 2 * qch * Qch[2] - Wstr * (Qin[2] - d_Qo[indx * 6 + 2]));
			d_Qold[d_Nvector_index[indx] * 6 + 3] = Qin[3] + dt * (L1 * Qelas[3] + L2 * Qelas2[3] + chiral * 2 * qch * Qch[3] - Wstr * (Qin[3] - d_Qo[indx * 6 + 3]));
			d_Qold[d_Nvector_index[indx] * 6 + 4] = Qin[4] + dt * (L1 * Qelas[4] + L2 * Qelas2[4] + chiral * 2 * qch * Qch[4] - Wstr * (Qin[4] - d_Qo[indx * 6 + 4]));
			d_Qold[d_Nvector_index[indx] * 6 + 5] = Qin[5] + dt * (L1 * Qelas[5] + L2 * Qelas2[5] + chiral * 2 * qch * Qch[5] - Wstr * (Qin[5] - d_Qo[indx * 6 + 5]));

		}
		
	}
}

__device__ void relax_degen(double Qin[6], double loc_nu[3], double Qdiff[6], double S){
	double Qtemp[3][3] = {{ 0. }};
	double ptemp[3][3] = {{ 0. }};
	double Qp[3][3] = {{ 0. }};
	double nuQnu = 0;
	Qtemp[0][0] = Qin[0] + devThird * S;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + devThird * S;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + devThird * S;
	
	ptemp[0][0] = 1. - loc_nu[0] * loc_nu[0];
	ptemp[0][1] = - loc_nu[0] * loc_nu[1];
	ptemp[0][2] = - loc_nu[0] * loc_nu[2];
	ptemp[1][1] = 1. - loc_nu[1] * loc_nu[1];
	ptemp[1][0] = - loc_nu[1] * loc_nu[0];
	ptemp[1][2] = - loc_nu[1] * loc_nu[2];
	ptemp[2][2] = 1. - loc_nu[2] * loc_nu[2];
	ptemp[2][0] = - loc_nu[2] * loc_nu[0];
	ptemp[2][1] = - loc_nu[2] * loc_nu[1];

	Qp[0][0] = ptemp[0][0] * Qtemp[0][0] * ptemp[0][0] + ptemp[0][0] * Qtemp[0][1] * ptemp[1][0] + ptemp[0][0] * Qtemp[0][2] * ptemp[2][0]\
				+ ptemp[0][1] * Qtemp[1][0] * ptemp[0][0] + ptemp[0][1] * Qtemp[1][1] * ptemp[1][0] + ptemp[0][1] * Qtemp[1][2] * ptemp[2][0]\
				+ ptemp[0][2] * Qtemp[2][0] * ptemp[0][0] + ptemp[0][2] * Qtemp[2][1] * ptemp[1][0] + ptemp[0][2] * Qtemp[2][2] * ptemp[2][0];

	Qp[0][1] = ptemp[0][0] * Qtemp[0][0] * ptemp[0][1] + ptemp[0][0] * Qtemp[0][1] * ptemp[1][1] + ptemp[0][0] * Qtemp[0][2] * ptemp[2][1]\
				+ ptemp[0][1] * Qtemp[1][0] * ptemp[0][1] + ptemp[0][1] * Qtemp[1][1] * ptemp[1][1] + ptemp[0][1] * Qtemp[1][2] * ptemp[2][1]\
				+ ptemp[0][2] * Qtemp[2][0] * ptemp[0][1] + ptemp[0][2] * Qtemp[2][1] * ptemp[1][1] + ptemp[0][2] * Qtemp[2][2] * ptemp[2][1];

	Qp[0][2] = ptemp[0][0] * Qtemp[0][0] * ptemp[0][2] + ptemp[0][0] * Qtemp[0][1] * ptemp[1][2] + ptemp[0][0] * Qtemp[0][2] * ptemp[2][2]\
				+ ptemp[0][1] * Qtemp[1][0] * ptemp[0][2] + ptemp[0][1] * Qtemp[1][1] * ptemp[1][2] + ptemp[0][1] * Qtemp[1][2] * ptemp[2][2]\
				+ ptemp[0][2] * Qtemp[2][0] * ptemp[0][2] + ptemp[0][2] * Qtemp[2][1] * ptemp[1][2] + ptemp[0][2] * Qtemp[2][2] * ptemp[2][2];

	Qp[1][0] = ptemp[1][0] * Qtemp[0][0] * ptemp[0][0] + ptemp[1][0] * Qtemp[0][1] * ptemp[1][0] + ptemp[1][0] * Qtemp[0][2] * ptemp[2][0]\
				+ ptemp[1][1] * Qtemp[1][0] * ptemp[0][0] + ptemp[1][1] * Qtemp[1][1] * ptemp[1][0] + ptemp[1][1] * Qtemp[1][2] * ptemp[2][0]\
				+ ptemp[1][2] * Qtemp[2][0] * ptemp[0][0] + ptemp[1][2] * Qtemp[2][1] * ptemp[1][0] + ptemp[1][2] * Qtemp[2][2] * ptemp[2][0];

	Qp[1][1] = ptemp[1][0] * Qtemp[0][0] * ptemp[0][1] + ptemp[1][0] * Qtemp[0][1] * ptemp[1][1] + ptemp[1][0] * Qtemp[0][2] * ptemp[2][1]\
				+ ptemp[1][1] * Qtemp[1][0] * ptemp[0][1] + ptemp[1][1] * Qtemp[1][1] * ptemp[1][1] + ptemp[1][1] * Qtemp[1][2] * ptemp[2][1]\
				+ ptemp[1][2] * Qtemp[2][0] * ptemp[0][1] + ptemp[1][2] * Qtemp[2][1] * ptemp[1][1] + ptemp[1][2] * Qtemp[2][2] * ptemp[2][1];

	Qp[1][2] = ptemp[1][0] * Qtemp[0][0] * ptemp[0][2] + ptemp[1][0] * Qtemp[0][1] * ptemp[1][2] + ptemp[1][0] * Qtemp[0][2] * ptemp[2][2]\
				+ ptemp[1][1] * Qtemp[1][0] * ptemp[0][2] + ptemp[1][1] * Qtemp[1][1] * ptemp[1][2] + ptemp[1][1] * Qtemp[1][2] * ptemp[2][2]\
				+ ptemp[1][2] * Qtemp[2][0] * ptemp[0][2] + ptemp[1][2] * Qtemp[2][1] * ptemp[1][2] + ptemp[1][2] * Qtemp[2][2] * ptemp[2][2];

	Qp[2][0] = ptemp[2][0] * Qtemp[0][0] * ptemp[0][0] + ptemp[2][0] * Qtemp[0][1] * ptemp[1][0] + ptemp[2][0] * Qtemp[0][2] * ptemp[2][0]\
				+ ptemp[2][1] * Qtemp[1][0] * ptemp[0][0] + ptemp[2][1] * Qtemp[1][1] * ptemp[1][0] + ptemp[2][1] * Qtemp[1][2] * ptemp[2][0]\
				+ ptemp[2][2] * Qtemp[2][0] * ptemp[0][0] + ptemp[2][2] * Qtemp[2][1] * ptemp[1][0] + ptemp[2][2] * Qtemp[2][2] * ptemp[2][0];

	Qp[2][1] = ptemp[2][0] * Qtemp[0][0] * ptemp[0][1] + ptemp[2][0] * Qtemp[0][1] * ptemp[1][1] + ptemp[2][0] * Qtemp[0][2] * ptemp[2][1]\
				+ ptemp[2][1] * Qtemp[1][0] * ptemp[0][1] + ptemp[2][1] * Qtemp[1][1] * ptemp[1][1] + ptemp[2][1] * Qtemp[1][2] * ptemp[2][1]\
				+ ptemp[2][2] * Qtemp[2][0] * ptemp[0][1] + ptemp[2][2] * Qtemp[2][1] * ptemp[1][1] + ptemp[2][2] * Qtemp[2][2] * ptemp[2][1];

	Qp[2][2] = ptemp[2][0] * Qtemp[0][0] * ptemp[0][2] + ptemp[2][0] * Qtemp[0][1] * ptemp[1][2] + ptemp[2][0] * Qtemp[0][2] * ptemp[2][2]\
				+ ptemp[2][1] * Qtemp[1][0] * ptemp[0][2] + ptemp[2][1] * Qtemp[1][1] * ptemp[1][2] + ptemp[2][1] * Qtemp[1][2] * ptemp[2][2]\
				+ ptemp[2][2] * Qtemp[2][0] * ptemp[0][2] + ptemp[2][2] * Qtemp[2][1] * ptemp[1][2] + ptemp[2][2] * Qtemp[2][2] * ptemp[2][2];
	
	nuQnu = (loc_nu[0]*Qtemp[0][0]*loc_nu[0] + loc_nu[0]*Qtemp[0][1]*loc_nu[1] + loc_nu[0]*Qtemp[0][2]*loc_nu[2]\
			+ loc_nu[1]*Qtemp[1][0]*loc_nu[0] + loc_nu[1]*Qtemp[1][1]*loc_nu[1] + loc_nu[1]*Qtemp[1][2]*loc_nu[2]\
			+ loc_nu[2]*Qtemp[2][0]*loc_nu[0] + loc_nu[2]*Qtemp[2][1]*loc_nu[1] + loc_nu[2]*Qtemp[2][2]*loc_nu[2]) * devThird;

	Qdiff[0] =  Qtemp[0][0]- Qp[0][0] - nuQnu;
	Qdiff[1] =  Qtemp[0][1]- Qp[0][1];
	Qdiff[2] =  Qtemp[0][2]- Qp[0][2];
	Qdiff[3] =  Qtemp[1][1]- Qp[1][1] - nuQnu;
	Qdiff[4] =  Qtemp[1][2]- Qp[1][2];
	Qdiff[5] =  Qtemp[2][2]- Qp[2][2] - nuQnu;
}

__device__ void relax_conic(double Qin[6], double loc_nu[3], double Qdiff[6], double S, double tiltAngle){
	
	double Qtemp[3][3] = {{ 0. }};
	double ptemp[3][3] = {{ 0. }};
	double Qp[3][3] = {{ 0. }};
	double trace = 0.;
	double cosTiltAngleSq = 0.;

	Qtemp[0][0] = Qin[0] + devThird * S;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + devThird * S;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + devThird * S;

	ptemp[0][0] = loc_nu[0] * loc_nu[0];
	ptemp[0][1] = loc_nu[0] * loc_nu[1];
	ptemp[0][2] = loc_nu[0] * loc_nu[2];
	ptemp[1][1] = loc_nu[1] * loc_nu[1];
	ptemp[1][0] = loc_nu[1] * loc_nu[0];
	ptemp[1][2] = loc_nu[1] * loc_nu[2];
	ptemp[2][2] = loc_nu[2] * loc_nu[2];
	ptemp[2][0] = loc_nu[2] * loc_nu[0];
	ptemp[2][1] = loc_nu[2] * loc_nu[1];

	Qp[0][0] = ptemp[0][0] * Qtemp[0][0] * ptemp[0][0] + ptemp[0][0] * Qtemp[0][1] * ptemp[1][0] + ptemp[0][0] * Qtemp[0][2] * ptemp[2][0]\
				+ ptemp[0][1] * Qtemp[1][0] * ptemp[0][0] + ptemp[0][1] * Qtemp[1][1] * ptemp[1][0] + ptemp[0][1] * Qtemp[1][2] * ptemp[2][0]\
				+ ptemp[0][2] * Qtemp[2][0] * ptemp[0][0] + ptemp[0][2] * Qtemp[2][1] * ptemp[1][0] + ptemp[0][2] * Qtemp[2][2] * ptemp[2][0];

	Qp[0][1] = ptemp[0][0] * Qtemp[0][0] * ptemp[0][1] + ptemp[0][0] * Qtemp[0][1] * ptemp[1][1] + ptemp[0][0] * Qtemp[0][2] * ptemp[2][1]\
				+ ptemp[0][1] * Qtemp[1][0] * ptemp[0][1] + ptemp[0][1] * Qtemp[1][1] * ptemp[1][1] + ptemp[0][1] * Qtemp[1][2] * ptemp[2][1]\
				+ ptemp[0][2] * Qtemp[2][0] * ptemp[0][1] + ptemp[0][2] * Qtemp[2][1] * ptemp[1][1] + ptemp[0][2] * Qtemp[2][2] * ptemp[2][1];

	Qp[0][2] = ptemp[0][0] * Qtemp[0][0] * ptemp[0][2] + ptemp[0][0] * Qtemp[0][1] * ptemp[1][2] + ptemp[0][0] * Qtemp[0][2] * ptemp[2][2]\
				+ ptemp[0][1] * Qtemp[1][0] * ptemp[0][2] + ptemp[0][1] * Qtemp[1][1] * ptemp[1][2] + ptemp[0][1] * Qtemp[1][2] * ptemp[2][2]\
				+ ptemp[0][2] * Qtemp[2][0] * ptemp[0][2] + ptemp[0][2] * Qtemp[2][1] * ptemp[1][2] + ptemp[0][2] * Qtemp[2][2] * ptemp[2][2];

	Qp[1][0] = ptemp[1][0] * Qtemp[0][0] * ptemp[0][0] + ptemp[1][0] * Qtemp[0][1] * ptemp[1][0] + ptemp[1][0] * Qtemp[0][2] * ptemp[2][0]\
				+ ptemp[1][1] * Qtemp[1][0] * ptemp[0][0] + ptemp[1][1] * Qtemp[1][1] * ptemp[1][0] + ptemp[1][1] * Qtemp[1][2] * ptemp[2][0]\
				+ ptemp[1][2] * Qtemp[2][0] * ptemp[0][0] + ptemp[1][2] * Qtemp[2][1] * ptemp[1][0] + ptemp[1][2] * Qtemp[2][2] * ptemp[2][0];

	Qp[1][1] = ptemp[1][0] * Qtemp[0][0] * ptemp[0][1] + ptemp[1][0] * Qtemp[0][1] * ptemp[1][1] + ptemp[1][0] * Qtemp[0][2] * ptemp[2][1]\
				+ ptemp[1][1] * Qtemp[1][0] * ptemp[0][1] + ptemp[1][1] * Qtemp[1][1] * ptemp[1][1] + ptemp[1][1] * Qtemp[1][2] * ptemp[2][1]\
				+ ptemp[1][2] * Qtemp[2][0] * ptemp[0][1] + ptemp[1][2] * Qtemp[2][1] * ptemp[1][1] + ptemp[1][2] * Qtemp[2][2] * ptemp[2][1];

	Qp[1][2] = ptemp[1][0] * Qtemp[0][0] * ptemp[0][2] + ptemp[1][0] * Qtemp[0][1] * ptemp[1][2] + ptemp[1][0] * Qtemp[0][2] * ptemp[2][2]\
				+ ptemp[1][1] * Qtemp[1][0] * ptemp[0][2] + ptemp[1][1] * Qtemp[1][1] * ptemp[1][2] + ptemp[1][1] * Qtemp[1][2] * ptemp[2][2]\
				+ ptemp[1][2] * Qtemp[2][0] * ptemp[0][2] + ptemp[1][2] * Qtemp[2][1] * ptemp[1][2] + ptemp[1][2] * Qtemp[2][2] * ptemp[2][2];

	Qp[2][0] = ptemp[2][0] * Qtemp[0][0] * ptemp[0][0] + ptemp[2][0] * Qtemp[0][1] * ptemp[1][0] + ptemp[2][0] * Qtemp[0][2] * ptemp[2][0]\
				+ ptemp[2][1] * Qtemp[1][0] * ptemp[0][0] + ptemp[2][1] * Qtemp[1][1] * ptemp[1][0] + ptemp[2][1] * Qtemp[1][2] * ptemp[2][0]\
				+ ptemp[2][2] * Qtemp[2][0] * ptemp[0][0] + ptemp[2][2] * Qtemp[2][1] * ptemp[1][0] + ptemp[2][2] * Qtemp[2][2] * ptemp[2][0];

	Qp[2][1] = ptemp[2][0] * Qtemp[0][0] * ptemp[0][1] + ptemp[2][0] * Qtemp[0][1] * ptemp[1][1] + ptemp[2][0] * Qtemp[0][2] * ptemp[2][1]\
				+ ptemp[2][1] * Qtemp[1][0] * ptemp[0][1] + ptemp[2][1] * Qtemp[1][1] * ptemp[1][1] + ptemp[2][1] * Qtemp[1][2] * ptemp[2][1]\
				+ ptemp[2][2] * Qtemp[2][0] * ptemp[0][1] + ptemp[2][2] * Qtemp[2][1] * ptemp[1][1] + ptemp[2][2] * Qtemp[2][2] * ptemp[2][1];

	Qp[2][2] = ptemp[2][0] * Qtemp[0][0] * ptemp[0][2] + ptemp[2][0] * Qtemp[0][1] * ptemp[1][2] + ptemp[2][0] * Qtemp[0][2] * ptemp[2][2]\
				+ ptemp[2][1] * Qtemp[1][0] * ptemp[0][2] + ptemp[2][1] * Qtemp[1][1] * ptemp[1][2] + ptemp[2][1] * Qtemp[1][2] * ptemp[2][2]\
				+ ptemp[2][2] * Qtemp[2][0] * ptemp[0][2] + ptemp[2][2] * Qtemp[2][1] * ptemp[1][2] + ptemp[2][2] * Qtemp[2][2] * ptemp[2][2];

	cosTiltAngleSq = pow( cos(tiltAngle / 180.0 * devPI) , 2);

	Qdiff[0] = Qp[0][0] - cosTiltAngleSq * S * ptemp[0][0];
	Qdiff[1] = Qp[0][1] - cosTiltAngleSq * S * ptemp[0][1];
	Qdiff[2] = Qp[0][2] - cosTiltAngleSq * S * ptemp[0][2];
	Qdiff[3] = Qp[1][1] - cosTiltAngleSq * S * ptemp[1][1];
	Qdiff[4] = Qp[1][2] - cosTiltAngleSq * S * ptemp[1][2];
	Qdiff[5] = Qp[2][2] - cosTiltAngleSq * S * ptemp[2][2];

	trace = trace_f(Qdiff);

	Qdiff[0] -= trace;
	Qdiff[3] -= trace;
	Qdiff[5] -= trace;

}