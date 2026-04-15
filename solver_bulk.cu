#include "solver_kernels.cuh"
#include "solver_common.cuh"

template <bool UseL2, bool UseChiral>
__global__ void relax_bulk(const double* d_Qin_old, double* d_Qout, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, unsigned char* d_Qtensor_signal,
	double U, double U2, int chiral, double qch, double L1, double L2, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt) 
	{
	
	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;


	if (indx < bulk) {
		
		double Qin[6];
		double QQ[6];
		double Qldg[6];
		double trQQ;
		int xm, xp, ym, yp, zm, zp;
		double Qelas[6] = { 0 };

		unsigned int Q_indx = d_Qtensor_index[indx];
		unsigned char Q_signal = d_Qtensor_signal[indx];

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
		Qin[0] = d_Qin_old[Q_indx * 6 + 0];
		Qin[1] = d_Qin_old[Q_indx * 6 + 1];
		Qin[2] = d_Qin_old[Q_indx * 6 + 2];
		Qin[3] = d_Qin_old[Q_indx * 6 + 3];
		Qin[4] = d_Qin_old[Q_indx * 6 + 4];
		Qin[5] = d_Qin_old[Q_indx * 6 + 5];

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

		Qelas[0] = (d_Qin_old[xp * 6 + 0] + d_Qin_old[xm * 6 + 0] - 2. * Qin[0]) * iddx
			+ (d_Qin_old[yp * 6 + 0] + d_Qin_old[ym * 6 + 0] - 2. * Qin[0]) * iddy
			+ (d_Qin_old[zp * 6 + 0] + d_Qin_old[zm * 6 + 0] - 2. * Qin[0]) * iddz;
		Qelas[1] = (d_Qin_old[xp * 6 + 1] + d_Qin_old[xm * 6 + 1] - 2. * Qin[1]) * iddx
			+ (d_Qin_old[yp * 6 + 1] + d_Qin_old[ym * 6 + 1] - 2. * Qin[1]) * iddy
			+ (d_Qin_old[zp * 6 + 1] + d_Qin_old[zm * 6 + 1] - 2. * Qin[1]) * iddz;
		Qelas[2] = (d_Qin_old[xp * 6 + 2] + d_Qin_old[xm * 6 + 2] - 2. * Qin[2]) * iddx
			+ (d_Qin_old[yp * 6 + 2] + d_Qin_old[ym * 6 + 2] - 2. * Qin[2]) * iddy
			+ (d_Qin_old[zp * 6 + 2] + d_Qin_old[zm * 6 + 2] - 2. * Qin[2]) * iddz;
		Qelas[3] = (d_Qin_old[xp * 6 + 3] + d_Qin_old[xm * 6 + 3] - 2. * Qin[3]) * iddx
			+ (d_Qin_old[yp * 6 + 3] + d_Qin_old[ym * 6 + 3] - 2. * Qin[3]) * iddy
			+ (d_Qin_old[zp * 6 + 3] + d_Qin_old[zm * 6 + 3] - 2. * Qin[3]) * iddz;
		Qelas[4] = (d_Qin_old[xp * 6 + 4] + d_Qin_old[xm * 6 + 4] - 2. * Qin[4]) * iddx
			+ (d_Qin_old[yp * 6 + 4] + d_Qin_old[ym * 6 + 4] - 2. * Qin[4]) * iddy
			+ (d_Qin_old[zp * 6 + 4] + d_Qin_old[zm * 6 + 4] - 2. * Qin[4]) * iddz;
		Qelas[5] = (d_Qin_old[xp * 6 + 5] + d_Qin_old[xm * 6 + 5] - 2. * Qin[5]) * iddx
			+ (d_Qin_old[yp * 6 + 5] + d_Qin_old[ym * 6 + 5] - 2. * Qin[5]) * iddy
			+ (d_Qin_old[zp * 6 + 5] + d_Qin_old[zm * 6 + 5] - 2. * Qin[5]) * iddz;

		//printf("L1:%lf L2:%lf L3:%lf L4:%lf", L1_dev, L2_dev, L3_dev, L4_dev);

		double qch0 = 0.;
		double qch1 = 0.;
		double qch2 = 0.;
		double qch3 = 0.;
		double qch4 = 0.;
		double qch5 = 0.;

		if (UseChiral) {
			double dQ[3][6];

				dQ[0][0] = (d_Qin_old[xp * 6 + 0] - d_Qin_old[xm * 6 + 0]) * 0.5 * idx;
			dQ[1][0] = (d_Qin_old[yp * 6 + 0] - d_Qin_old[ym * 6 + 0]) * 0.5 * idy;
			dQ[2][0] = (d_Qin_old[zp * 6 + 0] - d_Qin_old[zm * 6 + 0]) * 0.5 * idz;

			dQ[0][1] = (d_Qin_old[xp * 6 + 1] - d_Qin_old[xm * 6 + 1]) * 0.5 * idx;
			dQ[1][1] = (d_Qin_old[yp * 6 + 1] - d_Qin_old[ym * 6 + 1]) * 0.5 * idy;
			dQ[2][1] = (d_Qin_old[zp * 6 + 1] - d_Qin_old[zm * 6 + 1]) * 0.5 * idz;

			dQ[0][2] = (d_Qin_old[xp * 6 + 2] - d_Qin_old[xm * 6 + 2]) * 0.5 * idx;
			dQ[1][2] = (d_Qin_old[yp * 6 + 2] - d_Qin_old[ym * 6 + 2]) * 0.5 * idy;
			dQ[2][2] = (d_Qin_old[zp * 6 + 2] - d_Qin_old[zm * 6 + 2]) * 0.5 * idz;

			dQ[0][3] = (d_Qin_old[xp * 6 + 3] - d_Qin_old[xm * 6 + 3]) * 0.5 * idx;
			dQ[1][3] = (d_Qin_old[yp * 6 + 3] - d_Qin_old[ym * 6 + 3]) * 0.5 * idy;
			dQ[2][3] = (d_Qin_old[zp * 6 + 3] - d_Qin_old[zm * 6 + 3]) * 0.5 * idz;

			dQ[0][4] = (d_Qin_old[xp * 6 + 4] - d_Qin_old[xm * 6 + 4]) * 0.5 * idx;
			dQ[1][4] = (d_Qin_old[yp * 6 + 4] - d_Qin_old[ym * 6 + 4]) * 0.5 * idy;
			dQ[2][4] = (d_Qin_old[zp * 6 + 4] - d_Qin_old[zm * 6 + 4]) * 0.5 * idz;

			dQ[0][5] = (d_Qin_old[xp * 6 + 5] - d_Qin_old[xm * 6 + 5]) * 0.5 * idx;
			dQ[1][5] = (d_Qin_old[yp * 6 + 5] - d_Qin_old[ym * 6 + 5]) * 0.5 * idy;
			dQ[2][5] = (d_Qin_old[zp * 6 + 5] - d_Qin_old[zm * 6 + 5]) * 0.5 * idz;

			qch0 = 2. * (dQ[1][2] - dQ[2][1]);
			qch3 = 2. * (dQ[2][1] - dQ[0][4]);
			qch5 = 2. * (dQ[0][4] - dQ[1][2]);
			qch1 = dQ[1][4] - dQ[2][3] + dQ[2][0] - dQ[0][2];
			qch2 = dQ[1][5] - dQ[2][4] + dQ[0][1] - dQ[1][0];
			qch4 = dQ[2][2] - dQ[0][5] + dQ[0][3] - dQ[1][1];
		}

		double qelas2_0 = 0.;
		double qelas2_1 = 0.;
		double qelas2_2 = 0.;
		double qelas2_3 = 0.;
		double qelas2_4 = 0.;
		double qelas2_5 = 0.;

		if(UseL2){
			const double dxx0 = (d_Qin_old[xp * 6 + 0] + d_Qin_old[xm * 6 + 0] - 2. * Qin[0]) * iddx;
			const double dxx1 = (d_Qin_old[xp * 6 + 1] + d_Qin_old[xm * 6 + 1] - 2. * Qin[1]) * iddx;
			const double dxx2 = (d_Qin_old[xp * 6 + 2] + d_Qin_old[xm * 6 + 2] - 2. * Qin[2]) * iddx;
			const double dyy1 = (d_Qin_old[yp * 6 + 1] + d_Qin_old[ym * 6 + 1] - 2. * Qin[1]) * iddy;
			const double dyy3 = (d_Qin_old[yp * 6 + 3] + d_Qin_old[ym * 6 + 3] - 2. * Qin[3]) * iddy;
			const double dyy4 = (d_Qin_old[yp * 6 + 4] + d_Qin_old[ym * 6 + 4] - 2. * Qin[4]) * iddy;
			const double dzz2 = (d_Qin_old[zp * 6 + 2] + d_Qin_old[zm * 6 + 2] - 2. * Qin[2]) * iddz;
			const double dzz4 = (d_Qin_old[zp * 6 + 4] + d_Qin_old[zm * 6 + 4] - 2. * Qin[4]) * iddz;
			const double dzz5 = (d_Qin_old[zp * 6 + 5] + d_Qin_old[zm * 6 + 5] - 2. * Qin[5]) * iddz;

			double dxy0 = 0.;
			double dxz0 = 0.;
			double dxy1 = 0.;
			double dxz1 = 0.;
			double dyz1 = 0.;
			double dxy2 = 0.;
			double dxz2 = 0.;
			double dyz2 = 0.;
			double dxy3 = 0.;
			double dxz3 = 0.;
			double dyz3 = 0.;
			double dxy4 = 0.;
			double dxz4 = 0.;
			double dyz4 = 0.;
			double dxy5 = 0.;
			double dxz5 = 0.;
			double dyz5 = 0.;

			if(Q_signal == 0){
				const int xpyp = d_neighbor[xp * 6 + 3];
				const int xpym = d_neighbor[xp * 6 + 2];
				const int xmyp = d_neighbor[xm * 6 + 3];
				const int xmym = d_neighbor[xm * 6 + 2];
				const int xpzp = d_neighbor[xp * 6 + 5];
				const int xpzm = d_neighbor[xp * 6 + 4];
				const int xmzp = d_neighbor[xm * 6 + 5];
				const int xmzm = d_neighbor[xm * 6 + 4];
				const int ypzp = d_neighbor[yp * 6 + 5];
				const int ypzm = d_neighbor[yp * 6 + 4];
				const int ymzp = d_neighbor[ym * 6 + 5];
				const int ymzm = d_neighbor[ym * 6 + 4];

				dxy0 = (d_Qin_old[xpyp * 6 + 0] + d_Qin_old[xmym * 6 + 0] - d_Qin_old[xmyp * 6 + 0] - d_Qin_old[xpym * 6 + 0]) * idx * idy * 0.25;
				dxz0 = (d_Qin_old[xpzp * 6 + 0] + d_Qin_old[xmzm * 6 + 0] - d_Qin_old[xmzp * 6 + 0] - d_Qin_old[xpzm * 6 + 0]) * idx * idz * 0.25;
				dxy1 = (d_Qin_old[xpyp * 6 + 1] + d_Qin_old[xmym * 6 + 1] - d_Qin_old[xmyp * 6 + 1] - d_Qin_old[xpym * 6 + 1]) * idx * idy * 0.25;
				dxz1 = (d_Qin_old[xpzp * 6 + 1] + d_Qin_old[xmzm * 6 + 1] - d_Qin_old[xmzp * 6 + 1] - d_Qin_old[xpzm * 6 + 1]) * idx * idz * 0.25;
				dyz1 = (d_Qin_old[ypzp * 6 + 1] + d_Qin_old[ymzm * 6 + 1] - d_Qin_old[ymzp * 6 + 1] - d_Qin_old[ypzm * 6 + 1]) * idy * idz * 0.25;
				dxy2 = (d_Qin_old[xpyp * 6 + 2] + d_Qin_old[xmym * 6 + 2] - d_Qin_old[xmyp * 6 + 2] - d_Qin_old[xpym * 6 + 2]) * idx * idy * 0.25;
				dxz2 = (d_Qin_old[xpzp * 6 + 2] + d_Qin_old[xmzm * 6 + 2] - d_Qin_old[xmzp * 6 + 2] - d_Qin_old[xpzm * 6 + 2]) * idx * idz * 0.25;
				dyz2 = (d_Qin_old[ypzp * 6 + 2] + d_Qin_old[ymzm * 6 + 2] - d_Qin_old[ymzp * 6 + 2] - d_Qin_old[ypzm * 6 + 2]) * idy * idz * 0.25;
				dxy3 = (d_Qin_old[xpyp * 6 + 3] + d_Qin_old[xmym * 6 + 3] - d_Qin_old[xmyp * 6 + 3] - d_Qin_old[xpym * 6 + 3]) * idx * idy * 0.25;
				dxz3 = (d_Qin_old[xpzp * 6 + 3] + d_Qin_old[xmzm * 6 + 3] - d_Qin_old[xmzp * 6 + 3] - d_Qin_old[xpzm * 6 + 3]) * idx * idz * 0.25;
				dyz3 = (d_Qin_old[ypzp * 6 + 3] + d_Qin_old[ymzm * 6 + 3] - d_Qin_old[ymzp * 6 + 3] - d_Qin_old[ypzm * 6 + 3]) * idy * idz * 0.25;
				dxy4 = (d_Qin_old[xpyp * 6 + 4] + d_Qin_old[xmym * 6 + 4] - d_Qin_old[xmyp * 6 + 4] - d_Qin_old[xpym * 6 + 4]) * idx * idy * 0.25;
				dxz4 = (d_Qin_old[xpzp * 6 + 4] + d_Qin_old[xmzm * 6 + 4] - d_Qin_old[xmzp * 6 + 4] - d_Qin_old[xpzm * 6 + 4]) * idx * idz * 0.25;
				dyz4 = (d_Qin_old[ypzp * 6 + 4] + d_Qin_old[ymzm * 6 + 4] - d_Qin_old[ymzp * 6 + 4] - d_Qin_old[ypzm * 6 + 4]) * idy * idz * 0.25;
				dxy5 = (d_Qin_old[xpyp * 6 + 5] + d_Qin_old[xmym * 6 + 5] - d_Qin_old[xmyp * 6 + 5] - d_Qin_old[xpym * 6 + 5]) * idx * idy * 0.25;
				dxz5 = (d_Qin_old[xpzp * 6 + 5] + d_Qin_old[xmzm * 6 + 5] - d_Qin_old[xmzp * 6 + 5] - d_Qin_old[xpzm * 6 + 5]) * idx * idz * 0.25;
				dyz5 = (d_Qin_old[ypzp * 6 + 5] + d_Qin_old[ymzm * 6 + 5] - d_Qin_old[ymzp * 6 + 5] - d_Qin_old[ypzm * 6 + 5]) * idy * idz * 0.25;
			}

			qelas2_0 = dxx0 + dxy1 + dxz2;
			qelas2_1 = 0.5 * (dxy0 + dxx1 + dxy3 + dyy1 + dxz4 + dyz2);
			qelas2_2 = 0.5 * (dxz0 + dxx2 + dxy4 + dyz1 + dxz5 + dzz2);
			qelas2_3 = dxy1 + dyy3 + dyz4;
			qelas2_4 = 0.5 * (dxy2 + dxz1 + dzz4 + dyz5 + dyy4 + dyz3);
			qelas2_5 = dxz2 + dyz4 + dzz5;
			const double qelas2_trace = devThird * (qelas2_0 + qelas2_3 + qelas2_5);
			qelas2_0 -= qelas2_trace;
			qelas2_3 -= qelas2_trace;
			qelas2_5 -= qelas2_trace;
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

		const double chiral_scale = UseChiral ? (-2. * qch * L1) : 0.;
		d_Qout[Q_indx * 6 + 0] = Qin[0] + dt * (-Qldg[0] + L1 * Qelas[0] + L2 * qelas2_0 + chiral_scale * qch0);
		d_Qout[Q_indx * 6 + 1] = Qin[1] + dt * (-Qldg[1] + L1 * Qelas[1] + L2 * qelas2_1 + chiral_scale * qch1);
		d_Qout[Q_indx * 6 + 2] = Qin[2] + dt * (-Qldg[2] + L1 * Qelas[2] + L2 * qelas2_2 + chiral_scale * qch2);
		d_Qout[Q_indx * 6 + 3] = Qin[3] + dt * (-Qldg[3] + L1 * Qelas[3] + L2 * qelas2_3 + chiral_scale * qch3);
		d_Qout[Q_indx * 6 + 4] = Qin[4] + dt * (-Qldg[4] + L1 * Qelas[4] + L2 * qelas2_4 + chiral_scale * qch4);
		d_Qout[Q_indx * 6 + 5] = Qin[5] + dt * (-Qldg[5] + L1 * Qelas[5] + L2 * qelas2_5 + chiral_scale * qch5);
		
	}
 

}

template __global__ void relax_bulk<false, false>(const double* d_Qin, double* d_Qout, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, unsigned char* d_Qtensor_signal,
	double U, double U2, int chiral, double qch, double L1, double L2, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt);
template __global__ void relax_bulk<false, true>(const double* d_Qin, double* d_Qout, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, unsigned char* d_Qtensor_signal,
	double U, double U2, int chiral, double qch, double L1, double L2, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt);
template __global__ void relax_bulk<true, false>(const double* d_Qin, double* d_Qout, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, unsigned char* d_Qtensor_signal,
	double U, double U2, int chiral, double qch, double L1, double L2, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt);
template __global__ void relax_bulk<true, true>(const double* d_Qin, double* d_Qout, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, unsigned char* d_Qtensor_signal,
	double U, double U2, int chiral, double qch, double L1, double L2, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt);
