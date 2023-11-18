#ifndef __RELAXAIONS_CUH__
#define __RELAXAIONS_CUH__

#include<stdio.h>

//Relaxation functions.
__global__ void relax_bulk(double* d_Qold, unsigned char* d_bulktype, signed int* d_neighbor, unsigned int* d_Qtensor_index, unsigned char* d_Qtensor_signal, double U,
    double U2, int chiral, double qch, double L1, double L2, unsigned int bulk, double idx, double idy, double idz, double iddx, double iddy, double iddz, double dt);

__global__ void relax_surf(double* d_Qold, signed int* d_neighbor, unsigned int* d_Nvector_index, unsigned char* d_Nvector_signal, double* d_Qo,
	int chiral, double qch, double L1, double L2, double L3, double L4, double tiltAngle, unsigned int surf, int degenerate, int infinite, double W, double Wp, 
	double* d_nu, double d_idx, double d_idy, double d_idz, double dt, double S);

/* Device Functions */
__device__ void relax_degen(double Qin[6], double loc_nu[3], double Qdiff[6], double S);
__device__ void relax_conic(double Qin[6], double loc_nu[3], double Qdiff[6], double S, double tiltAngle);

#endif