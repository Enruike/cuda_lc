#pragma once
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

//Definiendo variables
extern double dE;
extern double en_tot, old_en;
extern double en_ldg;
extern double en_surf[2];
extern double en_el[5];

extern double U, U2, W, Wp;
extern double tmin, tmax, dtime, dt;
extern double increment;
extern int total_points, Nx, Ny, Nz;
extern double Lx, Ly, Lz;
extern int cycle;
extern bool flag;
extern double S, S2;
extern double dV, dVi, dVo, dA, dApart;

extern int Np;
extern double Rp;
extern int rand_seed;

extern int uppersurf, lowersurf;
extern double dir1[3], dir2[3];

extern int surfdegen, ideal, DoubleU;

//Aproximaci�n a una constante el�stica L1
extern double L1, L2, L3, L4, qch;
extern double redshift;
extern double accuracy;
extern int chiral, geo, seed;
extern int degenerate, infinite;
extern double init_dir[3];
extern unsigned int droplet;
extern unsigned int surf;

extern double iRx, iRy, iRz, Rx, Ry, Rz;
extern int rx, ry, rz;
extern int save_every, check_every, stopat;
extern double dx, dy, dz;
extern double idx, idy, idz, iddx, iddy, iddz;

//Alojamiento de arrays.
//Arrays alocations.
extern bool* drop;
extern bool* boundary;
//extern unsigned char* bulktype; We need this variable just locally.
extern int* qindex;
extern double* nu;
extern double* Qo;
extern double* Qold;
extern double* Q_new;
extern signed char* signal;
extern int* neighbor;

extern unsigned char* bulktype;

//device variables
extern double* d_Qold;
extern unsigned char* d_bulktype;
extern signed char* d_signal;

extern double* d_nu;

extern int cycle;

//host variables
extern unsigned char* h_bulktype;

//extern bool calling_function();
extern bool read_param();
extern bool ellipsoid();
extern bool initial();

//Other Functions
extern double dir2ten(double* vec, int n, double S);
extern bool conf();
extern double trqq(double Qin[6]);
extern double trqqq(double Q[6]);
extern bool norm_v(double* vec);

__device__ double d_trqq(double Qin[6]);
//extern __device__ __constant__ double devThird;

//Energy functions
extern void free_energy();
extern double ldg_energy(double* Qold);
extern void elastic_energy(double ans[5]);
extern void surface_energy(double ans[2]);

extern bool checktr(double* Q);
extern void output();

// void relax_surface();
// void device_relax_bulk();

//Relaxation functions.
__global__ void relax_bulk(double* d_Qold, unsigned char* d_bulktype, signed int* d_neighbor,
	unsigned int* d_Qtensor_index, int chiral, double U, double U2, double qch, int L1, unsigned int bulk, double idx, double idy, double idz,
	double iddx, double iddy, double iddz, double dt);

__global__ void relax_surf(double* d_Qold, signed int* d_neighbor, unsigned int* d_Nvector_index, unsigned char* d_Nvector_signal, double* d_Qo,
	int chiral, double qch, int L1, unsigned int surf, int degenerate, int infinite, double W, double Wp, double* d_nu, double d_idx, double d_idy, double d_idz, double dt);
