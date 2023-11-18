#pragma once
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include "relaxations.cuh"

//Definiendo variables
extern double dE;
extern double en_tot, old_en;
extern double en_ldg[3];
extern double en_surf[2];
extern double en_el[5];
extern double en_el_in[5];
extern double en_el_out[5];

extern double U;
extern double U2;
extern double W, Wp;
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
extern unsigned int bulk, surf, nsurf;


extern double iRx, iRy, iRz, Rx, Ry, Rz;
extern int rx, ry, rz;
extern int save_every, check_every, stopat;
extern double dx, dy, dz;
extern double idx, idy, idz, iddx, iddy, iddz;
extern int pRx, pRy, pRz, interface;

extern bool* drop;
extern bool* boundary;
extern bool* ndrop;
extern bool* nboundary;

//extern unsigned char* bulktype; We need this variable just locally.
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
extern double tiltAngle;

//host variables
extern unsigned char* h_bulktype;

//extern bool calling_function();
extern bool read_param();
extern bool ellipsoid();
extern bool nanochannel();
extern bool initial();

extern double matr_mult(double* vec);
extern double dir2ten(double* vec, int n, double S);
//Función para iniciar la configuración o seed deseado.
bool conf();
extern double trqq(double* Q);
extern double trqqq(double* Q);
extern bool norm_v(double* vec);
extern int peri(int node, int dir);

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

/* __device__ int chiral_dev;
__device__ __constant__ double L1_dev;
__device__ float L2_dev;
__device__ float L3_dev;
__device__ float L4_dev; */