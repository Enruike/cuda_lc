#include<iostream>
#include<stdio.h>

bool read_param();

extern bool norm_v(double *vec);

int Nx, Ny, Nz, Np;
double Lx, Ly, Lz;
double Wp, Rp;
int save_every, check_every, stopat;
double W, qch, L1, L2, L3, L4;
//Región 1 o interna.
double U;
//Región 2 o externa.
double U2;
int degenerate, infinite;
int geo, seed, chiral;
double tmin, tmax , increment, accuracy;
double redshift;
double init_dir[3];
double dir1[3], dir2[3];
int DoubleU;
int uppersurf, lowersurf;
int rand_seed;
int surfdegen;
double iRx, iRy, iRz;