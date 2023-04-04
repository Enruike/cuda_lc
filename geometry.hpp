#include<stdio.h>
#include<math.h>
#include<iostream>


bool ellipsoid();
extern double dir2ten(double* vec, int n, double S);
extern bool norm_v(double *vec);
extern bool conf();

bool* drop;
bool* boundary;
int* qindex;
double* Qo;
double* nu;
double* Qold;
signed char* signal;
int* neighbor;
unsigned int droplet;
unsigned char* h_bulktype;
double dV = 0.0;
double Rx, Ry, Rz;
int rx, ry, rz;
unsigned char* bulktype;

double dx, dy, dz;
double idx, idy, idz;
double iddx, iddy, iddz;

unsigned int surf;
double dA = 0.0;
double dApart = 0.0;

//Extern variables
extern double Lx, Ly, Lz;
extern int Nx, Ny, Nz;
extern int total_points;
extern double iRx, iRy, iRz;
extern int degenerate, infinite;
extern int seed;
extern double S, S2;
