#pragma once
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>

bool ellipsoid();
bool nanochannel();
bool shell();

extern double dir2ten(double* vec, int n, double S);
extern bool norm_v(double *vec);
extern bool conf();
extern int peri(int node, int dir);

/* Array que guarda los nodos que pertenecen al droplet. */
bool* drop;
/* Array que guarda los nodos pertenecientes a la superficie. */
bool* boundary;
/* Array que guarda los nodos pertenecientes a la nanopartícula. */
bool* ndrop;
/* Array que guarda los nodos pertenencientes a la superficie
de la nanopartícula. */
bool* nboundary;
extern int rand_seed;
double* Qo;
double* nu;
//Aquí se guarda el tensor Q.
double* Qold;
signed char* signal;
int* neighbor;
unsigned int droplet, bulk;
unsigned int surf, nsurf;
//Guarda los tipos de nodo para ser transferidos al dispositivo.
unsigned char* h_bulktype;
//Delta de Volumen.
double dV, dVi, dVo, dVshell;
//Delta de Área
double dA;
//Radios o centros del sistema. No de la caja.
double Rx, Ry, Rz;
//Centros de la caja de simulación.
int rx, ry, rz;
int pRx, pRy, pRz;
double pU;
double alpha, beta, gama;
int interface, anchoring;
//Guarda los tipos de nodo de manera temporal antes de ser transferidos al vector del host.
unsigned char* bulktype;

double x_rot, y_rot, z_rot;
double distance;
int nanoparticle_nodes;
double pivotX, pivotY, pivotZ;
int posX, posY, posZ;
int pivotflag;
int count1;
int pdegenerate;

//Deltas
double dx, dy, dz;
//Inverso de Delta
double idx, idy, idz;
//Delta cuadrada
double iddx, iddy, iddz;
double dApart;

//Extern variables
extern double Lx, Ly, Lz;
extern int Nx, Ny, Nz;
extern int total_points;
extern double iRx, iRy, iRz;
extern int degenerate, infinite;
extern int seed;
extern double S, S2;
extern bool DoubleU;
extern double dir1[3];