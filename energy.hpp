#include<stdio.h>

extern double* Qold;
extern signed char* signal;
void free_energy();
void ldg_energy(double ans[3]);
void elastic_energy(double ans[5], double ans_in[5], double ans_out[5]);
void surface_energy(double ans[2]);
void en_degen(double* Qin, double* loc_nu, double* Qdiff);
void en_conic(double* Qin, double* loc_nu, double* Qdiff);

//Definiendo variables
double dE;
double en_tot, old_en;
/* Guarda la energía de Landau - de Gennes para región total, interna y externa.
En ese preciso orden. */
double en_ldg[3];
/* Guarda la energía de superficie. Únicamente hace referencia a la energía del
sistema y a la de las nanopartículas en ello. */
double en_surf[2];
/* Almacena la energía elástica referente al primer tipo de U.
Las energías almacenadas son L1, L2, L3, L4 y Quiral. */
double en_el[5];
/* Energía elástica interna */
double en_el_in[5];
/* Energía elástica externa */
double en_el_out[5];
extern double S0;

extern int cycle, check_every, droplet;
extern double trqq(double* Q);
extern double trqqq(double* Q);
extern double matr_mult(double* vec);
extern unsigned char* h_bulktype;
extern int* neighbor;
extern double U, U2, idx, idy, idz;
extern int chiral;
extern double dA, W, Wp, dApart, qch;
//Constantes  Elásticas.
extern double L1, L2;
/* Delta de volúmenes para las distintas regiones.
dV = volumen total
dVi = volumen interno o principal
dVo  = volumen externo o secundario */
extern double dV, dVi, dVo;
extern double tiltAngle;

extern int geo;
extern int DoubleU;
extern int degenerate, infinite;
extern double* nu;
extern double* Qo;