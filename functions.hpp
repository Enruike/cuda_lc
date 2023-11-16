#ifndef __FUNCTIONS_HPP__
#define __FUNCTIONS_HPP__

#include<iostream>
#include<stdlib.h>
#include<math.h>

extern int DoubleU, cycle, droplet;
extern double* Qold;
extern double dE, en_tot;
extern int Nx, Ny, Nz;
extern double en_ldg[3], en_el[5], en_surf[2], en_el_in[5], en_el_out[5];

double matr_mult(double* vec);
//change director to Qtensor
double dir2ten(double* vec, int n, double S);

double trqq(double* Q);

double trqqq(double* Q);

bool norm_v(double* vec);

bool checktr(double* Q);

int peri(int node, int dir);

void output();

#endif