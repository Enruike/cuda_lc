#include "functions.hpp"

double matr_mult(double* vec){
	double ans = 0;
	ans = vec[0] * vec[0] +vec[1] * vec[1] +vec[2] * vec[2] + 2 * vec[0] * vec[1] + 2 * vec[0] * vec[2] + 2 * vec[1] * vec[2]; 
	return ans;
}

double dir2ten(double* vec, int n, double S) {
	double third = 1.0 / 3.0;
	switch (n) {
	case 0:
		return S * (vec[0] * vec[0] - third);
	case 1:
		return S * (vec[0] * vec[1]);
	case 2:
		return S * (vec[0] * vec[2]);
	case 3:
		return S * (vec[1] * vec[1] - third);
	case 4:
		return S * (vec[1] * vec[2]);
	case 5:
		return S * (vec[2] * vec[2] - third);
	default:
		printf("Error with dir_to_ten!\n");
		break;
		return 0;
	}
	return 0;
}

double trqq(double* Q) {
	double ans = 0.;
	ans = Q[0] * Q[0] + Q[3] * Q[3] + Q[5] * Q[5]\
		+ 2. * (Q[1] * Q[1] + Q[2] * Q[2] + Q[4] * Q[4]);
	return ans;
}

double trqqq(double* Q) {
	double ans = 0;
	ans = Q[0] * Q[0] * Q[0] + Q[3] * Q[3] * Q[3] + Q[5] * Q[5] * Q[5]\
		+ 6. * Q[1] * Q[2] * Q[4] + 3 * Q[0] * (Q[1] * Q[1] + Q[2] * Q[2])\
		+ 3. * Q[3] * (Q[1] * Q[1] + Q[4] * Q[4]) + 3. * Q[5] * (Q[4] * Q[4] + Q[2] * Q[2]);
	return ans;
}

bool norm_v(double* vec) {
	double mod = 0;
	mod = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	if (mod == 0) {
		printf("Zero vector!\n");
		return false;
	}
	else {
		vec[0] = vec[0] / mod;
		vec[1] = vec[1] / mod;
		vec[2] = vec[2] / mod;
		return true;
	}
}

bool checktr(double* Q){
	double tr = 0;
	double third =  1.0 / 3.0;
	int n;
	tr = (Q[0] + Q[3] + Q[5]) * third;
	if(tr > 1e-5) {
		printf("%f %f %f %f %f %f\n", Q[0], Q[1], Q[2], Q[3], Q[4], Q[5]);
		Q[0] -= tr;
		Q[3] -= tr;
		Q[5] -= tr;
//		printf("Non-tracelss.\n");
		return false;				
	}
	if(trqq(Q) > 1){
//		for(n = 0; n < 6; n ++){
//			Q[n] /= 1.3;
//		}
		printf("Order parameter exceed 1.\n");
		return false;
	}
	return true;
}

int peri(int node, int dir){
	if(dir == 0){
		if(node >= 0 && node < Nx)	
			return node;
		else if(node >= Nx)
			return node - Nx;
		else
			return node + Nx;
	}
	else if(dir == 1){
		if(node >= 0 && node < Ny)	
			return node;
		else if(node >= Ny)
			return node - Ny;
		else 
			return node + Ny;
	}
	else if(dir == 2){
		if(node >= 0 && node < Nz)	
			return node;
		else if(node >= Nz)
			return node - Nz;
		else 
			return node + Nz;
	}
	else{
		printf("Wrong input for function periodic.\n");
	}
	return 0;
}

void output(){
	FILE* energy;
	FILE* Qtensor;
	int l, indx;

	//print energy
	energy = fopen("energy.out", "a");
	fprintf(energy,"%d\t%.9lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", cycle, dE, en_ldg[0], en_el[0], en_el[1], en_el[2], en_el[3], en_el[4], en_surf[0], en_surf[1], en_tot);
	fclose(energy);

	//Print Separated Energy.
	if(DoubleU){
		FILE* energy2;
		energy2 = fopen("separated_energy.out", "a");
		fprintf(energy2,"%d\t%.9lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", cycle, dE, en_ldg[0], en_ldg[1], en_ldg[2], en_el[0], en_el_in[0], en_el_out[0], en_el[4], en_el_in[4], en_el_out[4], en_surf[0], en_tot);
		fclose(energy2);
	}

	//print Qtensor
	Qtensor = fopen("Qtensor.bin", "wb");
	indx = 0;
	for(int l = 0; l < droplet; l++){
			fwrite(&Qold[l * 6], sizeof(double), 6, Qtensor);
	}
	fclose(Qtensor);
}