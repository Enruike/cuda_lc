#include"energy.hpp"

void free_energy() {


	ldg_energy(en_ldg);
	elastic_energy(en_el, en_el_in, en_el_out);
	surface_energy(en_surf);

	en_tot = en_ldg[0] + en_el[0] + en_el[1] + en_el[2] + en_el[3] + en_el[4] + en_surf[0] + en_surf[1];
	dE = en_tot - old_en;
	old_en = en_tot;

	if(DoubleU){
			if(cycle % check_every == 0){
			printf("LdG: %lf, En_L1: %lf", en_ldg[0], en_el[0]);
			if(L2 != 0){
				printf(", L2: %lf", en_el[1]);
			}
			if(chiral != 0){
				printf(", Chiral: %lf", en_el[4]);
			}
			printf(", En_Surf1: %lf, Cycle: %d", en_surf[0], cycle);
			if(geo != 10){
				printf("\n\t\t ~Inner Region Energy~\n");
			}
			else{
				printf("\n\t\t ~Energy~\n");
			}
			printf("LdG_in: %lf, L1: %lf", en_ldg[1], en_el_in[0]);
			if(L2 != 0){
				printf(", L2: %lf", en_el_in[1]);
			}
			if(chiral != 0){
				printf(", Chiral: %lf", en_el_in[4]);
			}
			if(geo != 10){
				printf("\n\t\t ~Outer Region Energy~\n");
			}
			else{
				printf("\n\t ~Energy Around the Interface~\n");
			}
			printf("LdG_out: %lf, L1: %lf", en_ldg[2], en_el_out[0]);
			if(L2 != 0){
				printf(", L2: %lf", en_el_out[1]);
			}
			if(chiral != 0){
				printf(", Chiral: %lf", en_el_out[4]);
			}
			printf("\nTotal Energy: %lf \n", en_tot);		
			printf("dE: %lf \n\n", dE);		
		}
	}
	else{

		if (cycle % check_every == 0) {
			printf("En_LDG: %lf, En_L1: %lf, ", en_ldg[0], en_el[0]);
			if(en_el[1] != 0){
				printf("En_L2: %lf, ", en_el[1]);
			}
			if(en_el[2] != 0){
				printf("En_L3: %lf, ", en_el[2]);
			}
			if(en_el[3] != 0){
				printf("En_L4: %lf, ", en_el[3]);
			}
			if(chiral != 0){
				printf("Chiral: %lf\n", en_el[4]);
			}
			printf("En_Surf1: %lf, ", en_surf[0]);
			if(en_surf[1] != 0){
				printf("En_Surf2: %lf ", en_surf[1]);
			}
			printf("Cycle: %d\n", cycle);
			printf("Total energy is: %lf \n", en_tot);
			printf("dE: %lf \n", dE);
		}

	}

}

void ldg_energy(double ldg_ans[3]) {

	double Qin[6] = { 0. };
	ldg_ans[3] = { 0. };
	double trace2 = 0.;
	double trace3 = 0.;

	if(DoubleU){
		for (int i = 0; i < droplet; i++) {
			if (signal[i] == 0 || signal[i] == 1) {

				Qin[0] = Qold[i * 6 + 0];
				Qin[1] = Qold[i * 6 + 1];
				Qin[2] = Qold[i * 6 + 2];
				Qin[3] = Qold[i * 6 + 3];
				Qin[4] = Qold[i * 6 + 4];
				Qin[5] = Qold[i * 6 + 5];

				//for(int j = 0; j < 6; j++){
				//	printf("value Q1 is %lf\n",  Qin[j]);
				//}

				trace2 = trqq(Qin);
				trace3 = trqqq(Qin);

				if (h_bulktype[i] == 1) {
					ldg_ans[0] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
					ldg_ans[1] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
				}
				else if (h_bulktype[i] == 2 || h_bulktype[i] == 3) {
					ldg_ans[0] += 0.5 * (1. - U2 / 3.) * trace2 - U2 / 3. * trace3 + U2 * 0.25 * trace2 * trace2;
					ldg_ans[2] += 0.5 * (1. - U2 / 3.) * trace2 - U2 / 3. * trace3 + U2 * 0.25 * trace2 * trace2;
				}
			}
		}

		ldg_ans[0] *= dV;
		ldg_ans[1] *= dVi;
		ldg_ans[2] *= dVo;
	}
	else{
		for (int i = 0; i < droplet; i++) {
			if (signal[i] == 0 || signal[i] == 1) {

				Qin[0] = Qold[i * 6 + 0];
				Qin[1] = Qold[i * 6 + 1];
				Qin[2] = Qold[i * 6 + 2];
				Qin[3] = Qold[i * 6 + 3];
				Qin[4] = Qold[i * 6 + 4];
				Qin[5] = Qold[i * 6 + 5];

				//for(int j = 0; j < 6; j++){
				//	printf("value Q1 is %lf\n",  Qin[j]);
				//}

				trace2 = trqq(Qin);
				trace3 = trqqq(Qin);
			
				ldg_ans[0] += 0.5 * (1. - U2 / 3.) * trace2 - U2 / 3. * trace3 + U2 * 0.25 * trace2 * trace2;
				
			}
		}

		ldg_ans[0] *= dV;

	}
}

void elastic_energy(double ans[5], double ans_in[5], double ans_out[5]) {
	
	double dQ[3][6] = {{0.}};
	double Qin[6] = { 0 };
	double vec[3] = { 0. };
	ans[5] = { 0. };
	ans_in[5] = { 0. };
	ans_out[5] = { 0. };

	int xm, xp, ym, yp, zm, zp;

	for (int i = 0; i < droplet; i++) {
		if (signal[i] == 0 || signal[i] == 1) {
			for (int n = 0; n < 6; n++)	Qin[n] = Qold[i * 6 + n];
			xm = neighbor[i * 6 + 0];
			xp = neighbor[i * 6 + 1];
			ym = neighbor[i * 6 + 2];
			yp = neighbor[i * 6 + 3];
			zm = neighbor[i * 6 + 4];
			zp = neighbor[i * 6 + 5];
			for (int n = 0; n < 6; n++) {
				//dQ is the first derivative with second order approximation
				//first index for direction: 0-x; 1-y; 2-z;
				//second index for qtensor index;
				dQ[0][n] = (Qold[xp * 6 + n] - Qold[xm * 6 + n]) * 0.5 * idx;
				dQ[1][n] = (Qold[yp * 6 + n] - Qold[ym * 6 + n]) * 0.5 * idy;
				dQ[2][n] = (Qold[zp * 6 + n] - Qold[zm * 6 + n]) * 0.5 * idz;
			}
			ans[0] += trqq(dQ[0]) + trqq(dQ[1]) + trqq(dQ[2]);
			if(DoubleU){
				if(h_bulktype[i] == 1){
					ans_in[0] += trqq(dQ[0])+trqq(dQ[1])+trqq(dQ[2]);
				}
				else if(h_bulktype[i] == 2 || h_bulktype[i] == 3){
					ans_out[0] += trqq(dQ[0])+trqq(dQ[1])+trqq(dQ[2]);
				}
			}

			if (L2 != 0){
				vec[0] = dQ[0][0];
				vec[1] = dQ[1][1];
				vec[2] = dQ[2][2];
				ans[1] += matr_mult(vec);

				if(DoubleU){
					if(h_bulktype[i] == 1){
						ans_in[1] += matr_mult(vec);
					}
					else if(h_bulktype[i] == 2 || h_bulktype[i] == 3){
						ans_out[1] += matr_mult(vec);
					}
				}

				vec[0] = dQ[0][1];
				vec[1] = dQ[1][3];
				vec[2] = dQ[2][4];
				ans[1] += matr_mult(vec);

				if(DoubleU){
					if(h_bulktype[i] == 1){
						ans_in[1] += matr_mult(vec);
					}
					else if(h_bulktype[i] == 2 || h_bulktype[i] == 3){
						ans_out[1] += matr_mult(vec);
					}
				}

				vec[0] = dQ[0][2];
				vec[1] = dQ[1][4];
				vec[2] = dQ[2][5];
				ans[1] += matr_mult(vec);

				if(DoubleU){
					if(h_bulktype[i] == 1){
						ans_in[1] += matr_mult(vec);
					}
					else if(h_bulktype[i] == 2 || h_bulktype[i] == 3){
						ans_out[1] += matr_mult(vec);
					}
				}
				
			}
			
			if (chiral == 1) {
				//Chiral elastic energy
				ans[4] += Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] + Qin[1] * dQ[2][0] + Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2] + Qin[2] * dQ[0][1]\
					+ Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4] - Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4] - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1]\
					- Qin[5] * dQ[1][2] - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];

				if(DoubleU){
					if(h_bulktype[i] == 1){
						ans_in[4] += Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] + Qin[1] * dQ[2][0] + Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2]  + Qin[2] * dQ[0][1] + Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4]  - Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4]  - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1] - Qin[5] * dQ[1][2] - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];					}
					else if(h_bulktype[i] == 2 || h_bulktype[i] == 3){
						ans_out[4] += Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] + Qin[1] * dQ[2][0] + Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2]  + Qin[2] * dQ[0][1] + Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4]  - Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4]  - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1] - Qin[5] * dQ[1][2] - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];					}
				}
			}
		}
	}

	if(DoubleU){
		ans_in[0] *= 0.5 * dVi * L1;
		ans_in[1] *= 0.5 * dVi * L2;
		ans_in[4] *= dVi * (double)chiral * 2. * L1 * qch;

		ans_out[0] *= 0.5 * dVo * L1;
		ans_out[1] *= 0.5 * dVo * L2;
		ans_out[4] *= dVo * (double)chiral * 2. * L1 * qch;
	}
	ans[0] *= 0.5 * dV * L1;
	ans[1] *= 0.5 * dV * L2;
	ans[4] *= dV * (double)chiral * 2. * L1 * qch;
}

void surface_energy(double ans[2]) {
	int nb = 0;
	double Qdiff[6] = { 0 };
	double Qin[6] = { 0 };
	double loc_nu[3] = { 0 };
	int degen = 0, inf = 1;
	double Wstr = 0.;
	bool npboundary = true;

	for (int i = 0; i < droplet; i++) {
		if (signal[i] >= 2 && signal[i] <= 8 || signal[i] == 12 || signal[i] == 13) {
			//for channel boundary
			if (signal[i] == 2 || signal[i] == 3) {
				degen = degenerate;
				inf = infinite;
				Wstr = W;
				dA = dA;
				npboundary = false;
			}
			//Para superficie inferior degenerada.
			else if (signal[i] == 12 || signal[i] == 13) {
				degen = 1;
				inf = infinite;
				Wstr = W;
				dA = dA;
				npboundary = false;
			}

			//for nanoparticle boundary
			else if (signal[i] == 4 || signal[i] == 5) {
				degen = 0;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if (signal[i] == 6 || signal[i] == 7) {
				degen = 1;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if (signal[i] == 8) {
				degen = 0;
				inf = 1;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			// else {
			// 	printf("Error in energy_surf.\n");
			// }

			if (Wstr != 0 && inf != 1) {
				//printf("Test. sign[i] = %d\n", sign[i]);
				if (degen == 1) {
					for (int n = 0; n < 6; n++)	Qin[n] = Qold[i * 6 + n];
					for (int n = 0; n < 3; n++)	loc_nu[n] = nu[nb * 3 + n];
					//en_degen(Qin, loc_nu, Qdiff);

					if (npboundary) {
						ans[1] += Wstr * trqq(Qdiff) * dApart;
					}
					else {
						ans[0] += Wstr * trqq(Qdiff) * dA;
					}

				}
				else if (degen == 0 && inf == 0) {
					for (int n = 0; n < 6; n++) {
						Qdiff[n] = Qold[i * 6 + n] - Qo[nb * 6 + n];
					}
					if (npboundary) {
						ans[1] += Wstr * trqq(Qdiff) * dApart;
					}
					else {
						ans[0] += Wstr * trqq(Qdiff) * dA;
					}
				}
			}
			nb++;
		}
	}
}