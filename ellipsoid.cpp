#include"geometry.hpp"

/**
* Ellipsoidal geometry for a distorted droplet in the Z axis.
* For the moment, the inner shell with U and S just suffers a
* deformation in the mentioned axis. It also keeps a delta distance
* in the X and Z axes; just Y changes.
*/

//bool* drop = (bool*)malloc(total_points * sizeof(bool));
//bool* boundary = (bool*)malloc(total_points * sizeof(bool));
//unsigned char* bulktype = (unsigned char*)malloc(total_points * sizeof(unsigned char));
//signed char* index = (signed char*)malloc(total_points * sizeof(signed char));
bool ellipsoid() {

	double x = 0.0, y = 0.0, z = 0.0, dis = 0.0;
	unsigned int l = 0, bulk = 0, innerbulk = 0, outerbulk = 0;
	double dVi = 0.0, dVo = 0.0;
	
	dx = Lx / (double)(Nx - 1);
	dy = Ly / (double)(Ny - 1);
	dz = Lz / (double)(Nz - 1);

	rx = lrint(Nx / 2);
	ry = lrint(Ny / 2);
	rz = lrint(Nz / 2);

	double Rx = Lx / 2. - 2.;
	double Ry = Ly / 2. - 2.;
	double Rz = Lz / 2. - 2.;

	surf = 0;
	
	idx = 1. / dx;
	idy = 1. / dy;
	idz = 1. / dz;

	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

	bulktype = (unsigned char*)malloc(total_points * sizeof(unsigned char));
	drop = (bool*)malloc(total_points * sizeof(bool));
	boundary = (bool*)malloc(total_points * sizeof(bool));
	qindex = (int*)malloc(total_points * sizeof(int));

	//boudaries variables
	int xm = 0, xp = 0, ym = 0, yp = 0, zm = 0, zp = 0;

	for (int i = 0; i < total_points; i++) {
		drop[i] = false;
		boundary[i] = false;
		bulktype[i] = 0;
		qindex[i] = -1;
	}

	if(DoubleU){
		//Defining the droplet. 
		for (int k = 0; k < Nz; k++) {
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {

					//We set dx, dy and dz to 1
					x = (double)(i - rx) * dx;
					y = (double)(j - ry) * dy;
					z = (double)(k - rz) * dz;

					if ( ((x * x) / ((Rx + 0.5) * (Rx + 0.5)) + (y * y) / ((Ry + 0.5) * (Ry + 0.5)) + (z * z) / ((Rz + 0.5) * (Rz + 0.5)) ) <= 1) {
						drop[l] = true;
						bulk++;

						//A�adimos un if anidado para determinar el bulk externo.
						if ( ((x * x) / ((iRx + 0.5) * (iRx + 0.5)) + (y * y) / ((iRy + 0.5) * (iRy + 0.5)) + (z * z) / ((iRz + 0.5) * (iRz + 0.5)) ) <= 1) {
							//for inner bulk
							bulktype[l] = 1; //bulk interno.
							innerbulk++;
						}
						else {
							//for outer bulk
							bulktype[l] = 2; //para bulk externo.
							outerbulk++;
						}

					}
					l++;
				}
			}
		}
	}
	else{
		for (int k = 0; k < Nz; k++) {
			for (int j = 0; j < Ny; j++) {
				for (int i = 0; i < Nx; i++) {

					//We set dx, dy and dz to 1
					x = (double)(i - rx) * dx;
					y = (double)(j - ry) * dy;
					z = (double)(k - rz) * dz;

					if ( ((x * x) / ((Rx + 0.5) * (Rx + 0.5)) + (y * y) / ((Ry + 0.5) * (Ry + 0.5)) + (z * z) / ((Rz + 0.5) * (Rz + 0.5)) ) <= 1) {
						
						drop[l] = true;
						bulk++;

					}
					l++;
				}
			}
		}
	}
	

	droplet = bulk;

	l = 0;
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				if (drop[l]) {
					xm = i - 1 + j * Nx + k * Nx * Ny;
					xp = i + 1 + j * Nx + k * Nx * Ny;
					ym = i + (j - 1) * Nx + k * Nx * Ny;
					yp = i + (j + 1) * Nx + k * Nx * Ny;
					zm = i + j * Nx + (k - 1) * Nx * Ny;
					zp = i + j * Nx + (k + 1) * Nx * Ny;
					if (!drop[xm] || !drop[xp] || !drop[ym] || !drop[yp] || !drop[zm] || !drop[zp]) {
						boundary[l] = true;
						surf++;
					}
				}
				l++;
			}
		}
	}

	bulk -= surf;
	outerbulk -= surf;

	if ((innerbulk + outerbulk) != bulk) {
		printf("Problems with bulk nodes!\n");
		return false;
	}
	else {
		printf("Innerbulk and outerbulk nodes match with total bulk count!\n");
	}
	if ((bulk + surf) != droplet) {
		printf("Problems with droplet nodes!\n");
		return false;
	}
	else {
		printf("Surf and bulk nodes match with total droplet count!\n");
	}

	printf("Droplet nodes number is % d\nBulk nodes number is % d\nSurface nodes number is % d\n", droplet, bulk, surf);

	//Redifining drop nodes for drop array.
	for (int i = 0; i < total_points; i++) {
		if (boundary[i]) drop[i] = false;
	}

	//for the moment, no nanoparticle section is needed. So, let's ignore it.
	dV = ((4.0 / 3.0) * (double)M_PI * (Rx * Ry * Rz)) / (double)bulk;
	dVi = ((4.0 / 3.0) * (double)M_PI * (iRx * iRy * iRz)) / (double)innerbulk;
	dVo = ((4.0 / 3.0) * (double)M_PI * ((Rx) * (Ry) * (Rz)) - (4.0 / 3.0) * (double)M_PI * (iRx * iRy * iRz)) / (double)outerbulk;
	dA = (4.0 * (double)M_PI * pow((pow(Rx * Ry, 1.6075) + pow(Rx * Rz, 1.6075) + pow(Ry * Rz, 1.6075)) / 3.0, 1.0000 / 1.6075)) / (double)surf;

	printf("Internal nodes count = %d\nExternal nodes count = %d\n", innerbulk, outerbulk + surf);
	printf("External bulk count contains surface nodes!!\n");
	printf("External count after substracting surface nodes is %d\n", outerbulk);
	printf("dV = %lf\n", dV);
	printf("dVi = %lf\n", dVi);
	printf("dVo = %lf\n", dVo);
	printf("dA = %lf\n", dA);
	

	//Allocating memory for surface vectors and tensors
	nu = (double*)malloc(surf * 3 * sizeof(double));
	for (int i = 0; i < surf; i++) {
		nu[i * 3 + 0] = 0.0;
		nu[i * 3 + 1] = 0.0;
		nu[i * 3 + 2] = 0.0;
	}

	if (infinite == 0 && degenerate == 0) {
		Qo = (double*)malloc(surf * 6 * sizeof(double));
		for (int i = 0; i < surf; i++) {
			Qo[i * 6 + 0] = 0.0;
			Qo[i * 6 + 1] = 0.0;
			Qo[i * 6 + 2] = 0.0;
			Qo[i * 6 + 3] = 0.0;
			Qo[i * 6 + 4] = 0.0;
			Qo[i * 6 + 5] = 0.0;
		}
	}

	/**Creating array for Qoldand neighbor index.
	* Neighbor's index must contain up, down, left, right, forward and backward neighbor for each node.
	* Signal array determines whether the node is bulk or surface node, or whethere it has undefined neighbors.
	* signal = -1 undefined
	* signal = 1  bulk
	* signal = 2 node bulk with surface
	*/
	Qold = (double*)malloc(droplet * 6 * sizeof(double));
	neighbor = (int*)malloc(droplet * 6 * sizeof(int));
	signal = (signed char*)malloc(droplet * sizeof(signed char));

	for (int i = 0; i < droplet; i++) {
		Qold[i * 6 + 0] = 0.0;
		Qold[i * 6 + 1] = 0.0;
		Qold[i * 6 + 2] = 0.0;
		Qold[i * 6 + 3] = 0.0;
		Qold[i * 6 + 4] = 0.0;
		Qold[i * 6 + 5] = 0.0;

		neighbor[i * 6 + 0] = -1;
		neighbor[i * 6 + 1] = -1;
		neighbor[i * 6 + 2] = -1;
		neighbor[i * 6 + 3] = -1;
		neighbor[i * 6 + 4] = -1;
		neighbor[i * 6 + 5] = -1;

		signal[i] = -1;
	}

	int nb = 0;
	unsigned int nd = 0;

	for (int l = 0; l < total_points; l++) {

		if (drop[l] || boundary[l]) {

			qindex[l] = nd;

			// signal = 0 for bulk
			if (drop[l]) {
				signal[nd] = 0;
			}
			// signal = 2 for surface. signal = 4 for nanosurface node. 
			else if (boundary[l]) {
				signal[nd] = 2;
				nb++;
			}
			nd++;
		}

	}

	if (nd != droplet) {
		printf("Problem in initialization of qtensor. nd is %d not equal to droplet %d.\n", nd, droplet);
		return false;
	}
	if (nb != surf) {
		printf("Problem in initialization of qtensor. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}

	//Initial configuration
	if (!conf()) {
		printf("Initial configuration didn't set!\n");
		return false;
	}
	else {
		printf("Initial configuration set for seed %d.\n", seed);
	}
	l = 0;
	nb = 0;
	nd = 0;
	//float rNx = rint(Nx/2), rNy = rint(Ny / 2), rNz = rint(Nz / 2);
	//defining neighbors and calculating normal vectors for nu.
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {

				nd = qindex[l];

				if (drop[l]) {
					neighbor[nd * 6 + 0] = qindex[i - 1 + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 1] = qindex[i + 1 + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 2] = qindex[i + (j - 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 3] = qindex[i + (j + 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 4] = qindex[i + j * Nx + (k - 1) * Nx * Ny];
					neighbor[nd * 6 + 5] = qindex[i + j * Nx + (k + 1) * Nx * Ny];
				}
				else if (boundary[l]) {
					x = (double)(i - rx) * dx;
					y = (double)(j - ry) * dy;
					z = (double)(k - rz) * dz;
					dis = sqrt(x * x + y * y + z * z);

					if (dis == 0.) {
						printf("Error in neighbors' boundary.\n");
						return false;
					}

					nu[nb * 3 + 0] = -2. * x / (Rx * Rx);
					nu[nb * 3 + 1] = -2. * y / (Ry * Ry);
					nu[nb * 3 + 2] = -2. * z / (Rz * Rz);
					norm_v(&nu[nb * 3]);

					//In this case, we omit the Double U code condition.
					//Surface is S2.

					if(DoubleU){
						if (infinite == 1) {
							Qold[nd * 6 + 0] = dir2ten(&nu[nb * 3], 0, S2);
							Qold[nd * 6 + 1] = dir2ten(&nu[nb * 3], 1, S2);
							Qold[nd * 6 + 2] = dir2ten(&nu[nb * 3], 2, S2);
							Qold[nd * 6 + 3] = dir2ten(&nu[nb * 3], 3, S2);
							Qold[nd * 6 + 4] = dir2ten(&nu[nb * 3], 4, S2);
							Qold[nd * 6 + 5] = dir2ten(&nu[nb * 3], 5, S2);
						}
						else if (degenerate == 0 && infinite == 0) {
							Qo[nb * 6 + 0] = dir2ten(&nu[nb * 3], 0, S2);
							Qo[nb * 6 + 1] = dir2ten(&nu[nb * 3], 1, S2);
							Qo[nb * 6 + 2] = dir2ten(&nu[nb * 3], 2, S2);
							Qo[nb * 6 + 3] = dir2ten(&nu[nb * 3], 3, S2);
							Qo[nb * 6 + 4] = dir2ten(&nu[nb * 3], 4, S2);
							Qo[nb * 6 + 5] = dir2ten(&nu[nb * 3], 5, S2);
						}
					}
					else{
						//infinite, define qtensor and don't evolve any more
						//homeotropic noninfinite, define qo
						if(infinite == 1){
			
							for(int n = 0; n < 6; n ++){
								Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S);
							}	

						}
						else if(degenerate == 0 && infinite == 0){
							for(int n = 0; n < 6; n ++){
								Qo[nb * 6 + n] = dir2ten(&nu[nb * 3], n, S);
							}		
						}
					}
					

					//Define boundaries.
					if (nu[nb * 3 + 0] >= 0) {
						neighbor[nd * 6 + 0] = qindex[i + 1 + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = qindex[i + 2 + j * Nx + k * Nx * Ny];
					}
					else if (nu[nb * 3 + 0] < 0) {
						neighbor[nd * 6 + 0] = qindex[i - 1 + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = qindex[i - 2 + j * Nx + k * Nx * Ny];
					}
					if (nu[nb * 3 + 1] >= 0) {
						neighbor[nd * 6 + 2] = qindex[i + (j + 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = qindex[i + (j + 2) * Nx + k * Nx * Ny];
					}
					else if (nu[nb * 3 + 1] < 0) {
						neighbor[nd * 6 + 2] = qindex[i + (j - 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = qindex[i + (j - 2) * Nx + k * Nx * Ny];
					}
					if (nu[nb * 3 + 2] >= 0) {
						neighbor[nd * 6 + 4] = qindex[i + j * Nx + (k + 1) * Nx * Ny];
						neighbor[nd * 6 + 5] = qindex[i + j * Nx + (k + 2) * Nx * Ny];
					}
					else if (nu[nb * 3 + 2] < 0) {
						neighbor[nd * 6 + 4] = qindex[i + j * Nx + (k - 1) * Nx * Ny];
						neighbor[nd * 6 + 5] = qindex[i + j * Nx + (k - 2) * Nx * Ny];
					}
					nb++;
				}
				l++;
			}
		}
	}

	if (nb != surf) {
		printf("Problem in initialization of share. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}

	int count1;
	for (int nd = 0; nd < droplet; nd++) {
		//for all Bulk point, if one of the neighbor is surface point
		count1 = 0;
		if (signal[nd] == 0) {
			for (int n = 0; n < 6; n++) {
				if (signal[neighbor[nd * 6 + n]] >= 2) {
					count1++;
				}
			}
			if (count1 > 1) {
				signal[nd] += 1;
			}
		}
		//for all surface point, if one of the neighbor is not defined
		else if (signal[nd] < 8 && signal[nd] >= 2) {
			for (int n = 0; n < 6; n++) {
				if (neighbor[nd * 6 + n] == -1) {
					count1++;
				}
			}
			if (count1 > 0) {
				signal[nd] += 1;
			}
		}
		//for all nodes with problem, share +1
	}

	//This part is to resize the typebulk array into the host array.
	h_bulktype = (unsigned char*)malloc(droplet * sizeof(unsigned char));

	nd = 0;
	for (int i = 0; i < total_points; i++) {
		if (bulktype[i] == 1) {
			h_bulktype[nd] = 1;
			nd++;
		}
		else if (bulktype[i] == 2) {
			h_bulktype[nd] = 2;
			nd++;
		}
	}

	//Size must match.
	if (nd != droplet) {
		printf("Error in transfer data to droplet bulktype!\n");
		return false;
	}

	
	printf("Ellipsoid initialized successfully!\n");
	free(qindex);
	free(bulktype);

	return true;
}
