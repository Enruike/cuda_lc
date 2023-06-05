#include"geometry.hpp"

bool read_nano(){

	FILE* param = fopen("nano.in", "r");

	if (param == (FILE*)NULL) {
		printf("No nano.in file found!\n");
		return false;
	}

    fscanf(param, "pRx %d #size of nanoparticle\n", &pRx);
    fscanf(param, "pRy %d\n", &pRy);
    fscanf(param, "pRz %d\n", &pRz);
    fscanf(param, "pU %lf #U for the interface layer\n", &pU);
    fscanf(param, "alpha %lf #Angles for rotate or tilt the nano particle\n", &alpha);
    fscanf(param, "beta %lf\n", &beta);
    fscanf(param, "gamma %lf\n", &gama);
    fscanf(param, "interface %d #thickness of the interface layer; 0: no interface\n", &interface);
    fscanf(param, "anchoring %d #0:random 1:homeotropic 2:planar\n", &anchoring);

    printf("~ Nanoparticle data ~\n");
    printf("pRx %d\n", pRx);
    printf("pRy %d\n", pRy);
    printf("pRz %d\n", pRz);
    printf("pU %lf\n", pU);
    printf("alpha: %lf beta: %lf gamma: %lf\n", alpha, beta, gama);
    printf("interface nodes %d\n", interface);
    if(anchoring == 0){
        printf("random anchoring\n");
    }
    else if(anchoring == 1){
        printf("homeotropic anchoring\n");
    }
    else if(anchoring == 2){
        printf("planar anchoring\n");
    }
    else{
        printf("unknonw anchoring\n");
        return false;
    }


    return true;

}

bool nanochannel(){

	//Aquí guardo los índices de los nodos.
	int* qindex = (int*)malloc(total_points * sizeof(int));

    double x = 0.0, y = 0.0, z = 0.0, dis = 0.0;
	unsigned int l = 0;
	double x_rot, y_rot, z_rot;
    double distance;
	
	dx = Lx / (double)(Nx - 1);
	dy = Ly / (double)(Ny - 1);
	dz = Lz / (double)(Nz - 1);

	rx = lrint(Nx / 2);
	ry = lrint(Ny / 2);
	rz = lrint(Nz / 2);
	
	Rx = Lx / 2. - 2.;
	Ry = Ly / 2. - 2.;
	Rz = Lz / 2. - 2.;

	surf = 0;
    bulk = Nx * Ny * Nz;
	
	idx = 1 / dx;
	idy = 1 / dy;
	idz = 1 / dz;

	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

	bulktype = (unsigned char*)malloc(total_points * sizeof(unsigned char));
	drop = (bool*)malloc(total_points * sizeof(bool));
	boundary = (bool*)malloc(total_points * sizeof(bool));
    ndrop = (bool*)malloc(total_points * sizeof(bool));
    nboundary = (bool*)malloc(total_points * sizeof(bool));

    //boudaries variables
	int xm = 0, xp = 0, ym = 0, yp = 0, zm = 0, zp = 0;

	for (int i = 0; i < total_points; i++) {
		drop[i] = true;
		boundary[i] = false;
        ndrop[i] = false;
        nboundary[i] = false;
		bulktype[i] = 1;
		qindex[i] = -1;
	}

    //Definir superficie del nanocanal
    for(int k = 0; k < Nz; k += (Nz - 1)){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){
                l = i + j * Nx + k * Nx * Ny;
                boundary[l] = true;
                drop[l] = false;
                bulktype[l] = 4;
                surf++;
                bulk--;
            }
        }
    }

    //Leer valores del archivo nano.in
    if(!read_nano()){
        return false;
    }

    alpha = (alpha * M_PI) / 180.;
    beta = (beta * M_PI) / 180.;
    gama = (gama * M_PI) / 180.;

    /*
    Nano particle is in the center of the channel.
    'til now, there's no need of defining a position.
    */

    l = 0;
    int nanoparticle_nodes = 0;

    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

				x = (double)(i - rx) * dx;
				y = (double)(j - ry) * dy;
				z = (double)(k - rz) * dz;

                x_rot = x * cos(alpha) * cos(beta) + y * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
					+ z * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
				y_rot = x * sin(alpha) * cos(beta) + y * (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
					+ z * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
				z_rot = x * -sin(beta) + y * cos(beta) * sin(gama) + z * cos(beta) * cos(gama);

				x = x_rot;
				y = y_rot;
				z = z_rot;

                distance = (x * x) / ((pRx + 0.5) * (pRx + 0.5))\
                    + (y * y) / ((pRy + 0.5) * (pRy + 0.5))\
                    + (z * z) / ((pRz + 0.5) * (pRz + 0.5));
                
                if(distance <= 1){
                    ndrop[l] = true;
                    nanoparticle_nodes++;
                    drop[l] = false;
                    bulktype[l] = 5;
                    bulk--;
                }

                l++;
            
            }
        }
    }

    //Defining Nano boundaries
    /*
    La superficie de la nanopartícula tiene un anclaje infinito.
    Pero la superficie del canal no necesariamente debe tener un
    anclaje infinito también.
    */

    l = 0;

    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                if(drop[l]){
                    xm = peri(i - 1, 0) + j * Nx + k * Nx * Ny;
                    xp = peri(i + 1, 0) + j * Nx + k * Nx * Ny;
                    ym = i + peri(j - 1, 1) * Nx + k * Nx * Ny;
                    yp = i + peri(j + 1, 1) * Nx + k * Nx * Ny;
                    zm = i + j * Nx + (k - 1) * Nx * Ny;
                    zp = i + j * Nx + (k + 1) * Nx * Ny;
                    
                    if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zm] || ndrop[zp]){
						nboundary[l] = true;
						drop[l] = false;
                        bulktype[l] = 6;
						bulk--;
						nsurf++;
					}
                }
                l++;          
            }
        }
    }

    /*
    Definiendo bordes de cristal líquido que no evolucionarán.
    Todo dependerá del grosor de la capa (variable interface)
    */

    l = 0;
    int interbulk = 0;

    if(interface != 0){
       
        int rebulker = bulk;

        
        for(int k = 0; k < Nz; k++){
            for(int j = 0; j < Ny; j++){
                for(int i = 0; i < Nx; i++){
                    if(nboundary[l]){

                        for(int node = 0; node < interface; node++){
                            xm = peri(i - (node + 1), 0) + j * Nx + k * Nx * Ny;
                            xp = peri(i + (node + 1), 0) + j * Nx + k * Nx * Ny;
                            ym = i + peri(j - (node + 1), 1) * Nx + k * Nx * Ny;
                            yp = i + peri(j + (node + 1), 1) * Nx + k * Nx * Ny;
                            zm = i + j * Nx + peri(k - (node + 1), 2) * Nx * Ny;
                            zp = i + j * Nx + peri(k + (node + 1), 2) * Nx * Ny;
                            
                            
                            if(drop[xm] && bulktype[xm] != 3){
                                bulktype[xm] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[xp] && bulktype[xp] != 3){
                                bulktype[xp] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[ym] && bulktype[ym] != 3){
                                bulktype[ym] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[yp] && bulktype[yp] != 3){
                                bulktype[yp] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[zm] && bulktype[zm] != 3){
                                bulktype[zm] = 3;
                                interbulk++;
                                rebulker--;
                            }
                            if(drop[zp] && bulktype[zp] != 3){
                                bulktype[zp] = 3;
                                interbulk++;
                                rebulker--;
                            }
                        }
                    }
                    l++;
                }
            }
        }

        if(bulk != (rebulker + interbulk++)){
            printf("Problems with interface nodes\n");
            return false;
        }

    }

    dV = (Lx * Ly * Lz - 4. / 3. * M_PI * pRx * pRy * pRz) / bulk;
    dVi = (Lx * Ly * Lz - 4. / 3. * M_PI * pRx * pRy * pRz) / (bulk - interbulk);
    dVo = (4 / 3 * M_PI * ((pRx + interface) * (pRy + interface) * (pRz + interface) - (pRx) * (pRy) * (pRz))) / interbulk;
    dA = (2 * Lx * Ly) / (surf);
    dApart = 4. * M_PI * pow((pow(pRx * pRy, 1.6075) + pow(pRx * pRz, 1.6075) + pow(pRy * pRz, 1.6075)) / 3.0, 1.0/1.6075) / nsurf;

    int dAinterface;

    droplet = bulk + surf + nsurf;

    printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dA, dApart); 
    printf("dVi = %lf\ndVo = %lf\n", dVi, dVo);
	printf("\nDroplet nodes number is %d.\nBulk nodes number is %d.\nDroplet surface nodes number is %d.\nNanoparticle nodes is %d.\nParticle surface nodes number is %d.\n", droplet, bulk, surf, nanoparticle_nodes, nsurf);
    printf("Nanoparticle interface is %d\n\n", interbulk);
    
    nu = (double*)malloc((surf + nsurf) * 3 * sizeof(double));
	for(int i = 0; i < (surf + nsurf) * 3; i ++){
		nu[i] = 0.;
	}

    if(degenerate == 0 && infinite == 0){
		Qo = (double*)malloc(6 * (surf + nsurf) * sizeof(double));
		for(int i = 0; i < (surf + nsurf) * 6; i ++){
			Qo[i] = 0.;
		}
	}

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
	int nd = 0;

	for (int l = 0; l < total_points; l++) {

		if(!ndrop[l]){
            qindex[l] = nd;
            //bulto: share/sign = 0
            //superficie del canal: share/sign = 2
            //superficie de la nanopartícula: share/sign = 4
            if(drop[l]){       
                signal[nd] = 0;      
            }
            else if(boundary[l]){
                signal[nd] = 2;
                nb++;
            }
            else if(nboundary[l]){
                signal[nd] = 4;
                nb++;
            }
            nd++;
        }
	}

    int nondrop = 0;

    for(int i = 0; i < total_points; i++){
        if(qindex[i] == -1){
            nondrop++;
        }
    }

    //printf("-1 qindex nodes is %d and nanoparticle nodes is %d\n", nondrop, nanoparticle_nodes);

    if (nd != droplet){
		printf("Problem in initialization of qtensor. nd is %d not equal to droplet %d.\n", nd, droplet);
		return false;
        exit(1);
	}
	if (nb != surf + nsurf){
		printf("Problem in initialization of qtensor. nb is %d not equal to surf %d.\n", nb, surf);
        printf("Channel surface nodes: %d; Nanoparticle surface nodes %d; Nanoboundary counter: %d\n", surf, nsurf, nb);
		return false;
		exit(1);
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
    nondrop = 0;
    // time_t t;
    // srand((unsigned) time(&t));
	srand(rand_seed);
    //Defininiendo los vecinos
    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                nd = qindex[i + j * Nx + k * Nx * Ny];

                if(nd == -1){
                    nondrop++;
                    continue;
                }

                if(drop[i + j * Nx + k * Nx * Ny]){
                    neighbor[nd * 6 + 0] = qindex[peri(i - 1, 0) + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 1] = qindex[peri(i + 1, 0) + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 2] = qindex[i + peri(j - 1, 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 3] = qindex[i + peri(j + 1, 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 4] = qindex[i + j * Nx + (k - 1) * Nx * Ny];
					neighbor[nd * 6 + 5] = qindex[i + j * Nx + (k + 1) * Nx * Ny];
                }
                else if(boundary[i + j * Nx + k * Nx * Ny] || nboundary[i + j * Nx + k * Nx * Ny]){
                    if(boundary[i + j * Nx + k * Nx * Ny]){
                        if(k == 0){
                            nu[nb * 3 + 0] = dir1[0];
                            nu[nb * 3 + 1] = dir1[1];
                            nu[nb * 3 + 2] = dir1[2];
                        }
                        else if(k == Nz - 1){
                            nu[nb * 3 + 0] = dir1[0];
                            nu[nb * 3 + 1] = dir1[1];
                            nu[nb * 3 + 2] = -dir1[2];
                        }
                        else{
                            printf("Error in channel surface.\n");
                            return false;
                        }
                        
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
                    else if(nboundary[i + j * Nx + k * Nx * Ny]){
                        x = (double)(i - rx) * dx;
                        y = (double)(j - ry) * dy;
                        z = (double)(k - rz) * dz;

                        x_rot = x * cos(alpha) * cos(beta) + y * (cos(alpha) * sin(beta) * sin(gama) - sin(alpha) * cos(gama))\
                            + z * (cos(alpha) * sin(beta) * cos(gama) + sin(alpha) * sin(gama));
                        y_rot = x * sin(alpha) * cos(beta) + y * (sin(alpha) * sin(beta) *sin(gama) + cos(alpha) * cos(gama))\
                            + z * (sin(alpha) * sin(beta) * cos(gama) - cos(alpha) * sin(gama));
                        z_rot = x * -sin(beta) + y * cos(beta) * sin(gama) + z * cos(beta) * cos(gama);

                        x = x_rot;
                        y = y_rot;
                        z = z_rot;

                        distance = (x * x) / ((pRx + 0.5) * (pRx + 0.5))\
                            + (y * y) / ((pRy + 0.5) * (pRy + 0.5))\
                            + (z * z) / ((pRz + 0.5) * (pRz + 0.5));

                        if (distance == 0){
                            printf("Error in neighbors on particle boundary.\n");
                            return false;
                        }
                        else {

                            if(anchoring == 0){
                                nu[nb * 3 + 0] = (rand() % pRx + 1);
						        nu[nb * 3 + 1] = (rand() % pRy + 1);
						        nu[nb * 3 + 2] = (rand() % pRz + 1);
                                norm_v(&nu[nb * 3]);

                            }
                            else if(anchoring == 1){
                                nu[nb * 3 + 0] = 2. * x / (pRx * pRx);
						        nu[nb * 3 + 1] = 2. * y / (pRy * pRy);
						        nu[nb * 3 + 2] = 2. * z / (pRz * pRz);
						        norm_v(&nu[nb * 3]);
                            }
                            else if(anchoring == 2){
                                printf("Not available yet\n");
                                exit(1);
                            }
                            //La superficie no evoluciona
                            signal[nd] = 8;
                            Qold[nd * 6 + 0] = dir2ten(&nu[nb * 3], 0, S);
                            Qold[nd * 6 + 1] = dir2ten(&nu[nb * 3], 1, S);
                            Qold[nd * 6 + 2] = dir2ten(&nu[nb * 3], 2, S);
                            Qold[nd * 6 + 3] = dir2ten(&nu[nb * 3], 3, S);
                            Qold[nd * 6 + 4] = dir2ten(&nu[nb * 3], 4, S);
                            Qold[nd * 6 + 5] = dir2ten(&nu[nb * 3], 5, S);
                            
                        }
                    }

                    if(nu[nb * 3 + 0] >= 0){
                        neighbor[nd * 6 + 0] = qindex[peri(i + 1, 0) + j * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 1] = qindex[peri(i + 2, 0) + j * Nx + k * Nx * Ny];
                    }
                    else if(nu[nb * 3 + 0] < 0){
                        neighbor[nd * 6 + 0] = qindex[peri(i - 1, 0) + j * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 1] = qindex[peri(i - 2, 0) + j * Nx + k * Nx * Ny];
                    }
                    if(nu[nb * 3 + 1] >= 0){
                        neighbor[nd * 6 + 2] = qindex[i + peri(j + 1, 1) * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 3] = qindex[i + peri(j + 2, 1) * Nx + k * Nx * Ny];
                    }
                    else if(nu[nb * 3 + 1] < 0){
                        neighbor[nd * 6 + 2] = qindex[i + peri(j - 1, 1) * Nx + k * Nx * Ny];
                        neighbor[nd * 6 + 3] = qindex[i + peri(j - 2, 1) * Nx + k * Nx * Ny];
                    }
                    if(nu[nb *3 + 2] >= 0){
                        neighbor[nd * 6 + 4] = qindex[i + j * Nx + (k + 1) * Nx * Ny];
                        neighbor[nd * 6 + 5] = qindex[i + j * Nx + (k + 2) * Nx * Ny];
                    }
                    else if(nu[nb * 3 + 2] < 0){
                        neighbor[nd * 6 + 4] = qindex[i + j * Nx + (k - 1) * Nx * Ny];
                        neighbor[nd * 6 + 5] = qindex[i + j * Nx + (k - 2) * Nx * Ny];
                    }
                    nb++;
                }
            }
        }
    }

    if(nondrop != nanoparticle_nodes){
        printf("Problems in nondrop %d and nanop nodes %d", nondrop, nanoparticle_nodes);
    }

    if (nb != surf + nsurf) {
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
    int t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;
	for (int i = 0; i < total_points; i++) {
		if (bulktype[i] == 1) {
			h_bulktype[nd] = 1;
			nd++;
            t1++;
		}
		else if (bulktype[i] == 2) {
			h_bulktype[nd] = 2;
			nd++;
            t2++;
		}
        else if(bulktype[i] == 3){
            h_bulktype[nd] = 3;
			nd++;
            t3++;
        }
        else if(bulktype[i] == 4 || bulktype[i] == 6){
            h_bulktype[nd] = bulktype[i];
            nd++;
            if(bulktype[i] == 4){
                t4++;
            }
            else{
                t6++;
            }
        }  
	}

	//Size must match.
	if (nd != droplet) {
		printf("Error in transfer data to droplet bulktype!\n");
        printf("Count is %d and droplet is %d + np_count %d!\n", nd, droplet, nanoparticle_nodes);
        exit(1);
        return false;
	}

	int boundaryc = 0;
	int dropc = 0;
	int nboundaryc = 0;

	// for(int i = 0; i < droplet; i++){
	// 	if(h_bulktype[i] == 1){
	// 		if(nboundaryc < 20){
	// 			printf("Nboundary %d : %lf", nboundaryc, Qold[nboundaryc]);
				
	// 		}
	// 		nboundaryc++;
	// 	} 
	// }
	free(qindex);
	free(bulktype);
    printf("Nanochannel initialized successfully!\n");

	return true;

}

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
	unsigned int l = 0, innerbulk = 0, outerbulk = 0;
	//Aquí guardo los índices de los nodos.
	int* qindex = (int*)malloc(total_points * sizeof(int));
	
	bulk = 0;
	droplet = 0;

	dx = Lx / (double)(Nx - 1);
	dy = Ly / (double)(Ny - 1);
	dz = Lz / (double)(Nz - 1);

	rx = lrint(Nx / 2);
	ry = lrint(Ny / 2);
	rz = lrint(Nz / 2);

	Rx = Lx / 2. - 2.;
	Ry = Ly / 2. - 2.;
	Rz = Lz / 2. - 2.;

	surf = 0;
	nsurf = 0;
	
	idx = 1. / dx;
	idy = 1. / dy;
	idz = 1. / dz;

	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

	//boudaries variables
	int xm = 0, xp = 0, ym = 0, yp = 0, zm = 0, zp = 0;

    bulktype = (unsigned char*)malloc(total_points * sizeof(unsigned char));
	drop = (bool*)malloc(total_points * sizeof(bool));
	boundary = (bool*)malloc(total_points * sizeof(bool));
	nboundary = (bool*)malloc(total_points * sizeof(bool));

	for (int i = 0; i < total_points; i++) {
		drop[i] = false;
		boundary[i] = false;
		nboundary[i] = false;
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
						bulktype[l] = 1; //bulk interno.
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
	if(DoubleU){
		if ((innerbulk + outerbulk) != bulk) {
			printf("Problems with bulk nodes!\n");
			return false;
			exit(1);
		}
		else {
			printf("Innerbulk and outerbulk nodes match with total bulk count!\n");
		}
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
	if(DoubleU){
		dVi = ((4.0 / 3.0) * (double)M_PI * (iRx * iRy * iRz)) / (double)innerbulk;
		dVo = ((4.0 / 3.0) * (double)M_PI * ((Rx) * (Ry) * (Rz)) - (4.0 / 3.0) * (double)M_PI * (iRx * iRy * iRz)) / (double)outerbulk;
	}
	dA = (4.0 * (double)M_PI * pow((pow(Rx * Ry, 1.6075) + pow(Rx * Rz, 1.6075) + pow(Ry * Rz, 1.6075)) / 3.0, 1.0000 / 1.6075)) / (double)surf;
	if(DoubleU){
		printf("Internal nodes count = %d\nExternal nodes count = %d\n", innerbulk, outerbulk + surf);
		printf("External bulk count contains surface nodes!!\n");
		printf("External count after substracting surface nodes is %d\n", outerbulk);
	}
	printf("dV = %lf\n", dV);
	if(DoubleU){
		printf("dVi = %lf\n", dVi);
		printf("dVo = %lf\n", dVo);
	}
	
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
		exit(1);
	}
	
	free(qindex);
	free(bulktype);
    printf("Ellipsoid initialized successfully!\n");
	return true;
}
