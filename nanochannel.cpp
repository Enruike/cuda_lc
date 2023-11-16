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

    printf("\n~ Nanoparticle data ~\n");
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

    double x = 0.0, y = 0.0, z = 0.0, dis = 0.0;
	unsigned int l = 0, bulk = 0;
    unsigned int surf = 0, nsurf;
	
	
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
    bulk = Nx * Ny * (Nz - 2);
	
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
	qindex = (int*)malloc(total_points * sizeof(int));

    //boudaries variables
	int xm = 0, xp = 0, ym = 0, yp = 0, zm = 0, zp = 0;

	for (int i = 0; i < total_points; i++) {
		drop[i] = false;
		boundary[i] = false;
        ndrop[i] = false;
		bulktype[i] = 0;
		qindex[i] = -1;
	}

    //Definir superficie del nanocanal
    for(int k = 0; k < Nz; k += (Nz - 1)){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){
                l = i + j * Nx + k * Nx * Ny;
                boundary[l] = true;
                surf++;
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

    double x_rot, y_rot, z_rot;
    double distance;
    l = 0;
    int nanoparticle_nodes = 0;

    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                x = i - rx;
                y = j - ry;
                z = k - rz;

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
                    xm = (i - 1) + j * Nx + k * Nx * Ny;
                    xp = (i + 1) + j * Nx + k * Nx * Ny;
                    ym = i + (j - 1) * Nx + k * Nx * Ny;
                    yp = i + (j + 1) * Nx + k * Nx * Ny;
                    zm = i + j * Nx + (k - 1) * Nx * Ny;
                    zp = i + j * Nx + (k + 1) * Nx * Ny;
                    
                    if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zm] || ndrop[zp]){
						nboundary[l] = true;
						drop[l] = false;
                        bulktype[l] = 6;
						bulk--;
						nsurf++;	
						surf++;
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
                            xm =(i - (node + 1)) + j * Nx + k * Nx * Ny;
                            xp = (i + (node + 1)) + j * Nx + k * Nx * Ny;
                            ym = i + (j - (node + 1)) * Nx + k * Nx * Ny;
                            yp = i + (j + (node + 1)) * Nx + k * Nx * Ny;
                            zm = i + j * Nx + (k - (node + 1)) * Nx * Ny;
                            zp = i + j * Nx + (k + (node + 1)) * Nx * Ny;
                            
                            
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
    dAdrop = (2 * Lx * Ly) / (surf - nsurf);
    dApart = 4. * M_PI * pow((pow(pRx * pRy, 1.6075) + pow(pRx * pRz, 1.6075) + pow(pRy * pRz, 1.6075)) / 3.0, 1.0/1.6075) / nsurf;

    int dAinterface;

    printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 
    printf("dVi = %lf\ndVo = %lf\n", dVi, dVo);
    droplet = bulk + surf;
	printf("\nDroplet nodes number is %d.\nBulk nodes number is %d.\nDroplet surface nodes number is %d. \nParticle surface nodes number is %d.\n", droplet, bulk, surf, nsurf); 
    printf("Nanoparticle interface is %d\n", interbulk);

    nu = (double*)malloc(surf * 3 * sizeof(double));
	for(int i = 0; i < surf * 3; i ++){
		nu[i] = 0;
	}

    if(degenerate == 0 && infinite == 0){
		Qo = (double*)malloc(6 * surf * sizeof(double));
		for(int i = 0; i < surf * 6; i ++){
			Qo[i] = 0;
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
	unsigned int nd = 0;

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

    if (nd != droplet){
		printf("Problem in initialization of qtensor. nd is %d not equal to droplet %d.\n", nd, droplet);
		return false;
	}
	if (nb != surf){
		printf("Problem in initialization of qtensor. nb is %d not equal to surf %d.\n", nb, surf);
        printf("Channel surface nodes: %d; Nanoparticle surface nodes %d; Nanoboundary counter: %d\n", surf, nsurf, nb);
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

    time_t t;
    srand((unsigned) time(&t));
    //Defininiendo los vecinos
    for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){

                nd = qindex[i + j * Nx + k * Nx * Ny];

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
                        x = i - rx;
                        y = j - ry;
                        z = k - rz;

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
                        nsurf--;
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
        else if(bulktype[i] == 3){
            h_bulktype[nd] = 3;
			nd++;
        }
        else if(bulktype[i] == 4 || bulktype[i] == 6){
                h_bulktype[nd] = bulktype[i];
                nd++;
            }  
	}

	//Size must match.
	if (nd != droplet) {
		printf("Error in transfer data to droplet bulktype!\n");
		return false;
	}

	
	printf("Nanochannel initialized successfully!\n");
	free(qindex);
	free(bulktype);

	return true;

}