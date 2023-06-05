#include "definitions.cuh"

int l = 0;
int nd = 0;
double x = 0.0, y = 0.0, z = 0.0;
double A = 0.2;
double cst = 0. , phi = 0.;
double omega, disxy, dis;
double theta = 45 / 180.0 * M_PI;
double xj, yj, zj;
double xi, yi, zi;
double isq2 = 1.0 / sqrt(2.);
double sq2 = sqrt(2.);
//double norm = 0.0;
double dirvec1[3] = {0};
double mod;
double dir[3] = { 0. };
double costhe, sinthe, cosphi, sinphi;
double dir_temp[3] = { 0. };

bool conf() {
	if(seed == -1){
		
		double Qini[6] = { 0. };
		for(int n = 0; n < 6; n ++){
        		Qini[n] = dir2ten(init_dir, n, 0.5);
        	}
		double a[6] = { 0. };
		
		FILE* qtensor;
		qtensor = fopen("Qtensor.bin", "rb");
		FILE* grid;
		grid = fopen("grid.bin", "rb");
		int signal;

		if(qtensor == (FILE*)NULL){
			printf("File Qtensor.bin not found.\n");
			return false;
		}
			
		if(grid == (FILE*)NULL){
			printf("File grid.bin not found.\n");
			return false;
		}
		nd = 0;
		for(int l = 0; l < Nx * Ny * Nz; l++){
			fread(&signal, sizeof(int), 1, grid);
			
			if(signal == 0 || signal == 1){
				fread(a, sizeof(double), 6, qtensor);
					a[5] = - a[0] - a[3];
			}				
			else{
				for (int n = 0; n < 6; n++) {
					a[n] = Qini[n]; 
				}
			}
			
			if(drop[l] || boundary[l] || nboundary[l]){
				for (int n = 0; n < 6; n++) {
					Qold[nd * 6 + n] = a[n];
				}
				nd ++;
			}
			else{
				for(int n = 0; n < 6; n ++){
					Qold[nd * 6 + n] = Qini[n];
				}
			}
		}
		
		fclose(qtensor);
		fclose(grid);
	}
	else if (seed == 4 || seed == 5) {
					
		cst = 2. * qch * redshift;
			
		l = 0;
		nd = 0;

		if(seed == 4){
			for(int k = 0; k < Nz; k++){
				for (int j = 0; j < Ny; j++){
					for (int i = 0; i < Nx; i++){
						
						if(drop[l] || boundary[l] || nboundary[l]){
							if(geo == -2){
								x = (i - rx) * cst * isq2;
								y = (j - 2) * cst * isq2;
								z = (k) * cst * isq2;
							}
							else if(geo == -3){
								x = (i - rx) * cst * isq2;
								y = (j - 2) * cst * isq2;
								z = (k - rz) * cst * isq2;
							}
							else{
								x = (double)(i - rx) * cst * isq2;
								y = (double)(j - ry) * cst * isq2;
								z = (double)(k - rz) * cst * isq2;
							}
										
							Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
							Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
							Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
							Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
							Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
							Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));

							nd ++;
						}
						l ++;

					}
				}
			}
		}
		else{

			//time_t t;
    		//srand((unsigned) time(&t));
			srand(rand_seed);
			double dir_temp[3] = { 0. };

			for(int k = 0; k < Nz; k++){
				for (int j = 0; j < Ny; j++){
					for (int i = 0; i < Nx; i++){
						
						if(drop[l] || boundary[l] || nboundary[l]){
							
							if(geo == -2){

								x = i - rx;
								y = j - 2;
								z = k;

							}
							else if(geo == -3){

								x = i - rx;
								y = j - 2;
								z = k - rz;

							}
							else{

								x = (double)(i - rx) * dx;
								y = (double)(j - ry) * dy;
								z = (double)(k - rz) * dz;

							}

							if(interface != 0 && geo == 10){
								if(bulktype[l] == 3){

									dir_temp[0] = (rand() % (pRx + interface) + 1);
						        	dir_temp[1] = (rand() % (pRy + interface) + 1);
						        	dir_temp[2] = (rand() % (pRz + interface) + 1);
                               		norm_v(dir_temp);

									Qold[nd * 6 + 0] = dir2ten(dir_temp, 0, S2);
									Qold[nd * 6 + 1] = dir2ten(dir_temp, 1, S2);
									Qold[nd * 6 + 2] = dir2ten(dir_temp, 2, S2);
									Qold[nd * 6 + 3] = dir2ten(dir_temp, 3, S2);
									Qold[nd * 6 + 4] = dir2ten(dir_temp, 4, S2);
									Qold[nd * 6 + 5] = dir2ten(dir_temp, 5, S2);
								}
								else{
									Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
									Qold[nd * 6 + 1] = A * sin(cst * z);
									Qold[nd * 6 + 2] = A * sin(cst * y);
									Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
									Qold[nd * 6 + 4] = A * sin(cst * x);
									Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
								}
							}
							else{
								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
								
							}
							nd ++;
						}
						l ++;
					}
				}
			}
		}
	}

    //seed=6; [110] BPI; 7: [110] BPII; 8: [111] BPI; 9: [111] BPII
	else if(seed == 6 || seed == 7 || seed==8 || seed==9){

        
        cst = 2 * qch * redshift;
        
        l = 0;
        nd = 0;

        for(int k = 0; k < Nz; k++){
            for (int j = 0; j < Ny; j++){
                for (int i = 0; i < Nx; i++){
                    if(drop[l] || boundary[l] || nboundary[l]){
                        if(seed == 6){

                            if(geo == -2){
                                xi = (i - rx) * cst * isq2;
                                yi = (j - 2) * cst * isq2;
                                zi = k * cst * isq2;	
                            }

                            else if(geo == -3){

                                xi = (i - rx) * cst * isq2;
                                yi = (j - 2) * cst * isq2;
                                zi = (k - rz) * cst * isq2;

                            }

                            else{

                                xi = (i - rx) * cst * isq2;
                                yi = (j - ry) * cst * isq2;
                                zi = (k - rz) * cst * isq2;

                            }
                            
                            x = xi;
                            y = cos(theta) * yi + sin(theta) * zi;
                            z = -sin(theta) * yi + cos(theta) * zi;

                            Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
                            Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
                            Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
                            Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
                            Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
                            Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
                        }
                        
                        else if(seed == 7){
                            xi = i - rx;
                            yi = j - ry;
                            zi = k - rz;

                            x = xi;
                            y = cos(theta) * yi + sin(theta) * zi;
                            z = -sin(theta) * yi + cos(theta) * zi;
                            Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
                            Qold[nd * 6 + 1] = A * sin(cst * z);
                            Qold[nd * 6 + 2] = A * sin(cst * y);
                            Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
                            Qold[nd * 6 + 4] = A * sin(cst * x);
                            Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
                        }



                        else if(seed == 8){    //8 and 9, (111) planes oriented

                            xi = (i - rx) * cst * isq2;
                            yi = (j - ry) * cst * isq2;
                            zi = (k - rz) * cst * isq2;
                            
                            theta=atan(1.0/sqrt(2.0));
                            //BPI_(211)
                            //Rotation around vector (-1,1,0)
                            x=xi*0.5*(1.0+cos(theta))-0.5*yi*(1.0-cos(theta))+zi*sin(theta)/sqrt(2.0); 
                            y=-0.5*xi*(1.0-cos(theta))+yi*0.5*(1.0+cos(theta))+zi*sin(theta)/sqrt(2.0);
                            z=-xi*sin(theta)/sqrt(2.0)-yi*sin(theta)/sqrt(2.0)+zi*cos(theta);

                            theta= 1.0*M_PI/12.0;
                            xi=x;
                            yi=y;
                            zi=z;
                            //Rotation around the vector (1,1,1) here it should be implemented rotation around (2,1,1)
                            x = xi*1.0/3.0*(2.0*cos(theta)+1.0) + yi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + zi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) ;
                            y = xi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + yi*1.0/3.0*(2.0*cos(theta)+1.0) + zi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) ;
                            z = xi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + yi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + zi*1.0/3.0*(2.0*cos(theta)+1.0) ;

                                
                            Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
                            Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
                            Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
                            Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
                            Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
                            Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
                        }
                        else if(seed == 9){

                            if(geo == -2){
                                xi = i - rx;
                                yi = j - 2;
                                zi = k;
                            }
                            
                            else if(geo == -3){
                                xi = i - rx;
                                yi = j - 2;
                                zi = k - rz;
                            }
                            else{
                                xi = i - rx;
                                yi = j - ry;
                                zi = k - rz;
                            }
                                                
                            theta=atan(sqrt(2.0));

                            //BPII_(111)
                            //Rotation around vector (-1,1,0)
                            x=xi*0.5*(1.0+cos(theta))-0.5*yi*(1.0-cos(theta))+zi*sin(theta)/sqrt(2.0); 
                            y=-0.5*xi*(1.0-cos(theta))+yi*0.5*(1.0+cos(theta))+zi*sin(theta)/sqrt(2.0);
                            z=-xi*sin(theta)/sqrt(2.0)-yi*sin(theta)/sqrt(2.0)+zi*cos(theta);

                            theta= 1.0*M_PI/12.0;
                            xi=x;
                            yi=y;
                            zi=z;

                            //Rotation around the vector (1,1,1) 
                            x = xi*1.0/3.0*(2.0*cos(theta)+1.0) + yi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + zi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) ;
                            y = xi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + yi*1.0/3.0*(2.0*cos(theta)+1.0) + zi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) ;
                            z = xi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + yi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + zi*1.0/3.0*(2.0*cos(theta)+1.0) ;
                                
                            Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
                            Qold[nd * 6 + 1] = A * sin(cst * z);
                            Qold[nd * 6 + 2] = A * sin(cst * y);
                            Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
                            Qold[nd * 6 + 4] = A * sin(cst * x);
                            Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));

                        }
                        nd ++;
                    }
                    l ++;
                }
            }
        }
    }

    // ********** 11 + x = condición radial + condición BP *********** //////////////
	// ***** For BPI seed = 114 and 116 with [110]; and for BPII seed = 115 [110] and 119 [111] ***** //////////
	// ***** Seeds 124, 126, 125, 129 and so on are combinations (Here 1 is a dummy index) of seed 2(DSS) o 3(RSS) (cholesteric)...
	// ***** ... and BPI's seeds 4 & 6, and BPII's seeds 5 & 9.
	// *** for seed 141, 142, 143 BPII 100, 110, 111 outter shells with BPIII random inner core //
	

	else if(seed == 114 || seed == 116 || seed == 115 || seed == 119 ||	
			seed == 124 || seed == 126 || seed == 125 || seed == 129 ||
			seed == 134 || seed == 136 || seed == 135 || seed == 139 || 
			seed == 141 || seed == 142 || seed == 143 || 
			seed == 874 || seed == 875 || seed == 876 || seed == 879 ||
			seed == 884 || seed == 885 || seed == 886 || seed == 879)
		{

		srand(rand_seed);

		
		cst = 2 * qch * redshift;
		

		l = 0;
		nd = 0;

		for(int k = 0; k < Nz; k++){
			for (int j = 0; j < Ny; j++){
				for (int i = 0; i < Nx; i++){

					if(bulktype[l] == 1){
						if(drop[l] || boundary[l]){ // || nboundary[l]){
							if(seed == 114 || seed == 116 || seed == 115 || seed == 119){
								
								
								x = (i - rx) * dx;
								y = (j - ry) * dy;//j * dy;
								z = (k - rz) * dz; //10; //k * dz;    //Modifiqué z = (k - rz) * dz por k * dz

								dir[0] = -x;
								dir[1] = -y;
								dir[2] = -z; //-z;
						
								mod = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]); 

								if (mod == 0){
									dir[0] = 0;
									dir[1] = 0;
									dir[2] = -1; // -10;
									//printf("i: %d, j: %d k: %d \n", i, j, k);
								}

								else{
									dir[0] = dir[0] / mod;
									dir[1] = dir[1] / mod;
									dir[2] = dir[2] / mod;						
								}												

								for (int n = 0; n < 6; n++) {
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}
								nd ++;
							
							}
							else if(seed == 141 || seed == 142 || seed == 143){
								
									for(int n = 0; n < 3; n ++){
										dir[n] = (double)rand() / (double)RAND_MAX * 2 - 1;
									}

									if(!norm_v(dir)){
										printf("Problems with random initialization\n");
										return false;
									}        	

									for(int n = 0; n < 6; n ++){
										Qold[nd * 6 + n] = dir2ten(dir, n, 0.5);
									}            
								 nd ++;
							}
							else if(seed == 124 || seed == 126 || seed == 125 || seed == 129 ||
									seed == 134 || seed == 136 || seed == 135 || seed == 139){
					
								x = (i-rx)*dx;
								y = (j-ry)*dy;
								z = (k-rz)*dz;
								dis = sqrt(x*x+y*y+z*z);

								if(seed == 124 || seed == 126 || seed == 125 || seed == 129){
									omega = dis * qch;
								}   
								
								else{
									omega = atan2(y, x) + dis * qch;
									disxy = sqrt(x * x + y * y);
								}
									
								if(disxy == 0){
									dir[2] = 1;
									dir[0] = dir[1] = 0;
								}

								else{
									costhe = z / dis;
									sinthe = disxy / dis;
									cosphi = x / disxy;
									sinphi = y / disxy;
									dir[0] = cos(omega) * costhe * cosphi - sin(omega) * sinphi;
									dir[1] = cos(omega) * costhe * sinphi + sin(omega) * cosphi;
									dir[2] = - cos(omega) * sinthe;
								}
				
								if(!norm_v(dir)){
									return false;
								}
									
								for (int n = 0; n < 6; n++) {
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}
						
								nd ++;
							}
							else if(seed == 874 || seed == 875 || seed == 876 || seed == 879){
								
								dir[0] = cos(qch * (k - rz));
								dir[1] = sin(qch * (k - rz));

								for(int n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}
								nd ++;								
							}
							else if(seed == 884 || seed == 885 || seed == 886 || seed == 889){
							
								dir[1] = cos(qch * (i - rx));
								dir[2] = sin(qch * (i - rx));

								for(int n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}
								nd ++;
							}
						}
					}

					else if(bulktype[l] == 2){
					
						if(drop[l] || boundary[l]){ // || nboundary[l]){
							if(seed == 114 || seed == 124 || seed == 134 || seed == 874 || seed == 884){

								x = (i - rx) * cst * isq2;
								y = (j - ry) * cst * isq2;
								z = (k - rz) * cst * isq2;
						
								Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
							}

							else if(seed == 115 || seed == 125 || seed == 135 || seed == 141 || seed == 875 || seed == 885){
										
								x = i - rx;
								y = j - ry;
								z = k - rz;

								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));
							}

							else if(seed == 116  || seed == 126 || seed == 136 || seed == 876 || seed == 886){

								xi = (i - rx) * cst * isq2;
								yi = (j - ry) * cst * isq2;
								zi = (k - rz) * cst * isq2;
								
								x = xi;
								y = cos(theta) * yi + sin(theta) * zi;
								z = -sin(theta) * yi + cos(theta) * zi;

								Qold[nd * 6 + 0] = A * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								Qold[nd * 6 + 3] = A * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								Qold[nd * 6 + 5] = A * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								Qold[nd * 6 + 1] = A * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								Qold[nd * 6 + 2] = A * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								Qold[nd * 6 + 4] = A * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
							}
													
							else if (seed == 142){
							
								xi = i - rx;
                        		yi = j - ry;
                            	zi = k - rz;

								x = xi;
								y = cos(theta) * yi + sin(theta) * zi;
								z = -sin(theta) * yi + cos(theta) * zi;

								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));}

							else if(seed == 119  || seed == 129 || seed == 139 || seed == 143 || seed == 879 || seed == 889){

								xi = i - rx;
								yi = j - ry;
								zi = k - rz;
													
								theta=atan(sqrt(2.0));

								//BPII_(111)
								//Rotation around vector (-1,1,0)
								x=xi*0.5*(1.0+cos(theta))-0.5*yi*(1.0-cos(theta))+zi*sin(theta)/sqrt(2.0); 
								y=-0.5*xi*(1.0-cos(theta))+yi*0.5*(1.0+cos(theta))+zi*sin(theta)/sqrt(2.0);
								z=-xi*sin(theta)/sqrt(2.0)-yi*sin(theta)/sqrt(2.0)+zi*cos(theta);

								theta= 1.0*M_PI/12.0;
								xi=x;
								yi=y;
								zi=z;

								//Rotation around the vector (1,1,1) 
								x = xi*1.0/3.0*(2.0*cos(theta)+1.0) + yi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + zi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) ;
								y = xi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + yi*1.0/3.0*(2.0*cos(theta)+1.0) + zi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) ;
								z = xi*((1.0-cos(theta))/3.0-sin(theta)/sqrt(3.0)) + yi*((1.0-cos(theta))/3.0+sin(theta)/sqrt(3.0)) + zi*1.0/3.0*(2.0*cos(theta)+1.0) ;
									
								Qold[nd * 6 + 0] = A * (cos(cst * z) - cos(cst * y));
								Qold[nd * 6 + 1] = A * sin(cst * z);
								Qold[nd * 6 + 2] = A * sin(cst * y);
								Qold[nd * 6 + 3] = A * (cos(cst * x) - cos(cst * z));
								Qold[nd * 6 + 4] = A * sin(cst * x);
								Qold[nd * 6 + 5] = A * (cos(cst * y) - cos(cst * x));

							}

							nd ++;
						}
					}

					l ++;

				}
			}
		}
	}
    
    return true;

}
