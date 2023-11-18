#include "read_parameters.hpp"

bool read_param() {
	FILE* param = fopen("param.in", "r");

	if (param == (FILE*)NULL) {
		return false;
	}

	else {

        fscanf(param, "Nx %d\n", &Nx);
        fscanf(param, "Ny %d\n", &Ny);
        fscanf(param, "Nz %d\n", &Nz);
        fscanf(param, "Lx %lf\n", &Lx);
        fscanf(param, "Ly %lf\n", &Ly);
        fscanf(param, "Lz %lf\n", &Lz);
        fscanf(param, "W %lf\n", &W);
        fscanf(param, "U %lf\n", &U);
        fscanf(param, "L1 %lf\n", &L1);
        fscanf(param, "L2 %lf\n", &L2);
        fscanf(param, "L3 %lf\n", &L3);
        fscanf(param, "L4 %lf\n", &L4);
        fscanf(param, "chiral %d\n", &chiral);
        fscanf(param, "qch %lf\n", &qch);
		fscanf(param, "redshift %lf #Default values BPI = 0.71 & BPII = 0.86.\n", &redshift);
        fscanf(param, "geo %d\n", &geo);
        fscanf(param, "degenerate %d #1:Degenerated & 2:Conic\n", &degenerate);
		fscanf(param, "tiltAngle %lf\n", &tiltAngle);
        fscanf(param, "infinite %d\n", &infinite);
        fscanf(param, "infinite %d\n", &infinite);
        fscanf(param, "Np %d\n", &Np);
        fscanf(param, "Rp %lf\n", &Rp);
        fscanf(param, "Wp %lf\n", &Wp);
        fscanf(param, "seed %d\n", &seed);
        fscanf(param, "rand_seed %d\n", &rand_seed);
        fscanf(param, "tmin tmax %lf %lf\n", &tmin, &tmax);
        fscanf(param, "increment %lf\n", &increment);
        fscanf(param, "accuracy %lf\n", &accuracy);
        fscanf(param, "init_dir %lf %lf %lf\n", &init_dir[0], &init_dir[1], &init_dir[2]);
        fscanf(param, "dir1 %lf %lf %lf\n", &dir1[0], &dir1[1], &dir1[2]);
        fscanf(param, "dir2 %lf %lf %lf\n", &dir2[0], &dir2[1], &dir2[2]);
        fscanf(param, "UpperSurface %d\n", &uppersurf);
        fscanf(param, "LowerSurface %d\n", &lowersurf);
        fscanf(param, "LowerSurfaceDegen %d\n", &surfdegen);
        fscanf(param, "DoubleU Mode %d\n", &DoubleU);
        fscanf(param, "U2 %lf\n", &U2);
        fscanf(param, "iRx %lf\n", &iRx);
        fscanf(param, "iRy %lf\n", &iRy);
        fscanf(param, "iRz %lf\n", &iRz);
        fscanf(param, "Save Every %d\n", &save_every);
        fscanf(param, "Check Every %d\n", &check_every);
		fscanf(param, "Stop At %d #For non-stop condition use 0.\n", &stopat);
		fscanf(param, "Check Trace At %d\n", &trace_checker);

		char* status;

		printf("Nx %d\n", Nx);
		printf("Ny %d\n", Ny);
		printf("Nz %d\n", Nz);
		printf("Lx %lf\n", Lx);
		printf("Ly %lf\n", Ly);
		printf("Lz %lf\n", Lz);
		printf("W %lf\n", W);
		printf("U %lf\n", U);
		printf("L1 %lf\n", L1);
		printf("L2 %lf\n", L2);
		printf("L3 %lf\n", L3);
		printf("L4 %lf\n", L4);
		printf("chiral %d\n", chiral);
		printf("qch %lf\n", qch);
		printf("redshift %lf\n", redshift);
		printf("geo %d\n", geo);
		printf("degenerate %d\n", degenerate);
		printf("tilt angle %.2lf°\n", tiltAngle);
		printf("infinite %d\n", infinite);
		printf("Np %d\n", Np);
		printf("Rp %lf\n", Rp);
		printf("Wp %lf\n", Wp);
		//	printf("pdegenerate %d\n", pdegenerate);
		//	printf("pinfinite %d\n", pinfinite);
		printf("seed %d\n", seed);
		printf("rand_seed %d\n", rand_seed);
		printf("tmin tmax %lf %lf\n", tmin, tmax);
		printf("increment %lf\n", increment);
		printf("accuracy %lf\n", accuracy);
		printf("init_dir is %lf,%lf,%lf\n", init_dir[0], init_dir[1], init_dir[2]);
		printf("dir1 is %lf,%lf,%lf\n", dir1[0], dir1[1], dir1[2]);
		printf("dir2 is %lf,%lf,%lf\n", dir2[0], dir2[1], dir2[2]);
		printf("Upper Surface is %d\n", uppersurf);
		printf("Lower Surface is %d\n", lowersurf);
		printf("Checkpoint every %d!\n", save_every);
		printf("Energy will be compared every %d!\n", check_every);
		printf("Trace will be checked every %d!\n", trace_checker);

		//My new variables for dynamic savings.
		printf("Checkpoint every %d!\n", save_every);
		printf("Energy will be compared every %d!\n", check_every);
		if(stopat != 0){
			printf("Job will be\033[1;31m STOPPED\033[0m after %d!\n", stopat);
		}
		else{
			printf("Job running until dE is reached!\n");
		}

		if (surfdegen) {
			status = "Activated";
		}
		else {
			status = "Deactivated";
		}

		printf("Lower Surface Degenerate: %s\n", status);

		if (infinite && degenerate) {
			printf("Degenerate planar anchoring cannot be infinite.\n");
			return false;
		}

		if (!norm_v(dir1) || !norm_v(dir2)) {
			printf("Invalid dir1 or dir2 input.\n");
			return false;
		}

		if (tmax < tmin) {
			printf("Error for time input.\n");
			return false;
		}

		if (!norm_v(init_dir)) {
			printf("Problems in initial direction!\n");
			return false;
		}

		if (lowersurf == 0 && dir2[2] != 0) {
			printf("Error in dir2! Planar can only be in x or z direction.\n");
			return false;
		}
		/*Por ahora dejaré el comando así, pero después podremos configurarlo para que pueda
		configuraciones inclinadas.		*/
		else if (lowersurf == 1 && (dir2[0] != 0 || dir2[1] != 0)) {
			printf("Error in dir2! Homeotropic can only be in z direction.\n");
			return false;
		}

		//por ahora mantendremos normales los vectores en la superficie superior.
		if (uppersurf == 1 && (dir1[0] != 0 && dir1[1] != 0)) {
			printf("Error in dir1! Homeotropic can't be in x or z direction.\n");
			return false;
		}
		else if (uppersurf == 0 && dir1[2] == 1) {
			printf("Error in dir1! Planar can't be in y direction.\n");
			return false;
		}

		if (DoubleU) {
			printf("Double U mode is activated.\n");
			printf("U for outer shell = %lf\n", U2);
			printf("Inner Radii iRx = %lf; iRy = %lf; iRz = %lf\n\n", iRx, iRy, iRz);
		}
		else {
			printf("Double U mode is NOT activated.\n\n");
		}
	

		fclose(param);
		return true;
	}
}
