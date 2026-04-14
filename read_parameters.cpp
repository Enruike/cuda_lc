#include "read_parameters.hpp"
#define READ_PARAM_LINE(expected_count, ...)            \
	do {                                                \
		char line[512];                                 \
		if (fgets(line, sizeof(line), param) == NULL) { \
			fclose(param);                              \
			printf("Error reading param.in.\n");        \
			return false;                               \
		}                                               \
		if (sscanf(line, __VA_ARGS__) != expected_count) { \
			fclose(param);                              \
			printf("Error reading param.in.\n");        \
			return false;                               \
		}                                               \
	} while (0)

bool read_param() {
	FILE* param = fopen("param.in", "r");

	if (param == (FILE*)NULL) {
		return false;
	}

	else {

		READ_PARAM_LINE(1, "Nx %d", &Nx);
		READ_PARAM_LINE(1, "Ny %d", &Ny);
		READ_PARAM_LINE(1, "Nz %d", &Nz);
		READ_PARAM_LINE(1, "Lx %lf", &Lx);
		READ_PARAM_LINE(1, "Ly %lf", &Ly);
		READ_PARAM_LINE(1, "Lz %lf", &Lz);
		READ_PARAM_LINE(1, "W %lf", &W);
		READ_PARAM_LINE(1, "U %lf", &U);
		READ_PARAM_LINE(1, "L1 %lf", &L1);
		READ_PARAM_LINE(1, "L2 %lf", &L2);
		READ_PARAM_LINE(1, "L3 %lf", &L3);
		READ_PARAM_LINE(1, "L4 %lf", &L4);
		READ_PARAM_LINE(1, "chiral %d", &chiral);
		READ_PARAM_LINE(1, "qch %lf", &qch);
		READ_PARAM_LINE(1, "redshift %lf", &redshift);
		READ_PARAM_LINE(1, "geo %d", &geo);
		READ_PARAM_LINE(1, "degenerate %d", &degenerate);
		READ_PARAM_LINE(1, "tiltAngle %lf", &tiltAngle);
		READ_PARAM_LINE(1, "infinite %d", &infinite);
		READ_PARAM_LINE(1, "Np %d", &Np);
		READ_PARAM_LINE(1, "Rp %lf", &Rp);
		READ_PARAM_LINE(1, "Wp %lf", &Wp);
		READ_PARAM_LINE(1, "seed %d", &seed);
		READ_PARAM_LINE(1, "rand_seed %d", &rand_seed);
		READ_PARAM_LINE(2, "tmin tmax %lf %lf", &tmin, &tmax);
		READ_PARAM_LINE(1, "increment %lf", &increment);
		READ_PARAM_LINE(1, "accuracy %lf", &accuracy);
		READ_PARAM_LINE(3, "init_dir %lf %lf %lf", &init_dir[0], &init_dir[1], &init_dir[2]);
		READ_PARAM_LINE(3, "dir1 %lf %lf %lf", &dir1[0], &dir1[1], &dir1[2]);
		READ_PARAM_LINE(3, "dir2 %lf %lf %lf", &dir2[0], &dir2[1], &dir2[2]);
		READ_PARAM_LINE(1, "UpperSurface %d", &uppersurf);
		READ_PARAM_LINE(1, "LowerSurface %d", &lowersurf);
		READ_PARAM_LINE(1, "LowerSurfaceDegen %d", &surfdegen);
		READ_PARAM_LINE(1, "DoubleU Mode %d", &DoubleU);
		READ_PARAM_LINE(1, "U2 %lf", &U2);
		READ_PARAM_LINE(1, "iRx %lf", &iRx);
		READ_PARAM_LINE(1, "iRy %lf", &iRy);
		READ_PARAM_LINE(1, "iRz %lf", &iRz);
		READ_PARAM_LINE(1, "Save Every %d", &save_every);
		READ_PARAM_LINE(1, "Check Every %d", &check_every);
		READ_PARAM_LINE(1, "Stop At %d", &stopat);
		READ_PARAM_LINE(1, "Check Trace At %d", &trace_checker);

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
