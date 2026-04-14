#ifndef __SOLVER_COMMON_CUH__
#define __SOLVER_COMMON_CUH__

__device__ double fabs(double x);
__device__ double cos(double x);
__device__ double pow(double x, double y);

static __device__ __constant__ double devThird = 1. / 3.;
static __device__ __constant__ double devPI = 3.14159265358979323846;
static __device__ double delta[6] = { 1., 0., 0., 1., 0., 1. };

static __device__ double trQQ_f(double Q[6]) {
	return Q[0] * Q[0] + Q[3] * Q[3] + Q[5] * Q[5]
		+ 2. * (Q[1] * Q[1] + Q[2] * Q[2] + Q[4] * Q[4]);
}

static __device__ double trace_f(double Q[6]) {
	return devThird * (Q[0] + Q[3] + Q[5]);
}

static __device__ double surface_derivative_component(const double* d_Qold, int minus_idx, int plus_idx, double q_center, int comp, double inv_d) {
	if (minus_idx == -1) {
		return 0.;
	}
	if (plus_idx == -1) {
		return (d_Qold[minus_idx * 6 + comp] - q_center) * inv_d;
	}
	return (-d_Qold[plus_idx * 6 + comp] + 4. * d_Qold[minus_idx * 6 + comp] - 3. * q_center) * 0.5 * inv_d;
}

#endif
