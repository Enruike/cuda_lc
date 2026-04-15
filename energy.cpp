#include"energy.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#ifndef DEBUG_ENERGY_COMPARE
#define DEBUG_ENERGY_COMPARE 0
#endif

extern unsigned int bulk, surf, nsurf;
extern double* d_Qold;
extern unsigned char* d_bulktype;
extern double* d_Qo;
extern int* d_neighbor;
extern unsigned int* d_Qtensor_index;
extern unsigned int* d_Nvector_index;
extern unsigned char* d_Nvector_signal;
extern double* d_nu;

namespace {

enum EnergyTermIndex {
	LDG_TOTAL = 0,
	LDG_IN = 1,
	LDG_OUT = 2,
	SURF_DROP = 3,
	SURF_PART = 4,
	EL_TOT_L1 = 5,
	EL_TOT_L2 = 6,
	EL_TOT_L3 = 7,
	EL_TOT_L4 = 8,
	EL_TOT_CHIRAL = 9,
	EL_IN_L1 = 10,
	EL_IN_L2 = 11,
	EL_IN_L3 = 12,
	EL_IN_L4 = 13,
	EL_IN_CHIRAL = 14,
	EL_OUT_L1 = 15,
	EL_OUT_L2 = 16,
	EL_OUT_L3 = 17,
	EL_OUT_L4 = 18,
	EL_OUT_CHIRAL = 19,
	ENERGY_TERMS_COUNT = 20
};

double* d_energy_terms = nullptr;
double* d_bulk_block_terms = nullptr;
double* d_surf_block_terms = nullptr;
unsigned int cached_bulk_blocks = 0;
unsigned int cached_surf_blocks = 0;

__device__ void reduce_block_terms(double* shared_terms, unsigned int tid, unsigned int block_size) {
	for (unsigned int offset = block_size / 2; offset > 0; offset >>= 1) {
		__syncthreads();
		if (tid < offset) {
			for (int term = 0; term < ENERGY_TERMS_COUNT; term++) {
				shared_terms[tid * ENERGY_TERMS_COUNT + term] += shared_terms[(tid + offset) * ENERGY_TERMS_COUNT + term];
			}
		}
	}
	__syncthreads();
}

__device__ void flush_block_terms(double* block_terms, double* shared_terms, unsigned int tid) {
	if (tid == 0) {
		for (int term = 0; term < ENERGY_TERMS_COUNT; term++) {
			block_terms[blockIdx.x * ENERGY_TERMS_COUNT + term] = shared_terms[term];
		}
	}
}

__device__ double trqq_dev(const double Q[6]) {
	return Q[0] * Q[0] + Q[3] * Q[3] + Q[5] * Q[5]
		+ 2. * (Q[1] * Q[1] + Q[2] * Q[2] + Q[4] * Q[4]);
}

__device__ double trqqq_dev(const double Q[6]) {
	return Q[0] * Q[0] * Q[0] + Q[3] * Q[3] * Q[3] + Q[5] * Q[5] * Q[5]
		+ 6. * Q[1] * Q[2] * Q[4] + 3. * Q[0] * (Q[1] * Q[1] + Q[2] * Q[2])
		+ 3. * Q[3] * (Q[1] * Q[1] + Q[4] * Q[4]) + 3. * Q[5] * (Q[4] * Q[4] + Q[2] * Q[2]);
}

__device__ double matr_mult_dev(const double vec[3]) {
	return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]
		+ 2. * vec[0] * vec[1] + 2. * vec[0] * vec[2] + 2. * vec[1] * vec[2];
}

__device__ void en_degen_dev(const double Qin[6], const double loc_nu[3], double Qdiff[6], double S0_val) {
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	const double third = 1.0 / 3.0;

	Qtemp[0][0] = Qin[0] + third * S0_val;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + third * S0_val;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + third * S0_val;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) ptemp[i][j] = 1. - loc_nu[i] * loc_nu[j];
			else ptemp[i][j] = -loc_nu[i] * loc_nu[j];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Qp[i][j] = 0.;
			for (int l = 0; l < 3; l++) {
				for (int m = 0; m < 3; m++) {
					Qp[i][j] += ptemp[i][l] * Qtemp[l][m] * ptemp[m][j];
				}
			}
		}
	}

	Qdiff[0] = Qtemp[0][0] - Qp[0][0];
	Qdiff[1] = Qtemp[0][1] - Qp[0][1];
	Qdiff[2] = Qtemp[0][2] - Qp[0][2];
	Qdiff[3] = Qtemp[1][1] - Qp[1][1];
	Qdiff[4] = Qtemp[1][2] - Qp[1][2];
	Qdiff[5] = Qtemp[2][2] - Qp[2][2];
}

__device__ void en_conic_dev(const double Qin[6], const double loc_nu[3], double Qdiff[6], double S0_val, double tiltAngleVal) {
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	const double third = 1.0 / 3.0;

	Qtemp[0][0] = Qin[0] + third * S0_val;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + third * S0_val;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + third * S0_val;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			ptemp[i][j] = loc_nu[i] * loc_nu[j];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Qp[i][j] = 0.;
			for (int l = 0; l < 3; l++) {
				for (int m = 0; m < 3; m++) {
					Qp[i][j] += ptemp[i][l] * Qtemp[l][m] * ptemp[m][j];
				}
			}
		}
	}

	const double cosTiltAngle = cos(tiltAngleVal / 180.0 * M_PI);
	const double cosTiltAngleSq = cosTiltAngle * cosTiltAngle;

	Qdiff[0] = Qp[0][0] - cosTiltAngleSq * S0_val * ptemp[0][0];
	Qdiff[1] = Qp[0][1] - cosTiltAngleSq * S0_val * ptemp[0][1];
	Qdiff[2] = Qp[0][2] - cosTiltAngleSq * S0_val * ptemp[0][2];
	Qdiff[3] = Qp[1][1] - cosTiltAngleSq * S0_val * ptemp[1][1];
	Qdiff[4] = Qp[1][2] - cosTiltAngleSq * S0_val * ptemp[1][2];
	Qdiff[5] = Qp[2][2] - cosTiltAngleSq * S0_val * ptemp[2][2];
}

__global__ void bulk_energy_kernel(double* block_terms, const double* d_Qold, const unsigned char* d_bulktype,
	const int* d_neighbor, const unsigned int* d_Qtensor_index, unsigned int bulk, int chiral, double U, double U2,
	double idx, double idy, double idz) {
	extern __shared__ double shared_terms[];
	const unsigned int tid = threadIdx.x;
	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;
	double* local_terms = &shared_terms[tid * ENERGY_TERMS_COUNT];
	for (int term = 0; term < ENERGY_TERMS_COUNT; term++) {
		local_terms[term] = 0.;
	}

	if (indx < bulk) {

		const unsigned int Q_indx = d_Qtensor_index[indx];
		const unsigned char bulk_type = d_bulktype[Q_indx];

		double Qin[6];
		for (int n = 0; n < 6; n++) {
			Qin[n] = d_Qold[Q_indx * 6 + n];
		}

		const double trace2 = trqq_dev(Qin);
		const double trace3 = trqqq_dev(Qin);
		const bool inner_region = (bulk_type == 1 || bulk_type == 13);
		const double local_u = inner_region ? U : U2;
		const double ldg_val = 0.5 * (1. - local_u / 3.) * trace2 - local_u / 3. * trace3 + local_u * 0.25 * trace2 * trace2;

		local_terms[LDG_TOTAL] += ldg_val;
		if (inner_region) local_terms[LDG_IN] += ldg_val;
		else local_terms[LDG_OUT] += ldg_val;

		const int xm = d_neighbor[Q_indx * 6 + 0];
		const int xp = d_neighbor[Q_indx * 6 + 1];
		const int ym = d_neighbor[Q_indx * 6 + 2];
		const int yp = d_neighbor[Q_indx * 6 + 3];
		const int zm = d_neighbor[Q_indx * 6 + 4];
		const int zp = d_neighbor[Q_indx * 6 + 5];

		double dQ[3][6];
		for (int n = 0; n < 6; n++) {
			dQ[0][n] = (d_Qold[xp * 6 + n] - d_Qold[xm * 6 + n]) * 0.5 * idx;
			dQ[1][n] = (d_Qold[yp * 6 + n] - d_Qold[ym * 6 + n]) * 0.5 * idy;
			dQ[2][n] = (d_Qold[zp * 6 + n] - d_Qold[zm * 6 + n]) * 0.5 * idz;
		}

		const double l1_val = trqq_dev(dQ[0]) + trqq_dev(dQ[1]) + trqq_dev(dQ[2]);
		local_terms[EL_TOT_L1] += l1_val;
		if (inner_region) local_terms[EL_IN_L1] += l1_val;
		else local_terms[EL_OUT_L1] += l1_val;

		double vec[3];
		vec[0] = dQ[0][0];
		vec[1] = dQ[1][1];
		vec[2] = dQ[2][2];
		double l2_val = matr_mult_dev(vec);
		vec[0] = dQ[0][1];
		vec[1] = dQ[1][3];
		vec[2] = dQ[2][4];
		l2_val += matr_mult_dev(vec);
		vec[0] = dQ[0][2];
		vec[1] = dQ[1][4];
		vec[2] = dQ[2][5];
		l2_val += matr_mult_dev(vec);
		local_terms[EL_TOT_L2] += l2_val;
		if (inner_region) local_terms[EL_IN_L2] += l2_val;
		else local_terms[EL_OUT_L2] += l2_val;

		if (chiral == 1) {
			const double chiral_val = Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] + Qin[1] * dQ[2][0]
				+ Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2] + Qin[2] * dQ[0][1] + Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4]
				- Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4] - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1]
				- Qin[5] * dQ[1][2] - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];
			local_terms[EL_TOT_CHIRAL] += chiral_val;
			if (inner_region) local_terms[EL_IN_CHIRAL] += chiral_val;
			else local_terms[EL_OUT_CHIRAL] += chiral_val;
		}
	}

	reduce_block_terms(shared_terms, tid, blockDim.x);
	flush_block_terms(block_terms, shared_terms, tid);
}

__global__ void surface_energy_kernel(double* block_terms, const double* d_Qold, const unsigned char* d_Nvector_signal,
	const unsigned int* d_Nvector_index, const double* d_Qo, const double* d_nu, unsigned int surface_nodes,
	int degenerate, int infinite, double W, double Wp, double dA, double dApart, double S0_val, double tiltAngleVal) {
	extern __shared__ double shared_terms[];
	const unsigned int tid = threadIdx.x;
	unsigned int indx = threadIdx.x + blockDim.x * blockIdx.x;
	double* local_terms = &shared_terms[tid * ENERGY_TERMS_COUNT];
	for (int term = 0; term < ENERGY_TERMS_COUNT; term++) {
		local_terms[term] = 0.;
	}

	if (indx < surface_nodes) {
		const unsigned char sig = d_Nvector_signal[indx];
		int degen = 0;
		int inf = 1;
		double Wstr = 0.;
		double area = dA;
		bool npboundary = false;

		if (sig == 2 || sig == 3) {
			degen = degenerate;
			inf = infinite;
			Wstr = W;
			area = dA;
		} else if (sig == 12 || sig == 13) {
			degen = 1;
			inf = infinite;
			Wstr = W;
			area = dA;
		} else if (sig == 4 || sig == 5) {
			degen = 0;
			inf = 0;
			Wstr = Wp;
			area = dApart;
			npboundary = true;
		} else if (sig == 6 || sig == 7) {
			degen = 1;
			inf = 0;
			Wstr = Wp;
			area = dApart;
			npboundary = true;
		} else if (sig == 20 || sig == 21) {
			degen = 1;
			inf = 0;
			Wstr = Wp;
			area = dApart;
			npboundary = true;
		} else if (sig == 22 || sig == 23) {
			degen = 2;
			inf = 0;
			Wstr = Wp;
			area = dApart;
			npboundary = true;
		} else if (sig == 8) {
			degen = 0;
			inf = 1;
			Wstr = Wp;
			area = dApart;
			npboundary = true;
		}

		if (Wstr != 0. && inf != 1) {
			const unsigned int Q_indx = d_Nvector_index[indx];
			double Qin[6];
			double Qdiff[6];
			double loc_nu[3];

			for (int n = 0; n < 6; n++) {
				Qin[n] = d_Qold[Q_indx * 6 + n];
			}
			for (int n = 0; n < 3; n++) {
				loc_nu[n] = d_nu[indx * 3 + n];
			}

			if (degen == 1) {
				en_degen_dev(Qin, loc_nu, Qdiff, S0_val);
			} else if (degen == 2) {
				en_conic_dev(Qin, loc_nu, Qdiff, S0_val, tiltAngleVal);
			} else {
				for (int n = 0; n < 6; n++) {
					Qdiff[n] = Qin[n] - d_Qo[indx * 6 + n];
				}
			}

			const double surf_val = Wstr * trqq_dev(Qdiff) * area;
			if (npboundary) local_terms[SURF_PART] += surf_val;
			else local_terms[SURF_DROP] += surf_val;
		}
	}

	reduce_block_terms(shared_terms, tid, blockDim.x);
	flush_block_terms(block_terms, shared_terms, tid);
}

} // namespace

#if DEBUG_ENERGY_COMPARE
static void debug_compare_cpu_gpu_energy() {
	cudaError_t status = cudaMemcpy(Qold, d_Qold, droplet * 6 * sizeof(double), cudaMemcpyDeviceToHost);
	if (status != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed for CPU/GPU energy comparison: %s\n", cudaGetErrorString(status));
		exit(1);
	}

	double cpu_ldg[3] = { 0. };
	double cpu_el[5] = { 0. };
	double cpu_el_in[5] = { 0. };
	double cpu_el_out[5] = { 0. };
	double cpu_surf[2] = { 0. };

	ldg_energy(cpu_ldg);
	elastic_energy(cpu_el, cpu_el_in, cpu_el_out);
	surface_energy(cpu_surf);

	const double cpu_tot = cpu_ldg[0] + cpu_el[0] + cpu_el[1] + cpu_el[2] + cpu_el[3] + cpu_el[4] + cpu_surf[0] + cpu_surf[1];

	printf("[energy cmp] cycle %d\n", cycle);
	printf("[energy cmp] GPU total=%0.12lf CPU total=%0.12lf diff=%+.12e\n", en_tot, cpu_tot, en_tot - cpu_tot);
	printf("[energy cmp] GPU LdG=%0.12lf CPU LdG=%0.12lf diff=%+.12e\n", en_ldg[0], cpu_ldg[0], en_ldg[0] - cpu_ldg[0]);
	printf("[energy cmp] GPU L1=%0.12lf CPU L1=%0.12lf diff=%+.12e\n", en_el[0], cpu_el[0], en_el[0] - cpu_el[0]);
	printf("[energy cmp] GPU L2=%0.12lf CPU L2=%0.12lf diff=%+.12e\n", en_el[1], cpu_el[1], en_el[1] - cpu_el[1]);
	printf("[energy cmp] GPU Surf1=%0.12lf CPU Surf1=%0.12lf diff=%+.12e\n", en_surf[0], cpu_surf[0], en_surf[0] - cpu_surf[0]);
	printf("[energy cmp] GPU SurfP=%0.12lf CPU SurfP=%0.12lf diff=%+.12e\n", en_surf[1], cpu_surf[1], en_surf[1] - cpu_surf[1]);
}
#endif

void free_energy() {
	if (d_energy_terms == nullptr) {
		cudaError_t status = cudaMalloc((void**)&d_energy_terms, sizeof(double) * ENERGY_TERMS_COUNT);
		if (status != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed for energy terms: %s\n", cudaGetErrorString(status));
			exit(1);
		}
	}

	cudaError_t status = cudaMemset(d_energy_terms, 0, sizeof(double) * ENERGY_TERMS_COUNT);
	if (status != cudaSuccess) {
		fprintf(stderr, "cudaMemset failed for energy terms: %s\n", cudaGetErrorString(status));
		exit(1);
	}

	const unsigned int threads = 256;
	const size_t shared_bytes = sizeof(double) * ENERGY_TERMS_COUNT * threads;
	unsigned int bulk_blocks = 0;
	unsigned int surf_blocks = 0;
	if (bulk > 0) {
		bulk_blocks = (bulk + threads - 1) / threads;
		if (d_bulk_block_terms == nullptr || cached_bulk_blocks != bulk_blocks) {
			if (d_bulk_block_terms != nullptr) cudaFree(d_bulk_block_terms);
			status = cudaMalloc((void**)&d_bulk_block_terms, sizeof(double) * ENERGY_TERMS_COUNT * bulk_blocks);
			if (status != cudaSuccess) {
				fprintf(stderr, "cudaMalloc failed for bulk energy partials: %s\n", cudaGetErrorString(status));
				exit(1);
			}
			cached_bulk_blocks = bulk_blocks;
		}
		bulk_energy_kernel<<<bulk_blocks, threads, shared_bytes>>>(d_bulk_block_terms, d_Qold, d_bulktype, d_neighbor, d_Qtensor_index, bulk,
			chiral, U, U2, idx, idy, idz);
	}
	if ((surf + nsurf) > 0) {
		surf_blocks = (surf + nsurf + threads - 1) / threads;
		if (d_surf_block_terms == nullptr || cached_surf_blocks != surf_blocks) {
			if (d_surf_block_terms != nullptr) cudaFree(d_surf_block_terms);
			status = cudaMalloc((void**)&d_surf_block_terms, sizeof(double) * ENERGY_TERMS_COUNT * surf_blocks);
			if (status != cudaSuccess) {
				fprintf(stderr, "cudaMalloc failed for surface energy partials: %s\n", cudaGetErrorString(status));
				exit(1);
			}
			cached_surf_blocks = surf_blocks;
		}
		surface_energy_kernel<<<surf_blocks, threads, shared_bytes>>>(d_surf_block_terms, d_Qold, d_Nvector_signal, d_Nvector_index, d_Qo,
			d_nu, surf + nsurf, degenerate, infinite, W, Wp, dA, dApart, S0, tiltAngle);
	}

	status = cudaGetLastError();
	if (status != cudaSuccess) {
		fprintf(stderr, "Energy kernel launch failed: %s\n", cudaGetErrorString(status));
		exit(1);
	}

	status = cudaDeviceSynchronize();
	if (status != cudaSuccess) {
		fprintf(stderr, "Energy kernel execution failed: %s\n", cudaGetErrorString(status));
		exit(1);
	}

	double host_energy_terms[ENERGY_TERMS_COUNT] = { 0. };
	if (bulk_blocks > 0) {
		double* host_bulk_terms = (double*)malloc(sizeof(double) * ENERGY_TERMS_COUNT * bulk_blocks);
		if (host_bulk_terms == nullptr) {
			fprintf(stderr, "malloc failed for bulk energy partials.\n");
			exit(1);
		}
		status = cudaMemcpy(host_bulk_terms, d_bulk_block_terms, sizeof(double) * ENERGY_TERMS_COUNT * bulk_blocks, cudaMemcpyDeviceToHost);
		if (status != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for bulk energy partials: %s\n", cudaGetErrorString(status));
			exit(1);
		}
		for (unsigned int block = 0; block < bulk_blocks; block++) {
			for (int term = 0; term < ENERGY_TERMS_COUNT; term++) {
				host_energy_terms[term] += host_bulk_terms[block * ENERGY_TERMS_COUNT + term];
			}
		}
		free(host_bulk_terms);
	}
	if (surf_blocks > 0) {
		double* host_surf_terms = (double*)malloc(sizeof(double) * ENERGY_TERMS_COUNT * surf_blocks);
		if (host_surf_terms == nullptr) {
			fprintf(stderr, "malloc failed for surface energy partials.\n");
			exit(1);
		}
		status = cudaMemcpy(host_surf_terms, d_surf_block_terms, sizeof(double) * ENERGY_TERMS_COUNT * surf_blocks, cudaMemcpyDeviceToHost);
		if (status != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed for surface energy partials: %s\n", cudaGetErrorString(status));
			exit(1);
		}
		for (unsigned int block = 0; block < surf_blocks; block++) {
			for (int term = 0; term < ENERGY_TERMS_COUNT; term++) {
				host_energy_terms[term] += host_surf_terms[block * ENERGY_TERMS_COUNT + term];
			}
		}
		free(host_surf_terms);
	}

	en_ldg[0] = host_energy_terms[LDG_TOTAL] * dV;
	en_ldg[1] = host_energy_terms[LDG_IN] * dVi;
	en_ldg[2] = host_energy_terms[LDG_OUT] * dVo;
	en_surf[0] = host_energy_terms[SURF_DROP];
	en_surf[1] = host_energy_terms[SURF_PART];

	for (int i = 0; i < 5; i++) {
		en_el[i] = 0.;
		en_el_in[i] = 0.;
		en_el_out[i] = 0.;
	}

	en_el[0] = 0.5 * dV * L1 * host_energy_terms[EL_TOT_L1];
	en_el[1] = 0.5 * dV * L2 * host_energy_terms[EL_TOT_L2];
	en_el[4] = dV * (double)chiral * 2. * L1 * qch * host_energy_terms[EL_TOT_CHIRAL];

	en_el_in[0] = 0.5 * dVi * L1 * host_energy_terms[EL_IN_L1];
	en_el_in[1] = 0.5 * dVi * L2 * host_energy_terms[EL_IN_L2];
	en_el_in[4] = dVi * (double)chiral * 2. * L1 * qch * host_energy_terms[EL_IN_CHIRAL];

	en_el_out[0] = 0.5 * dVo * L1 * host_energy_terms[EL_OUT_L1];
	en_el_out[1] = 0.5 * dVo * L2 * host_energy_terms[EL_OUT_L2];
	en_el_out[4] = dVo * (double)chiral * 2. * L1 * qch * host_energy_terms[EL_OUT_CHIRAL];

	en_tot = en_ldg[0] + en_el[0] + en_el[1] + en_el[2] + en_el[3] + en_el[4] + en_surf[0] + en_surf[1];
	dE = en_tot - old_en;
	old_en = en_tot;

#if DEBUG_ENERGY_COMPARE
	if (cycle % check_every == 0) {
		debug_compare_cpu_gpu_energy();
	}
#endif
	
	if(DoubleU){
		if(cycle % check_every == 0){
			printf("LdG: %lf, En_L1: %lf", en_ldg[0], en_el[0]);
			if(L2 != 0){
				printf(", L2: %lf", en_el[1]);
			}
			if(chiral != 0){
				printf(", Chiral: %lf", en_el[4]);
			}
			printf(", En_Surf1: %lf, Par_Surf: %lf, Cycle: %d", en_surf[0], en_surf[1], cycle);
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

				if (h_bulktype[i] == 1 || h_bulktype[i] == 13){
					ldg_ans[0] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
					ldg_ans[1] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
				}
				else if (h_bulktype[i] == 2 || h_bulktype[i] == 3){

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
			
				ldg_ans[0] += 0.5 * (1. - U / 3.) * trace2 - U / 3. * trace3 + U * 0.25 * trace2 * trace2;
				
			}
		}
		ldg_ans[0] *= dV;
	}
}

void elastic_energy(double ans[5], double ans_in[5], double ans_out[5]) {
	
	double dQ[3][6] = {{0.}};
	double Qin[6] = { 0 };
	double vec[3] = { 0. };

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
				if(h_bulktype[i] == 1 || h_bulktype[i] == 13){
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
					if(h_bulktype[i] == 1 || h_bulktype[i] == 13){
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
					if(h_bulktype[i] == 1 || h_bulktype[i] == 13){
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
					if(h_bulktype[i] == 1 || h_bulktype[i] == 13){
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
					if(h_bulktype[i] == 1 || h_bulktype[i] == 13){
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
		if ((signal[i] >= 2 && signal[i] <= 8) || signal[i] == 12 || signal[i] == 13 || (signal[i] >= 20 && signal[i] <= 23)) {
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
			else if (signal[i] == 20 || signal[i] == 21) {
				degen = 1;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if (signal[i] == 22 || signal[i] == 23) {
				degen = 2;
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
					en_degen(Qin, loc_nu, Qdiff);

					if (npboundary) {
						ans[1] += Wstr * trqq(Qdiff) * dApart;
					}
					else {
						ans[0] += Wstr * trqq(Qdiff) * dA;
					}

				}
				else if(degen == 2){
					for (int n = 0; n < 6; n++)	Qin[n] = Qold[i * 6 + n];
					for (int n = 0; n < 3; n++)	loc_nu[n] = nu[nb * 3 + n];
					en_conic(Qin, loc_nu, Qdiff);

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

void en_degen(double* Qin, double* loc_nu, double* Qdiff){
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	double third = 1.0 / 3;
	Qtemp[0][0] = Qin[0] + third * S0;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + third * S0;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + third * S0;
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			if(i == j) ptemp[i][j] = 1 - loc_nu[i] * loc_nu[j];
			else ptemp[i][j] = - loc_nu[i] * loc_nu[j];
		}
	}
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			Qp[i][j] = 0;
			for(int l = 0; l < 3; l++){
				for(int m = 0; m < 3; m++){
					Qp[i][j] += ptemp[i][l]*Qtemp[l][m]*ptemp[m][j];
				}
			}
		}
	}
	Qdiff[0] = Qtemp[0][0] - Qp[0][0];
	Qdiff[1] = Qtemp[0][1] - Qp[0][1];
	Qdiff[2] = Qtemp[0][2] - Qp[0][2];
	Qdiff[3] = Qtemp[1][1] - Qp[1][1];
	Qdiff[4] = Qtemp[1][2] - Qp[1][2];
	Qdiff[5] = Qtemp[2][2] - Qp[2][2];
}


void en_conic(double* Qin, double* loc_nu, double* Qdiff){
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	double third = 1. / 3.;
	double cosTiltAngle;
	double cosTiltAngleSq;
	
	Qtemp[0][0] = Qin[0] + third * S0;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + third * S0;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + third * S0;

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			ptemp[i][j] = loc_nu[i] * loc_nu[j];
		}
	}

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			Qp[i][j] = 0;
			for(int l = 0; l < 3; l++){
				for(int m = 0; m < 3; m++){
					Qp[i][j] += ptemp[i][l] * Qtemp[l][m] * ptemp[m][j];
				}
			}
		}
	}
	
	cosTiltAngle = cos(tiltAngle / 180.0 * M_PI);
	cosTiltAngleSq = pow(cosTiltAngle, 2);
	
	Qdiff[0] =  Qp[0][0] - cosTiltAngleSq * S0 * ptemp[0][0];
	Qdiff[1] =  Qp[0][1] - cosTiltAngleSq * S0 * ptemp[0][1];
	Qdiff[2] =  Qp[0][2] - cosTiltAngleSq * S0 * ptemp[0][2];
	Qdiff[3] =  Qp[1][1] - cosTiltAngleSq * S0 * ptemp[1][1];
	Qdiff[4] =  Qp[1][2] - cosTiltAngleSq * S0 * ptemp[1][2];
	Qdiff[5] =  Qp[2][2] - cosTiltAngleSq * S0 * ptemp[2][2];
}
