#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <cuda.h>
#include <input.h>
#include <output.h>

#define _index_macro(a, b, c) a[M * (b) + (c)]
#define _index_macro_pat(a, b, c) a[(M + 2) * (b) + (c)]

static void checkCudaCall(cudaError_t result) {
        if (result != cudaSuccess) {
                printf("cuda error \n");
                exit(1);
        }
}

__device__ static double atomicMaxf(double* address, double val)
{
        unsigned long long int* address_as_i = (unsigned long long int*) address;
        unsigned long long int old = *address_as_i, assumed;
        do {
                assumed = old;
                old = atomicCAS(address_as_i, assumed,
                                __double_as_longlong(fmaxf(val,__longlong_as_double(assumed))));
        } while (assumed != old);
        return __longlong_as_double(old);
}


__global__ void iteration(double* temp_init, double* temp_tmp, double* conductivity, int N, int M, double* maxdiff) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned j = i % M + 1;
        double dir_nc = sqrtf(2)/(sqrtf(2) + 1) / 4;
        double dig_nc = 1 /(sqrtf(2) + 1) / 4;
        double diff;
        i = i / M + 1;
        if(i > N)
                return;

        *maxdiff = 0;
        for(int i = 0; i < N; ++i) {
                _index_macro_pat(temp_init, i + 1, 0) = _index_macro_pat(temp_init, i + 1, M); // move last column to 0's
                _index_macro_pat(temp_init, i + 1,M + 1) = _index_macro_pat(temp_init, i + 1, 1); // move first column to (M+1)'s
                _index_macro_pat(temp_tmp, i + 1,0) = _index_macro_pat(temp_init, i + 1, M); // move last column to 0's
                _index_macro_pat(temp_tmp, i + 1, M + 1) = _index_macro_pat(temp_init, i + 1, 1); // move first column to (M+1)'s
        }
        double weighted_neighb = dir_nc *
                                 ( // Direct neighbors
                _index_macro_pat(temp_init, i - 1, j) + _index_macro_pat(temp_init, i, j - 1) +
                _index_macro_pat(temp_init, i + 1, j) + _index_macro_pat(temp_init, i, j + 1)
                                 ) + dig_nc *
                                 ( // Diagonal neighbors
                _index_macro_pat(temp_init, i - 1, j - 1) +_index_macro_pat(temp_init, i + 1, j - 1) +
                _index_macro_pat(temp_init, i - 1, j + 1) + _index_macro_pat(temp_init, i + 1, j + 1)
                                 );
        weighted_neighb *= (1 - _index_macro(conductivity, i - 1, j - 1));
        _index_macro_pat(temp_tmp, i, j) = _index_macro_pat(temp_init, i, j) * _index_macro(conductivity, i - 1, j - 1);
        _index_macro_pat(temp_tmp, i, j) += weighted_neighb;

        diff = fabs(_index_macro_pat(temp_tmp, i, j) - _index_macro_pat(temp_init, i, j));
        atomicMaxf(maxdiff, diff);
}
extern "C"
void cuda_do_compute(const struct parameters* p, struct results *r) {
        struct timeval before, after;
        int threadBlockSize = 512;
        double (*temp_init)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2])malloc((p->N + 2) * (p->M + 2) * sizeof(double));
        double maxdiff = p->threshold + 1.;
        unsigned int iter = 0;
        double* deviceMaxdiff;
        double* deviceA;
        double* deviceB;
        double* deviceConductivity;
        double local_sum = 0;
        int N = p->N, M = p->M;
        cudaMalloc((void**)&deviceMaxdiff, sizeof(uint64_t));
        cudaMemset(deviceMaxdiff, 0, sizeof(double));

        checkCudaCall(cudaMalloc((void **) &deviceA, ((p->N + 2) * (p->M + 2) * sizeof(double) )));
        checkCudaCall(cudaMalloc((void **) &deviceB, ((p->N + 2) * (p->M + 2) * sizeof(double) )));
        checkCudaCall(cudaMalloc((void **) &deviceConductivity, ((p->N) * (p->M) * sizeof(double) )));

        for(int i = 0; i < N; ++i)
                memcpy(&(*temp_init)[i + 1][1], &_index_macro(p->tinit, i, 0), M * sizeof(double));

        for(int i = 0; i < M; ++i) {
                (*temp_init)[0][i + 1] = (*temp_init)[1][i + 1];
                (*temp_init)[N + 1][i + 1] = (*temp_init)[N][i + 1];
        }
// Filling [0][0], [0][M + 1], [N + 1][0], [N + 1][M + 1] elems
        (*temp_init)[0][0] = (*temp_init)[0][M];
        (*temp_init)[0][M + 1] = (*temp_init)[0][1];
        (*temp_init)[N + 1][0] = (*temp_init)[N][M];
        (*temp_init)[N + 1][M + 1] = (*temp_init)[N + 1][1];


// copy the original vectors to the GPU
        checkCudaCall(cudaMemcpy(deviceA, temp_init, (p->N + 2) * (p->M + 2) * sizeof(double), cudaMemcpyHostToDevice));
        checkCudaCall(cudaMemcpy(deviceB, temp_init, (p->N + 2) * (p->M + 2) * sizeof(double), cudaMemcpyHostToDevice));
        checkCudaCall(cudaMemcpy(deviceConductivity, p->conductivity, (p->N) * (p->M) * sizeof(double), cudaMemcpyHostToDevice));


        gettimeofday(&before, NULL);

        while(iter++ < (p->maxiter) && maxdiff > p->threshold) {
// execute kernel
                maxdiff = 0;
                {double* tmp = deviceA; deviceA = deviceB; deviceB = tmp;}
                dim3 gridSize(N/threadBlockSize + 1, M / threadBlockSize + 1);
                iteration<<<(N * M / threadBlockSize + 1), threadBlockSize>>>(deviceA, deviceB, deviceConductivity, N, M, deviceMaxdiff);
                cudaMemcpy(&maxdiff, deviceMaxdiff, sizeof(uint64_t), cudaMemcpyDeviceToHost);
                if((iter % p->period) == 0){
            			local_sum = 0;
            			gettimeofday(&end, 0);
            			r->tmin = r->tmax = (*temp_tmp)[1][1];
            			for(int i = 1; i <= N; ++i){
            				for(int j = 1; j <= M ; ++j){
            					if((*temp_init)[i][j] > r->tmax)
            						r->tmax = (*temp_init)[i][j];
            					if((*temp_init)[i][j] < r->tmin)
            						r->tmin = (*temp_init)[i][j];
            					local_sum += (*temp_init)[i][j];
            				}
            			}
        }
        gettimeofday(&after, NULL);

// check whether the kernel invocation was successful
        checkCudaCall(cudaGetLastError());

// copy result back
        checkCudaCall(cudaMemcpy(temp_init, deviceB, (N + 2) * (M + 2) * sizeof(double), cudaMemcpyDeviceToHost));
        r->tmin = r->tmax = (*temp_init)[1][1];
        local_sum = 0;
        for(int i = 1; i <= N; ++i) {
                for(int j = 1; j <= M; ++j) {
                        if((*temp_init)[i][j] > r->tmax)
                                r->tmax = (*temp_init)[i][j];
                        if((*temp_init)[i][j] < r->tmin)
                                r->tmin = (*temp_init)[i][j];
                        local_sum += (*temp_init)[i][j];
                }
        }
        r->niter = iter - 1;
        r->tavg = local_sum /(N * M);
        r->maxdiff = maxdiff;
        r->time = (double)(after.tv_sec - before.tv_sec) +
                  (double)(after.tv_usec - before.tv_usec) / 1e6;
        checkCudaCall(cudaFree(deviceA));
        checkCudaCall(cudaFree(deviceB));
        checkCudaCall(cudaFree(deviceConductivity));
        checkCudaCall(cudaFree(deviceMaxdiff));

}
