#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include <string.h>
#include <output.h>

static unsigned int N,M;

#define _index_macro(a, b, c) a[M * (b) + (c)]
#define max(a,b) \
({\
	__typeof(a) _a = (a);\
	__typeof(b) _b = (b);\
	(_a) > (_b) ? (_a) : (_b); })

void do_compute(const struct parameters* p, struct results *r)
{
	struct timeval start, end;
	double (* restrict temp_tmp)[p->N][p->M] = (double (*)[p->N][p->M]) malloc((p->N) * (p->M) * sizeof(double));
	double (* restrict temp_init)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2]) malloc((p->N + 2) * (p->M + 2) * sizeof(double));
	temp_tmp = (__typeof(temp_tmp)) (((unsigned long)temp_tmp + 15) & ~0xF);
	temp_init = (__typeof(temp_init)) (((unsigned long)temp_init + 15) & ~0xF);
	temp_init = __builtin_assume_aligned(temp_init, 16);
	temp_tmp = __builtin_assume_aligned(temp_tmp, 16);
	double maxdiff = p->threshold + 1.;
	double sum = 0.;
	unsigned int iter = 0;
	double dir_nc = sqrt(2)/(sqrt(2) + 1) / 4;
	double dig_nc = 1 /(sqrt(2) + 1) / 4;
	N = p->N;
	M = p->M;
	//copying init matrix
	for(int i = 0; i < N; ++i)
		memcpy(&(*temp_init)[i + 1][1], &_index_macro(p->tinit, i, 0), M * sizeof(double));
	// copying top and bottom rows to new rows
	#pragma GCC ivdep
	for(int i = 0; i < M; ++i){
		(*temp_init)[0][i + 1] = (*temp_init)[1][i + 1];
		(*temp_init)[N + 1][i + 1] = (*temp_init)[N][i + 1];
	}
	// Filling [0][0], [0][M + 1], [N + 1][0], [N + 1][M + 1] elems
	(*temp_init)[0][0] = (*temp_init)[0][M];
	(*temp_init)[0][M + 1] = (*temp_init)[0][1];
	(*temp_init)[N + 1][0] = (*temp_init)[N][M];
	(*temp_init)[N + 1][M + 1] = (*temp_init)[N + 1][1];

	gettimeofday(&start, 0);
	while(iter++ < p->maxiter && maxdiff > p->threshold){
		// do computations;
		maxdiff = 0;
		// update most left and most right columns( cache suffers )
		for(int i = 0; i < N; ++i){
			(*temp_init)[i + 1][0] = (*temp_init)[i + 1][M]; // move last column to 0's
			(*temp_init)[i + 1][M + 1] = (*temp_init)[i + 1][1]; // move first column to (M+1)'s
		}

		// finally start computations
		for(int i = 1; i <= N; ++i){
			#pragma clang loop vectorize_width(4)
			#pragma clang loop interleave_count(4)
			for(int j = 1; j <= M ; ++j){
				double weighted_neighb = dir_nc *
				( // Direct neighbors
					(*temp_init)[i - 1][j] + (*temp_init)[i][j - 1] +
					(*temp_init)[i + 1][j] + (*temp_init)[i][j + 1]
				) + dig_nc *
				( // Diagonal neighbors
					(*temp_init)[i - 1][j - 1] + (*temp_init)[i + 1][j - 1] +
					(*temp_init)[i - 1][j + 1] + (*temp_init)[i + 1][j + 1]
				);
				weighted_neighb *= (1 - _index_macro(p->conductivity, i - 1, j - 1));
				(*temp_tmp)[i - 1][j - 1] = (*temp_init)[i][j] * _index_macro(p->conductivity, i - 1, j - 1);
				(*temp_tmp)[i - 1][j - 1] += weighted_neighb;
				maxdiff = max(maxdiff, fabs((*temp_init)[i][j] - (*temp_tmp)[i - 1][j - 1]));

			}
		}
		// syncrhonizing init matrix and temporary one
		for(int i = 0; i < N; ++i)
			memcpy( &(*temp_init)[i + 1][1], &(*temp_tmp)[i][0], M * sizeof(double));
		if((iter % p->period) == 0){
			double local_sum = 0;
			gettimeofday(&end, 0);
			r->tmin = r->tmax = (*temp_tmp)[0][0];
			for(int i = 1; i <= N; ++i){
				for(int j = 1; j <= M ; ++j){
					if((*temp_init)[i][j] > r->tmax)
						r->tmax = (*temp_init)[i][j];
					if((*temp_init)[i][j] < r->tmin)
						r->tmin = (*temp_init)[i][j];
					local_sum += (*temp_init)[i][j];
				}
			}
			r->time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
			r->niter = iter;
			r->tavg = local_sum /(N * M);
			r->maxdiff = maxdiff;
			report_results(p, r);
		}
	}

	r->tmin = r->tmax = (*temp_tmp)[0][0];
	#pragma clang loop vectorize_width(4)
	#pragma clang loop interleave_count(4)
	for(int i = 1; i <= N; ++i){
		for(int j = 1; j <= M ; ++j){
			if((*temp_init)[i][j] > r->tmax)
				r->tmax = (*temp_init)[i][j];
			if((*temp_init)[i][j] < r->tmin)
				r->tmin = (*temp_init)[i][j];
			sum += (*temp_init)[i][j];
		}
	}
	gettimeofday(&end, 0);
	r->time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
	r->niter = iter - 1;
	r->tavg = sum /(N * M);
	r->maxdiff = maxdiff;
	return;
}
