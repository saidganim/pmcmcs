#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include <string.h>
#include <output.h>
#include <smmintrin.h>

static unsigned int N,M;

#define _index_macro(a, b, c) a[M * (b) + (c)]

void do_compute(const struct parameters* p, struct results *r)
{
	struct timeval start, end;
	double (*temp_tmp)[p->N][p->M] = (double (*)[p->N][p->M]) malloc((p->N) * (p->M) * sizeof(double) + 15);
	temp_tmp = (__typeof(temp_tmp)) ((unsigned long)temp_tmp & ~0xF);
	double (*temp_init)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2]) malloc((p->N + 2) * (p->M + 2) * sizeof(double) + 15);
	temp_init = (__typeof(temp_init)) ((unsigned long)temp_init & ~0xF);
	double maxdiff = p->threshold + 1.;
	double sum = 0.;
	unsigned int iter = 0;
	double dir_nc = sqrt(2)/(sqrt(2) + 1) / 4;
	double dig_nc = 1 /(sqrt(2) + 1) / 4;
	__m128d weighted_neighb_reg;
	__m128d temp_init_reg;
	__m128d temp_init_reg_left;
	__m128d temp_init_reg_right;
	__m128d temp_init_reg_top;
	__m128d temp_init_reg_d1;
	__m128d temp_init_reg_d2;
	__m128d temp_init_reg_d3;
	__m128d temp_init_reg_d4;
	__m128d temp_init_reg_bottom;
	__m128d temp_tmp_reg;
	__m128d tmp_reg;
	__m128d direct_sum;
	__m128d diag_sum;
	__m128d const4_reg;
	__m128d const1_reg;
	__m128d coef1_reg;
	__m128d coef2_reg;
	__m128d cond_reg;
	__m128d rev_cond_reg;
	//const4_reg = _mm128_set_pd(1./4,1./4,1./4,1./4);
	const1_reg = _mm_set_pd(1.,1.);
	coef1_reg = _mm_set_pd(dir_nc, dir_nc);
	coef2_reg = _mm_set_pd(dig_nc, dig_nc);
	N = p->N;
	M = p->M;
	//copying init matrix
	for(int i = 0; i < N; ++i)
		memcpy(&(*temp_init)[i + 1][1], &_index_macro(p->tinit, i, 0), M * sizeof(double));
	// copying top and bottom rows to new rows
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
		int elems_per_iter = (128 / sizeof(double) / 8);
		for(int i = 1; i <= N; ++i){
			for(int j = 1; j <= (M - M % elems_per_iter) ; j += elems_per_iter){
				cond_reg = _mm_load_pd(&_index_macro(p->conductivity, i - 1, j - 1));
				rev_cond_reg = _mm_sub_pd(const1_reg, cond_reg);
				temp_init_reg = _mm_loadu_pd(&(*temp_init)[i][j]);
				temp_init_reg_left = _mm_loadu_pd(&(*temp_init)[i][j - 1]);
				temp_init_reg_right = _mm_loadu_pd(&(*temp_init)[i][j + 1]);
				temp_init_reg_top = _mm_loadu_pd(&(*temp_init)[i - 1][j]);
				temp_init_reg_bottom = _mm_loadu_pd(&(*temp_init)[i + 1][j]);

				temp_init_reg_d1 = _mm_loadu_pd(&(*temp_init)[i - 1][j - 1]);
				temp_init_reg_d2 = _mm_loadu_pd(&(*temp_init)[i + 1][j - 1]);
				temp_init_reg_d3 = _mm_loadu_pd(&(*temp_init)[i - 1][j + 1]);
				temp_init_reg_d4 = _mm_loadu_pd(&(*temp_init)[i + 1][j + 1]);

				direct_sum = _mm_add_pd(temp_init_reg_left, temp_init_reg_right);
				direct_sum = _mm_add_pd(direct_sum, temp_init_reg_top);
				direct_sum = _mm_add_pd(direct_sum, temp_init_reg_bottom);
				direct_sum = _mm_mul_pd(direct_sum, coef1_reg);

				diag_sum = _mm_add_pd(temp_init_reg_d1, temp_init_reg_d2);
				diag_sum = _mm_add_pd(diag_sum, temp_init_reg_d3);
				diag_sum = _mm_add_pd(diag_sum, temp_init_reg_d4);
				diag_sum = _mm_mul_pd(diag_sum, coef2_reg);

				diag_sum = _mm_add_pd(diag_sum, direct_sum);
				diag_sum = _mm_mul_pd(diag_sum, rev_cond_reg);
				temp_init_reg = _mm_mul_pd(temp_init_reg, cond_reg);
				temp_init_reg = _mm_add_pd(temp_init_reg, diag_sum);
				_mm_storeu_pd(&(*temp_tmp)[i - 1][j - 1], temp_init_reg);
				for(int maxdiff_i = 0; maxdiff_i < elems_per_iter; ++maxdiff_i)
					if(fabs((*temp_tmp)[i - 1][j - 1 + maxdiff_i] - (*temp_init)[i][j + maxdiff_i]) >  maxdiff)
						maxdiff = fabs((*temp_tmp)[i - 1][j - 1 + maxdiff_i] - (*temp_init)[i][j + maxdiff_i]);
			}


			for(int j = (M - M % elems_per_iter) + 1; j <= M ; ++j){
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
				if(fabs((*temp_init)[i][j] - (*temp_tmp)[i - 1][j - 1]) >  maxdiff)
					maxdiff = fabs((*temp_init)[i][j] - (*temp_tmp)[i - 1][j - 1]);

			}
		}
		// syncrhonizing init matrix and temporary one
		for(int i = 0; i < N; ++i)
			memcpy( &(*temp_init)[i + 1][1], &(*temp_tmp)[i][0], M * sizeof(double));
	}

	gettimeofday(&end, 0);
	r->tmin = r->tmax = (*temp_tmp)[0][0];
	for(int i = 1; i <= N; ++i){
		for(int j = 1; j <= M ; ++j){
			if((*temp_init)[i][j] > r->tmax)
				r->tmax = (*temp_init)[i][j];
			if((*temp_init)[i][j] < r->tmin)
				r->tmin = (*temp_init)[i][j];
			sum += (*temp_init)[i][j];
		}
	}
	r->time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
	r->niter = iter - 1;
	r->tavg = sum /(N * M);
	r->maxdiff = maxdiff;
	return;
}
