#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include <string.h>
#include <output.h>
#include <immintrin.h>
#include <omp.h>

static unsigned int N,M;

#define max(a,b) \
({\
	__typeof(a) _a = (a);\
	__typeof(b) _b = (b);\
	(_a) > (_b) ? (_a) : (_b); })

#define min(a,b) \
({\
	__typeof(a) _a = (a);\
	__typeof(b) _b = (b);\
	(_a) < (_b) ? (_a) : (_b); })

#define _index_macro(a, b, c) a[M * (b) + (c)]

static inline __m256d max_from_vect(__m256d array[], int n){
		__m256d res = array[0];
		for(int i = 1; i < n; ++i)
			res = _mm256_max_pd(res, array[i]);
		return res;
}

void do_compute(const struct parameters* p, struct results *r)
{
	struct timeval start, end;
	double (*temp_tmp)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2]) malloc((p->N + 2) * (p->M + 2) * sizeof(double) + 15);
	double (*temp_init)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2]) malloc((p->N + 2) * (p->M + 2) * sizeof(double) + 15);
	temp_tmp = (__typeof(temp_tmp)) (((unsigned long)temp_tmp + 15) & ~0xF);
	temp_init = (__typeof(temp_init)) (((unsigned long)temp_init + 15) & ~0xF);
	double maxdiff = p->threshold + 1.;
	double sum = 0.;
	unsigned int iter = 0;
	double dir_nc = sqrt(2)/(sqrt(2) + 1) / 4;
	double dig_nc = 1 /(sqrt(2) + 1) / 4;
	int elems_per_iter = (256 / sizeof(double) / 8);
	double maxdiff_vect_mem[elems_per_iter / 2];
	int thread_num = 2; //omp_get_max_threads() / 2;
	__m256d weighted_neighb_reg;
	__m256d temp_init_reg;
	__m256d temp_init_reg_clone;
	__m256d maxdiff_vect[thread_num];
	__m256d maxdiff_vect2;
	__m128d tmp_reg;
	__m256d direct_sum;
	__m256d diag_sum;
	__m256d cond_reg;
	__m256d rev_cond_reg;
	__m256d double_temp_reg;
	__m256d const1_reg;
	__m256d coef1_reg;
	__m256d coef2_reg;
	omp_set_num_threads(thread_num);
	const1_reg = _mm256_set1_pd(1.);
	coef1_reg = _mm256_set1_pd(dir_nc);
	coef2_reg = _mm256_set1_pd(dig_nc);
	N = p->N;
	M = p->M;
	//copying init matrix
	//#pragma omp parallel for
	for(int i = 0; i < N; ++i)
		memcpy(&(*temp_init)[i + 1][1], &_index_macro(p->tinit, i, 0), M * sizeof(double));
	// copying top and bottom rows to new rows
	#pragma omp parallel for
	for(int i = 0; i < M; ++i){
		(*temp_init)[0][i + 1] = (*temp_init)[1][i + 1];
		(*temp_init)[N + 1][i + 1] = (*temp_init)[N][i + 1];
		(*temp_tmp)[0][i + 1] = (*temp_init)[1][i + 1];
		(*temp_tmp)[N + 1][i + 1] = (*temp_init)[N][i + 1];
	}
	// Filling [0][0], [0][M + 1], [N + 1][0], [N + 1][M + 1] elems
	(*temp_tmp)[0][0] = (*temp_init)[0][0] = (*temp_init)[0][M];
	(*temp_tmp)[0][M + 1] = (*temp_init)[0][M + 1] = (*temp_init)[0][1];
	(*temp_tmp)[N + 1][0] = (*temp_init)[N + 1][0] = (*temp_init)[N][M];
	(*temp_tmp)[N + 1][M + 1] = (*temp_init)[N + 1][M + 1] = (*temp_init)[N + 1][1];

	gettimeofday(&start, 0);
	#pragma omp parallel private(\
		weighted_neighb_reg, \
		temp_init_reg, \
		temp_init_reg_clone, \
		tmp_reg, \
		direct_sum, \
		diag_sum, \
		cond_reg, \
		rev_cond_reg, \
		double_temp_reg\
	)
	{
		__m256d temp_init_reg_left;
		__m256d temp_init_reg_right;
		__m256d temp_init_reg_top;
		__m256d temp_init_reg_d1;
		__m256d temp_init_reg_d2;
		__m256d temp_init_reg_d3;
		__m256d temp_init_reg_d4;
		__m256d temp_init_reg_bottom;
		__m256d localmaxdiff;


		while(iter < p->maxiter && maxdiff > p->threshold){
			int thread_id = omp_get_thread_num();
			// do computations;
			#pragma omp barrier
			maxdiff = 0.0;
			localmaxdiff = _mm256_set1_pd(0.0);
			// update most left and most right columns( cache suffers )
			#pragma omp for
			for(int i = 0; i < N; ++i){
				(*temp_init)[i + 1][0] = (*temp_init)[i + 1][M]; // move last column to 0's
				(*temp_init)[i + 1][M + 1] = (*temp_init)[i + 1][1]; // move first column to (M+1)'s
				(*temp_tmp)[i + 1][0] = (*temp_init)[i + 1][M]; // move last column to 0's
				(*temp_tmp)[i + 1][M + 1] = (*temp_init)[i + 1][1]; // move first column to (M+1)'s
			}

			// finally start computations
			#pragma omp for reduction (max: maxdiff)
			for(int i = 1; i <= N; ++i){
				int thread_id = omp_get_thread_num();
				for(int j = 1; j <= (M - M % elems_per_iter) ; j += elems_per_iter){
					cond_reg = _mm256_loadu_pd(&_index_macro(p->conductivity, i - 1, j - 1));
					rev_cond_reg = _mm256_sub_pd(const1_reg, cond_reg);

					temp_init_reg = _mm256_loadu_pd(&(*temp_init)[i][j]);
					temp_init_reg_clone = temp_init_reg; //_mm256_loadu_pd(&(*temp_init)[i][j]);
					temp_init_reg_left = _mm256_loadu_pd(&(*temp_init)[i][j - 1]);
					temp_init_reg_right = _mm256_loadu_pd(&(*temp_init)[i][j + 1]);
					temp_init_reg_top = _mm256_loadu_pd(&(*temp_init)[i - 1][j]);

					temp_init_reg_d1 = _mm256_loadu_pd(&(*temp_init)[i - 1][j - 1]);
					temp_init_reg_d3 = _mm256_loadu_pd(&(*temp_init)[i - 1][j + 1]);

					temp_init_reg_bottom = _mm256_loadu_pd(&(*temp_init)[i + 1][j]);

					temp_init_reg_d4 = _mm256_loadu_pd(&(*temp_init)[i + 1][j + 1]);
					temp_init_reg_d2 = _mm256_loadu_pd(&(*temp_init)[i + 1][j - 1]);

					direct_sum = _mm256_add_pd(temp_init_reg_left, temp_init_reg_right);
					direct_sum = _mm256_add_pd(direct_sum, temp_init_reg_top);
					direct_sum = _mm256_add_pd(direct_sum, temp_init_reg_bottom);
					direct_sum = _mm256_mul_pd(direct_sum, coef1_reg);

					diag_sum = _mm256_add_pd(temp_init_reg_d1, temp_init_reg_d2);
					diag_sum = _mm256_add_pd(diag_sum, temp_init_reg_d3);
					diag_sum = _mm256_add_pd(diag_sum, temp_init_reg_d4);
					diag_sum = _mm256_mul_pd(diag_sum, coef2_reg);

					diag_sum = _mm256_add_pd(diag_sum, direct_sum);
					diag_sum = _mm256_mul_pd(diag_sum, rev_cond_reg);
					temp_init_reg = _mm256_mul_pd(temp_init_reg, cond_reg);
					temp_init_reg = _mm256_add_pd(temp_init_reg, diag_sum);
					_mm256_storeu_pd(&(*temp_tmp)[i][j], temp_init_reg);
					double_temp_reg = _mm256_sub_pd(temp_init_reg, temp_init_reg_clone);
					localmaxdiff = _mm256_max_pd(localmaxdiff, _mm256_max_pd(_mm256_sub_pd(_mm256_set1_pd(0.0), double_temp_reg), double_temp_reg));
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
					(*temp_tmp)[i][j] = (*temp_init)[i][j] * _index_macro(p->conductivity, i - 1, j - 1);
					(*temp_tmp)[i][j] += weighted_neighb;
					if(fabs((*temp_init)[i][j] - (*temp_tmp)[i][j]) >  maxdiff)
						maxdiff = fabs((*temp_init)[i][j] - (*temp_tmp)[i][j]);
				}
			}
			// syncrhonizing init matrix and temporary one
			// for(int i = 0; i < N; ++i)
			// 	memcpy( &(*temp_init)[i + 1][1], &(*temp_tmp)[i][0], M * sizeof(double));
			maxdiff_vect[thread_id] = localmaxdiff;
			#pragma omp barrier
			#pragma omp single
			{
				maxdiff_vect2 =max_from_vect(maxdiff_vect, thread_num)  ;
				tmp_reg = _mm_max_pd(_mm256_extractf128_pd(maxdiff_vect2, 0), _mm256_extractf128_pd(maxdiff_vect2, 1));
				_mm_storeu_pd(maxdiff_vect_mem, tmp_reg);
				maxdiff = max(maxdiff, max(maxdiff_vect_mem[0], maxdiff_vect_mem[1]));
				++iter;
			}

			#pragma omp single nowait
			{
				double* tmp_tmp = temp_init;
				temp_init = temp_tmp;
				temp_tmp = tmp_tmp;
			}

				if(iter % p->period == 0){
					double local_sum = 0;
					r->tmin = r->tmax = (*temp_tmp)[1][1];
					maxdiff_vect2 = _mm256_set1_pd((*temp_tmp)[1][1]);
					temp_init_reg_clone = _mm256_set1_pd((*temp_tmp)[1][1]);
					direct_sum = _mm256_set1_pd(0);
					#pragma omp for
					for(int i = 1; i <= N; ++i){
						for(int j = 1; j <= M - (M % elems_per_iter) ; j += elems_per_iter){
							temp_init_reg =_mm256_loadu_pd(&(*temp_init)[i][j]);
							direct_sum = _mm256_add_pd(direct_sum, temp_init_reg);
							maxdiff_vect2 = _mm256_max_pd(maxdiff_vect2, temp_init_reg);
							temp_init_reg_clone = _mm256_min_pd(temp_init_reg_clone, temp_init_reg);
						}

						for(int j = M - (M % elems_per_iter) + 1; j <= M ; ++j){
							if((*temp_init)[i][j] > r->tmax)
								r->tmax = (*temp_init)[i][j];
							if((*temp_init)[i][j] < r->tmin)
								r->tmin = (*temp_init)[i][j];
							local_sum += (*temp_init)[i][j];
						}
					}

					tmp_reg = _mm_max_pd(_mm256_extractf128_pd(maxdiff_vect2, 0), _mm256_extractf128_pd(maxdiff_vect2, 1));
					_mm_storeu_pd(maxdiff_vect_mem, tmp_reg);
					r->tmax = max(r->tmax, max(maxdiff_vect_mem[0], maxdiff_vect_mem[1]));

					tmp_reg = _mm_min_pd(_mm256_extractf128_pd(temp_init_reg_clone, 0), _mm256_extractf128_pd(temp_init_reg_clone, 1));
					_mm_storeu_pd(maxdiff_vect_mem, tmp_reg);
					r->tmin = min(r->tmin, min(maxdiff_vect_mem[0], maxdiff_vect_mem[1]));

					tmp_reg = _mm_add_pd(_mm256_extractf128_pd(direct_sum, 0), _mm256_extractf128_pd(direct_sum, 1));
					_mm_storeu_pd(maxdiff_vect_mem, tmp_reg);

					gettimeofday(&end, 0);
					r->time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
					r->niter = iter;
					r->tavg = (local_sum + maxdiff_vect_mem[0] + maxdiff_vect_mem[1]) /(N * M);
					r->maxdiff = maxdiff;
					#pragma omp single nowait
					report_results(p, r);
				}

		}

	}


	r->tmin = r->tmax = (*temp_tmp)[1][1];
	maxdiff_vect2 = _mm256_set1_pd((*temp_tmp)[1][1]);
	temp_init_reg_clone = _mm256_set1_pd((*temp_tmp)[1][1]);
	direct_sum = _mm256_set1_pd(0);
	for(int i = 1; i <= N; ++i){
		for(int j = 1; j <= M - (M % elems_per_iter) ; j += elems_per_iter){
			temp_init_reg =_mm256_loadu_pd(&(*temp_init)[i][j]);
			direct_sum = _mm256_add_pd(direct_sum, temp_init_reg);
			maxdiff_vect2 = _mm256_max_pd(maxdiff_vect2, temp_init_reg);
			temp_init_reg_clone = _mm256_min_pd(temp_init_reg_clone, temp_init_reg);
		}

		for(int j = M - (M % elems_per_iter) + 1; j <= M ; ++j){
			if((*temp_init)[i][j] > r->tmax)
				r->tmax = (*temp_init)[i][j];
			if((*temp_init)[i][j] < r->tmin)
				r->tmin = (*temp_init)[i][j];
			sum += (*temp_init)[i][j];
		}
	}

	tmp_reg = _mm_max_pd(_mm256_extractf128_pd(maxdiff_vect2, 0), _mm256_extractf128_pd(maxdiff_vect2, 1));
	_mm_storeu_pd(maxdiff_vect_mem, tmp_reg);
	r->tmax = max(r->tmax, max(maxdiff_vect_mem[0], maxdiff_vect_mem[1]));

	tmp_reg = _mm_min_pd(_mm256_extractf128_pd(temp_init_reg_clone, 0), _mm256_extractf128_pd(temp_init_reg_clone, 1));
	_mm_storeu_pd(maxdiff_vect_mem, tmp_reg);
	r->tmin = min(r->tmin, min(maxdiff_vect_mem[0], maxdiff_vect_mem[1]));

	tmp_reg = _mm_add_pd(_mm256_extractf128_pd(direct_sum, 0), _mm256_extractf128_pd(direct_sum, 1));
	_mm_storeu_pd(maxdiff_vect_mem, tmp_reg);

	gettimeofday(&end, 0);
	r->time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
	r->niter = iter;
	r->tavg = (sum + maxdiff_vect_mem[0] + maxdiff_vect_mem[1]) /(N * M);
	r->maxdiff = maxdiff;
	free(temp_tmp);
	free(temp_init);
	return;
}
