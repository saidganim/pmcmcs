#define _GNU_SOURCE

#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include <string.h>
#include <output.h>
#include <pthread.h>

static unsigned int N,M;

#define PTHREAD_NUM 2

#define _index_macro(a, b, c) a[M * (b) + (c)]

pthread_spinlock_t maxdiff_lock;

typedef struct params{
  pthread_barrier_t *barriers;
  double (*conductivity);
  double (*temp_tmp);
  double (*temp_init);
  double* maxdiff;
  double threshold;
  int maxiter;
  int period;
  int* iter;
  unsigned int ptid;
  unsigned int N;
  unsigned int M;
} params_t;

void* run_job(void* a){
  params_t *p = (params_t*)a;
  double (*temp_tmp)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2])p->temp_tmp;
  double (*temp_init)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2])p->temp_init;
  double (*conductivity)[p->N][p->M] = (double (*)[p->N][p->M])p->conductivity;
  double maxdiff = *p->maxdiff;
	unsigned int iter = 0;
	double dir_nc = sqrt(2)/(sqrt(2) + 1) / 4;
	double dig_nc = 1 /(sqrt(2) + 1) / 4;
  printf("PTHREAD %d WORK CHUNK IS (%d, %f)\n", p->ptid,  1 + p->ptid * N/PTHREAD_NUM, p->ptid * ((double)N)/PTHREAD_NUM + ((double)N) /PTHREAD_NUM);
  while(iter++ < p->maxiter && *p->maxdiff > p->threshold){
    // do computations;
		// update most left and most right columns( cache suffers )
		for(int i = p->ptid * N/PTHREAD_NUM; i < p->ptid * ((double)N)/PTHREAD_NUM + ((double)N) /PTHREAD_NUM; ++i){
			(*temp_init)[i + 1][0] = (*temp_init)[i + 1][M]; // move last column to 0's
			(*temp_init)[i + 1][M + 1] = (*temp_init)[i + 1][1]; // move first column to (M+1)'s
			(*temp_tmp)[i + 1][0] = (*temp_init)[i + 1][M]; // move last column to 0's
			(*temp_tmp)[i + 1][M + 1] = (*temp_init)[i + 1][1]; // move first column to (M+1)'s
		}
    if(pthread_barrier_wait(&p->barriers[0]) == PTHREAD_BARRIER_SERIAL_THREAD){
      pthread_barrier_destroy(&p->barriers[0]);
      pthread_barrier_init(&p->barriers[0], NULL, PTHREAD_NUM);
    }
    maxdiff = 0.;
    *p->maxdiff = maxdiff;

		// finally start computations
		for(int i = 1 + p->ptid * N/PTHREAD_NUM; i <= p->ptid * ((double)N)/PTHREAD_NUM + ((double)N) /PTHREAD_NUM; ++i){
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
				(*temp_tmp)[i][j] = (*temp_init)[i][j] * _index_macro(p->conductivity, i - 1, j - 1);
				(*temp_tmp)[i][j] += weighted_neighb;
				if(fabs((*temp_init)[i][j] - (*temp_tmp)[i][j]) >  maxdiff)
					maxdiff = fabs((*temp_init)[i][j] - (*temp_tmp)[i][j]);
			}
		}
		// syncrhonizing init matrix and temporary one
		double* tmp_tmp = temp_init;
		temp_init = temp_tmp;
		temp_tmp = tmp_tmp;
    {
      pthread_spin_lock(&maxdiff_lock);
      if(maxdiff > *p->maxdiff)
        *p->maxdiff = maxdiff;
      pthread_spin_unlock(&maxdiff_lock);
    }
    // if(
    // }
    if(pthread_barrier_wait(&p->barriers[1]) == PTHREAD_BARRIER_SERIAL_THREAD){
      pthread_barrier_destroy(&p->barriers[1]);
      pthread_barrier_init(&p->barriers[1], NULL, PTHREAD_NUM);
    }
		// if((iter % p->period) == 0){
		// 	double local_sum = 0;
		// 	gettimeofday(&end, 0);
		// 	r->tmin = r->tmax = (*temp_tmp)[1][1];
		// 	for(int i = 1; i <= N; ++i){
		// 		for(int j = 1; j <= M ; ++j){
		// 			if((*temp_init)[i][j] > r->tmax)
		// 				r->tmax = (*temp_init)[i][j];
		// 			if((*temp_init)[i][j] < r->tmin)
		// 				r->tmin = (*temp_init)[i][j];
		// 			local_sum += (*temp_init)[i][j];
		// 		}
		// 	}
		// 	r->time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
		// 	r->niter = iter;
		// 	r->tavg = local_sum /(N * M);
		// 	r->maxdiff = maxdiff;
		// 	report_results(p, r);
		// }
	}
  *p->iter = iter;
  return NULL;
}

void do_compute(const struct parameters* p, struct results *r)
{
	struct timeval start, end;
	double (*temp_tmp)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2]) malloc((p->N + 2) * (p->M + 2) * sizeof(double));
	double (*temp_init)[p->N + 2][p->M + 2] = (double (*)[p->N + 2][p->M + 2]) malloc((p->N + 2) * (p->M + 2) * sizeof(double));
	double maxdiff = p->threshold + 1.;
	double sum = 0;
  unsigned int iter = 0;
  pthread_t threads[PTHREAD_NUM];
  pthread_barrier_t barriers[2];
  params_t params[PTHREAD_NUM];
	N = p->N;
	M = p->M;
  //
  // pthread_attr_t ptattr;
  // pthread_Attr_init(&*ptattr);
  //

	//copying init matrix
	for(int i = 0; i < N; ++i)
		memcpy(&(*temp_init)[i + 1][1], &_index_macro(p->tinit, i, 0), M * sizeof(double));
	// copying top and bottom rows to new rows
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

  pthread_barrier_init(barriers, NULL, PTHREAD_NUM);
  pthread_barrier_init(barriers + 1, NULL, PTHREAD_NUM);
  pthread_spin_init(&maxdiff_lock, 0);
	gettimeofday(&start, 0);
  for( int i = 0; i < PTHREAD_NUM; ++i){
    params[i].barriers = barriers;
    params[i].temp_tmp = temp_tmp;
    params[i].temp_init = temp_init;
    params[i].conductivity = p->conductivity;
    params[i].maxdiff = &maxdiff;
    params[i].threshold =p->threshold;
    params[i].maxiter = p->maxiter;
    params[i].period = p->period;
    params[i].iter = &iter;
    params[i].ptid = i;
    params[i].N = p->N;
    params[i].M = p->M;
    pthread_create(&threads[i], NULL, run_job, &params[i]);
  }


  for( int i = 0; i < PTHREAD_NUM; ++i){
    void* res = NULL;
    pthread_join(threads[i], &res);
  }

	gettimeofday(&end, 0);
	r->tmin = r->tmax = (*temp_tmp)[1][1];
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
	free(temp_tmp);
	free(temp_init);

	return;
}
