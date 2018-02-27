
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>

/* Ordering of the vector */
typedef enum Ordering {ASCENDING, DESCENDING, RANDOM} Order;

int debug = 0;

int* __mrg(int *v, long l){
    int* tmp = (int*) malloc(l * sizeof(int));
    int* tmps = tmp;
    int *left = v,  *right = v + l/2;
    int *lefts = left, *rights = right;
    for(int i = 0; i < l; ++i){
      if( right - rights == l/2)
        *(tmp++) = *(left++);
      else if (left - lefts == l/2)
        *(tmp++) = *(right++);
      else if(*left < *right)
        *(tmp++) = *(left++);
      else
        *(tmp++) = *(right++);
    }
    memcpy(v, tmps, l * sizeof(int));
    free(tmps);
    return v;
}

int* __merge(int *v, long l){
  if(l == 2){
    if(v[0] > v[1]){
      v[0] += v[1];
      v[1] = v[0] - v[1];
      v[0] = v[0] - v[1];
    }
    return v;
  } else {
    int* tmp = (int*) malloc(l * sizeof(int));
    int* tmps = tmp;
    int *left, *right;
    // #pragma omp task
      left = __merge(v, l / 2);
    // #pragma omp task
      right = __merge(v + l/2, l/2);

    // #pragma omp taskwait
    int *lefts = left, *rights = right;
    for(int i = 0; i < l; ++i){
      if( right - rights == l/2)
        *(tmp++) = *(left++);
      else if (left - lefts == l/2)
        *(tmp++) = *(right++);
      else if(*left < *right)
        *(tmp++) = *(left++);
      else
        *(tmp++) = *(right++);
    }
    memcpy(v, tmps, l * sizeof(int));
    free(tmps);
    return v;
  }

}

/* Sort vector v of l elements using mergesort */
void msort(int *v, long l){
  unsigned int thread_num = omp_get_max_threads();
  // #pragma omp parallel
  {
    // #pragma omp single
    {
      for(int i = 0; i < thread_num; ++i){
        __merge(v + i * (l / thread_num), l / thread_num);
      }
    }

  }

  __mrg(v, l);
  // __merge(v, l);

}

void print_v(int *v, long l) {
  printf("\n");
  for(long i = 0; i < l; i++) {
    if(i != 0 && (i % 10 == 0)) {
      printf("\n");
    }
    printf("%d ", v[i]);
  }
  printf("\n");
}

int main(int argc, char **argv) {

  int c;
  int seed = 42;
  long length = 1e4 * 1e4;
  Order order = DESCENDING;
  int *vector;
  struct timeval start, end;
  double time;

  /* Read command-line options. */
  while((c = getopt(argc, argv, "adrgl:s:")) != -1) {
    switch(c) {
      case 'a':
        order = ASCENDING;
        break;
      case 'd':
        order = DESCENDING;
        break;
      case 'r':
        order = RANDOM;
        break;
      case 'l':
        length = atol(optarg);
        break;
      case 'g':
        debug = 1;
        break;
      case 's':
        seed = atoi(optarg);
	break;
      case '?':
        if(optopt == 'l' || optopt == 's') {
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        }
        else if(isprint(optopt)) {
          fprintf(stderr, "Unknown option '-%c'.\n", optopt);
        }
        else {
          fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
        }
        return -1;
      default:
        return -1;
      }
  }

  /* Seed such that we can always reproduce the same random vector */
  srand(seed);

  /* Allocate vector. */
  vector = (int*)malloc(length*sizeof(int));
  if(vector == NULL) {
    fprintf(stderr, "Malloc failed...\n");
    return -1;
  }

  /* Fill vector. */
  switch(order){
    case ASCENDING:
      for(long i = 0; i < length; i++) {
        vector[i] = (int)i;
      }
      break;
    case DESCENDING:
      for(long i = 0; i < length; i++) {
        vector[i] = (int)(length - i);
      }
      break;
    case RANDOM:
      for(long i = 0; i < length; i++) {
        vector[i] = rand();
      }
      break;
  }

  if(debug) {
    print_v(vector, length);
  }

  /* Sort */
  gettimeofday(&start, 0);
  msort(vector, length);
  gettimeofday(&end, 0);
  time = (end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0));
  printf("TIME SPENT ON SORTING IS %f\n", time);


  if(debug) {
    print_v(vector, length);
  }

  return 0;
}
