#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>
#include <immintrin.h>
#include <string.h>

static const int thread_num = 2;
static int N, M;
unsigned int *hist;
sem_t sem;

void* hist_thread(void* a){
  unsigned char (*img)[N/thread_num][M] = (unsigned char (*)[N/thread_num][M])a;

  for(int i = 0; i < N/thread_num; ++i){
    for(int j = 0; j < M; ++j){
      unsigned int index = (*img)[i][j];
      // _xbegin();
      // ++hist[index];
      // _xend();
      __transaction_atomic{++hist[index];}
    }
  }
  return NULL;
}

int main(int argc, char *argv[]){
  sem_init(&sem, 0, 1);
  N = atoi(argv[1]);
  M = atoi(argv[2]);
  pthread_t threads[thread_num];
  unsigned char (*image)[N][M] = malloc(sizeof(char) * N * M);
  hist = (unsigned int*)malloc(sizeof(int) * 256);
  memset(hist, 0x0, sizeof(int) * 256);

  for(int i = 0; i < N; ++i)
    for(int j = 0; j < M; ++j){
      (*image)[i][j] = rand() % 256;
    }

  for(int thread_i = 0; thread_i < thread_num; ++thread_i){
    pthread_create(&threads[thread_i], NULL, hist_thread, &(*image)[thread_i * N/thread_num][0]);
  }
  for(int thread_i = 0; thread_i < thread_num; ++thread_i){
    pthread_join(threads[thread_i], NULL);
  }

  for(int i = 0; i < 256; ++i){
    printf("%d:", i);
    for(int j = 0; j < hist[i]; ++j)
      printf("|");
    printf("\n");
  }
}
