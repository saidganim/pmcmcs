#include "pipesort.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

static int thread_num = 2;
pthread_spinlock_t spin;
static pthread_mutex_t mutex;
static pthread_cond_t cond;
static cond_counter = 0;

void __pipe_barrier(){
  pthread_mutex_lock(&mutex);
  ++cond_counter;
  if(cond_counter == thread_num){
    pthread_cond_broadcast(&cond);
    cond_counter = 0;
  } else{
    pthread_cond_wait(&cond, &mutex);
  }
  pthread_mutex_unlock(&mutex);
}

void __pseudo_pipe_barrier(){
  pthread_mutex_lock(&mutex);
  if(cond_counter == thread_num){
    pthread_cond_broadcast(&cond);
    cond_counter = 0;
  }
  pthread_mutex_unlock(&mutex);
}

void* pipe(void* a){
  int** param = (int**)a;
  pthread_t* successor = NULL;
  int own_value;
  int* value = malloc(sizeof(int));
  __pipe_barrier();
  own_value = **param;
  __pipe_barrier();
  while(TRUE){
    // Everything inside two barriers work with consistency memory. All writes
    // and variable changes happen outside this section
    __pipe_barrier();
      if(*param == NULL){
        __pipe_barrier();
        value = NULL;
        break;
      }
      int stable_value = **param;
    __pipe_barrier();

    if(stable_value > own_value){
      *value = own_value;
      own_value = stable_value;
    } else
      *value = stable_value;

    if(successor == NULL){
      pthread_spin_lock(&spin);
        ++thread_num;
      pthread_spin_unlock(&spin);

      successor = (pthread_t*)malloc(sizeof(successor));
      pthread_create(successor, NULL, pipe, &value);
    }
  }
  pthread_spin_lock(&spin);
      --thread_num;
  pthread_spin_unlock(&spin);
  __pseudo_pipe_barrier();

  if(successor){
    pthread_join(*successor, NULL);
  }
  printf("PIPE IS FINISHED WITH VALUE %d\n", own_value);

  return NULL;
}

int main(int argc, char** argv){
  int arg = atoi(argv[1]);
  int* value = malloc(sizeof(int));
  int own_value;
  pthread_t successor;
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init(&mutex, NULL);
  pthread_spin_init(&spin, 0);
  pthread_create(&successor, NULL, pipe, &value);
  *value = rand() % 100;
  for(int i = 0; i < arg; ++i){
    printf("ITERATION %d\n", i);
    __pipe_barrier();
    __pipe_barrier();
    if(i == (arg - 1)){
      value = NULL;
      pthread_spin_lock(&spin);
        --thread_num;
      pthread_spin_unlock(&spin);
      __pseudo_pipe_barrier();
    } else
      *value = rand() % 100;
  }
  pthread_join(successor, NULL);

}
