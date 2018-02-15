CURR_DIR=$(notdir $(basename $(shell pwd)))
PRJ=$(CURR_DIR)
SRC=$(wildcard *.c)
OBJ=$(patsubst %.c,%.o,$(SRC))

CC=clang
INCLUDES=-I../../include
ifndef DEBUG
CFLAGS=-lm -O3 -mllvm -force-vector-interleave=2
LIB=
else
CFLAGS=-O0 -g3 -std=c99 -fopt-info-vec -ftree-vectorize -march=native
LIB=
endif

all: $(PRJ)

$(PRJ): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJ) -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@ $(LIB)

clean:
	-rm -f $(OBJ) $(PRJ)
