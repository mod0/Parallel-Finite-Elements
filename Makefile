
OBJS = grid.o main.o vector.o matrix.o domain.o element.o error.o fe.o
CC = gcc
CFLAGS = -Wall -g

all: $(OBJS)
	echo "Making all objects"

clean:
	rm -rf *.o *.out
