
OBJS = grid.o main.o vector.o matrix.o fe.o domain.o elements.o error.o
CC = gcc
CFLAGS = -Wall

all: $(OBJS)
	echo "Making all objects"

clean:
	rm -rf *.o
