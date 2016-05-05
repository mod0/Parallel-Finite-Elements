OBJS = output.o grid.o main.o vector.o matrix.o domain.o elements.o assemble.o error.o fe.o
CC = gcc
CFLAGS = -Wall -g

all: $(OBJS)
	echo "Making all objects"

clean:
	rm -rf *.o *.out
