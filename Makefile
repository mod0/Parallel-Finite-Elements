OBJS = output.o grid.o main.o vector.o matrix.o domain.o elements.o assemble.o fep.o
CC = gcc
CFLAGS = -Wall -g

all: $(OBJS)
	echo "Making all objects"
	$(CC) grid.o main.o vector.o matrix.o domain.o elements.o -o main -g  -lm

clean:
	rm -rf *.o *.out
