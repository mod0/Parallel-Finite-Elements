OBJS = output.o grid.o main.o vector.o matrix.o domain.o elements.o assemble.o fep.o
CC = gcc
CFLAGS = -Wall -g

all: $(OBJS)
	$(CC) $(OBJS) -o main -g  -lm

clean:
	rm -rf *.o *.out
