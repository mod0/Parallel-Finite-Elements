OBJS = output.o grid.o main.o vector.o matrix.o domain.o elements.o assemble.o fep.o solver.o
CC = gcc
CFLAGS = -Wall -g $(f)

all: clean $(OBJS)
	$(CC) $(OBJS) -o main -g -lm $(f)

clean:
	rm -rf *.o *.out
