OBJS = output.o grid.o vector.o matrix.o domain.o elements.o assemble.o fep.o solver.o
CC = gcc
CFLAGS = -Wall -g $(f)

all: elliptic parabolic

elliptic: clean $(OBJS) ellipticmain.o
	$(CC) $(OBJS) ellipticmain.o -o ellipticmain -g -lm $(f)

parabolic: clean $(OBJS) parabolicmain.o
	$(CC) $(OBJS) parabolicmain.o -o parabolicmain -g -lm $(f)


clean:
	rm -rf *.o *.out ellipticmain parabolicmain main
