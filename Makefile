
OBJS = grid.o main.o fe.o domain.o elements.o error.o
CC = gcc
CFLAGS = -Wall -Werror

all: $(OBJS)
	echo "Making all objects"

clean:
	rm -rf *.o
