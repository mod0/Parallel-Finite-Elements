
OBJS = grid.o main.o fe.o domain.o elements.o

all: $(OBJS)
	echo "Making all objects"

clean:
	rm -rf *.o
