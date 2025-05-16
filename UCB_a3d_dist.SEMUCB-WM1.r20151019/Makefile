CC	= gcc
CFLAGS	= -O2 -Wall
LDFLAGS =
LIBS	= -lm

O	= obj
S	= src

OBJ	= $O/ucb_A3d.o $O/main.o

BIN	= a3d_dist.x

$(BIN):	dirs $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) $(LIBS) -o $@

$O/%.o:	$S/%.c
	$(CC) -c $(CFLAGS) -Iinclude $< -o $@

.PHONY:	dirs
dirs:
	mkdir -p $O

.PHONY:	clean
clean:
	rm -rf $O

.PHONY:	distclean
distclean:
	rm -rf $O $(BIN)

