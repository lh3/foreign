CC=			gcc
CFLAGS=		-g -Wall -O2
OBJS=		regexp9.o
PROG=		test
INCLUDES=	-I.
LIBS=		
SUBDIRS=	.

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

test:regexp9.o test.o
		$(CC) test.o regexp9.o -o $@

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
