CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=
OBJS=		msGCHOT.o streecGCHOT.o rand1.o
PROG=		msHOT-lite
LIBS=		-lm

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

msHOT-lite:$(OBJS)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)

clean:
		rm -f *.o a.out *~ *.a $(PROG)
