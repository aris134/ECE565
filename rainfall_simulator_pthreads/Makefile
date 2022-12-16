
CC = gcc
OPT = -O2
LIBS = -pthread
EXECS = rainfall_seq rainfall_pt

all-redirect: $(EXECS)

rainfall_seq: rainfall_seq.o
	$(CC) -o rainfall_seq $(OPT) rainfall_seq.o

rainfall_pt: rainfall_pt.o
	$(CC) $(LIBS) $(OPT) -o rainfall_pt rainfall_pt.o

.c.c:
	$(CC) -c $(OPT) $<

clean:
	rm -rf rainfall_seq rainfall_pt *.o
