CFLAGS = -c -fPIC

OBJS = resampler.o dbopl.o st3play.o ft2play.o


OPTS = -O3

all: libmodplay.a

libmodplay.a : $(OBJS)
	$(AR) rcs $@ $^

.c.o:
	$(CC) $(CFLAGS) $(OPTS) -o $@ $*.c

clean:
	rm -f $(OBJS) libmodplay.a > /dev/null

