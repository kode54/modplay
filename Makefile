CFLAGS = -c -fPIC

KISSFFT_OBJS = kissfft/kiss_fft.o kissfft/kiss_fftr.o

OBJS = resampler.o simple_convolver.o dbopl.o st3play.o ft2play.o $(KISSFFT_OBJS)


OPTS = -O3

all: libmodplay.a

libmodplay.a : $(OBJS)
	$(AR) rcs $@ $^

.c.o:
	$(CC) $(CFLAGS) $(OPTS) -o $@ $*.c

clean:
	rm -f $(OBJS) libmodplay.a > /dev/null

