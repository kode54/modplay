CFLAGS = -c -fPIC

OBJS = fir_resample.o resampler.o dbopl.o st3play.o ft2play.o


OPTS = -O3

all: libmodplay.a

libmodplay.a : $(OBJS)
	$(AR) rcs $@ $^

.c.o:
	$(CC) $(CFLAGS) $(OPTS) -o $@ $*.c

.cpp.o:
	$(CXX) $(CFLAGS) $(OPTS) -o $@ $*.cpp

clean:
	rm -f $(OBJS) libmodplay.a > /dev/null

