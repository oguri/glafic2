CC	= gcc
CFLAGS	= -O2 -Wall
#CFLAGS	= 
LIBS    = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas -lglafic -L..
#LIBS	= -lm /usr/local/lib/libcfitsio.a /usr/local/lib/libfftw3.a /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -lcurl
TARGET	= test
OBJS	= test.o

CFLAGS2 = -static

ifeq ($(HOSTTYPE),intel-mac)
	CFLAGS2 = 
endif

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(OBJS) $(LIBS) 

%.o:%.c 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(TARGET) $(OBJS) *~ \#* core* 
