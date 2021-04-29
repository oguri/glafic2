CC	= gcc
CFLAGS	= -O2 -Wall
#LIBS    = -L.. -lglafic -lcfitsio -lfftw3 -lgsl -lgslcblas -lm
LIBS    = -L.. -lglafic /usr/local/lib/libcfitsio.a /usr/local/lib/libfftw3.a /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -lm
TARGET	= test
OBJS	= test.o

CFLAGS2 = -static

ifeq ($(HOSTTYPE),intel-mac)
	CFLAGS2 = -lcurl
else
        CFLAGS2 = 
endif

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(OBJS) $(LIBS) 

%.o:%.c 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(TARGET) $(OBJS) *~ \#* core* 
