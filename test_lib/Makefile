UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
        ifeq ($(UNAME_M),x86_64)
                CFLAGS2 = -lcurl
                LIBPATH = /usr/local/lib
                INCPATH = /usr/local/include
        else
                CFLAGS2 = -lcurl
                LIBPATH = /opt/homebrew/lib
                INCPATH = /opt/homebrew/include
        endif
else
        ifeq ($(UNAME_S),Linux)
                CFLAGS2 = 
                LIBPATH = /usr/local/lib
                INCPATH = /usr/local/include
        else
                CFLAGS2 = 
                LIBPATH =
                INCPATH =
        endif
endif

CC	= gcc
CFLAGS	= -O2 -Wall -I$(INCPATH)
#LIBS    = -L.. -lglafic -lcfitsio -lfftw3 -lgsl -lgslcblas -lm
LIBS    = -L.. -lglafic $(LIBPATH)/libcfitsio.a $(LIBPATH)/libfftw3.a $(LIBPATH)/libgsl.a $(LIBPATH)/libgslcblas.a -lm

all: test mock_siex

test: test.o
	$(CC) $(CFLAGS) -o $@ test.o $(CFLAGS2) $(LIBS) 

mock_siex: mock_siex.o
	$(CC) $(CFLAGS) -o $@ mock_siex.o $(CFLAGS2) $(LIBS) 

%.o:%.c 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f *.o *~ \#* core* 
