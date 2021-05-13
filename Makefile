CC	= gcc
CFLAGS	= -O2 -Wall -fPIC
CFLAGS2 = -static
#LIBS    = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
LIBS	= -lm /usr/local/lib/libcfitsio.a /usr/local/lib/libfftw3.a /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a 

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CFLAGS2 = -lcurl
else
	CFLAGS2 = 
endif

# for binary program
BIN	= glafic
OBJ_BIN	= glafic.o 
OBJS	= call.o ein_tab.o mass.o util.o fits.o \
	  init.o distance.o gsl_zbrent.o gsl_integration.o \
	  source.o extend.o point.o opt_extend.o opt_lens.o \
	  opt_point.o example.o mock.o calcein.o vary.o \
	  gnfw_tab.o commands.o mcmc.o amoeba_opt.o 

# for C library
LIB	= libglafic.a
AR	= ar
AFLAGS	= r

# for python interface
PY	= glafic.so
OBJ_PY	= python.o
CFLAGS3	= -Wall -shared 
PY_INC  := $(shell python3-config --includes)
PY_LDS  := $(shell python3-config --ldflags --embed)
PY_LIBS = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
#PY_LIBS = -lm /usr/local/lib/libcfitsio.so /usr/local/lib/libfftw3.so /usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.so

default: bin

all: bin lib python

bin: $(OBJ_BIN) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ_BIN) $(OBJS) $(CFLAGS2) $(LIBS) 

lib: $(OBJS)
	$(AR) $(AFLAGS) $(LIB) $(OBJS) 

python: $(OBJ_PY) $(OBJS)
	$(CC) $(CFLAGS3) -o $(PY) $(OBJ_PY) $(OBJS) $(PY_LDS) $(PY_LIBS) 

python.o:python.c glafic.h 
	$(CC) $(CFLAGS) $(PY_INC) -c $< -o $@

%.o:%.c glafic.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(BIN) $(OBJ_BIN) $(OBJS) $(LIB) $(PY) $(OBJ_PY) *~ \#* core* 
