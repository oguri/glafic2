UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
	ifeq ($(UNAME_M),x86_64)
		CFLAGS2 = -lcurl
		LIBPATH = /usr/local/lib
		INCPATH = /usr/local/include
		PY_LIBS = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
	else
		CFLAGS2 = -lcurl
		LIBPATH = /opt/homebrew/lib
		INCPATH = /opt/homebrew/include
		PY_LIBS = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
	endif
else
	ifeq ($(UNAME_S),Linux)
		CFLAGS2 = 
		LIBPATH = /usr/local/lib
		INCPATH = /usr/local/include
		PY_LIBS = -lm $(LIBPATH)/libcfitsio.so $(LIBPATH)/libfftw3.so $(LIBPATH)/libgsl.so $(LIBPATH)/libgslcblas.so
	else
		CFLAGS2 = 
		LIBPATH =
		INCPATH =
		PY_LIBS = 
	endif
endif

CC	= gcc
CFLAGS	= -O2 -Wall -fPIC -I$(INCPATH)
#LIBS    = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
LIBS	= -lm $(LIBPATH)/libcfitsio.a $(LIBPATH)/libfftw3.a $(LIBPATH)/libgsl.a $(LIBPATH)/libgslcblas.a 

# for binary program
BIN	= glafic
OBJ_BIN	= glafic.o 
OBJS	= call.o ein_tab.o mass.o util.o fits.o \
	  init.o distance.o gsl_zbrent.o gsl_integration.o \
	  source.o extend.o point.o opt_extend.o opt_lens.o \
	  opt_point.o example.o mock.o calcein.o vary.o \
	  gnfw_tab.o commands.o mcmc.o amoeba_opt.o app_ell.o

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
