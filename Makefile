CC	= gcc
CFLAGS	= -O2 -Wall
#CFLAGS	= 
CFLAGS2 = -static
#LIBS    = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
LIBS	= -lm /usr/local/lib/libcfitsio.a /usr/local/lib/libfftw3.a /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -lcurl

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
PYT	= glafic.so
OBJ_PYT	= python.o
CFLAGS3	= -Wall -shared 
PINCDIR	= /usr/local/Cellar/python@3.8/3.8.5/Frameworks/Python.framework/Versions/3.8/include/python3.8/
PLIBDIR	= /usr/local/Cellar/python@3.8/3.8.5/Frameworks/Python.framework/Versions/3.8/lib
PLIBS	= -lpython3.8 -lm -lcfitsio -lfftw3 -lgsl -lgslcblas

ifeq ($(HOSTTYPE),intel-mac)
	CFLAGS2 = 
endif

default: bin

all: bin lib python

bin: $(OBJ_BIN) $(OBJS)
	$(CC) $(CFLAGS) $(CFLAGS2) -o $(BIN) $(OBJ_BIN) $(OBJS) $(LIBS) 

lib: $(OBJS)
	$(AR) $(AFLAGS) $(LIB) $(OBJS) 

python: $(OBJ_PYT) $(OBJS)
	$(CC) $(CFLAGS3) -L$(PLIBDIR) $(PLIBS) -o $(PYT) $(OBJ_PYT) $(OBJS)

python.o:python.c glafic.h 
	$(CC) $(CFLAGS) -I$(PINCDIR) -c $< -o $@

%.o:%.c glafic.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(BIN) $(OBJ_BIN) $(OBJS) $(LIB) $(PYT) $(OBJ_PYT) *~ \#* core* 
