CC	= gcc
CFLAGS	= -O2 -Wall
#CFLAGS	= 
LIBS    = -lm -lcfitsio -lfftw3 -lgsl -lgslcblas
#LIBS	= -lm /usr/local/lib/libcfitsio.a /usr/local/lib/libfftw3.a /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -lcurl
TARGET	= glafic
OBJS	= glafic.o ein_tab.o mass.o util.o fits.o \
	  init.o distance.o gsl_zbrent.o gsl_integration.o \
	  source.o extend.o point.o opt_extend.o opt_lens.o \
	  opt_point.o example.o mock.o calcein.o vary.o \
	  gnfw_tab.o commands.o mcmc.o amoeba_opt.o

CFLAGS2 = -static

ifeq ($(HOSTTYPE),intel-mac)
	CFLAGS2 = 
endif

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(OBJS) $(LIBS) 

%.o:%.c glafic.h 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(TARGET) $(OBJS) *~ \#* core* 
