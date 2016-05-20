
CC=gcc-4.4
MEX = /usr/local/MATLAB/R2012a/bin/mex
# MEX=/MATLAB/bin/mex


# optimizaciones
OPT_1=-O3 -fomit-frame-pointer -malign-double -fstrict-aliasing 
OPT_2=-ffast-math -fPIC 

OPTIONS= $(OPT_1) $(OPT_2)

INCDIR=-I./ -I/usr/include/
#LIBDIR=-L/usr/lib64

#INCDIR=-I/usr/include 
LIBDIR=-L/usr/lib 


# cambia esto por el directorio adecuado
INCMATLAB=-I/usr/local/MATLAB/R2012a/extern/include

CFLAGS= $(OPTIONS) $(INCDIR) $(LIBDIR) 


# incluye otras librerias
LIB=-lm 

LIBMEX=$(LIBDIR) $(LIB)

###############################


TARGET= mex_ssd_simd
OBJ= $(TARGET).o
SRC= $(TARGET).c


$(TARGET): $(SRC)
	$(MEX) CC=$(CC) CFLAGS="$(OPTIONS)" $(INCDIR) $(LIBMEX) $(SRC)


clean:
	rm -rf *~
	rm -rf $(OBJ) 
	rm -rf $(TARGET)
