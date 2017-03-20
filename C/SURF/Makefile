CC=gcc

# optimizaciones
OPT_1= -Wall -O2 -O3 

# Uncomment for INTEL
OPT_2= -fomit-frame-pointer -fno-schedule-insns -malign-double 
OPT_3= -fstrict-aliasing -mpreferred-stack-boundary=4 -fPIC 
OPT_4= -ffast-math -ftree-vectorize -msse2 -D_MPSIDINTEL_

# Uncomment for CORTEXA8
#OPT_2= -mcpu=cortex-a8 -mfpu=neon -mfloat-abi=hard -D_CORTEXA8_
#OPT_3= -fsingle-precision-constant -Wno-unused-result
#OPT_4= -ffast-math -ftree-vectorize 

OPTIONS= $(OPT_1) $(OPT_2) $(OPT_3) $(OPT_4)


INCDIR=-I./  -I./../include  -I/usr/include
LIBDIR=-L/usr/lib


CFLAGS= $(OPTIONS) $(INCDIR) $(LIBDIR) 

LIB=-lm 

#-------------------------------------------------------
#-------------------------------------------------------

SRC= main.c
OBJ= main.o
TARGET= main

$(TARGET): $(OBJ) 
	$(CC) $(CFLAGS) -o $(TARGET)  $(OBJ) $(LIB)

$(OBJ): $(SRC) surf_opts.h 
	$(CC) $(CFLAGS) -c $(SRC) $(LIB)

clean:
	rm -rf *~
	rm -rf $(OBJ)
	rm -rf $(TARGET)
	rm -rf *.raw

