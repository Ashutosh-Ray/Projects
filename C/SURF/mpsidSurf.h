
#ifndef _MPSID_SURF_H_
#define _MPSID_SURF_H_

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <stdio.h> 
#include <string.h>
 
#include <mpsid.h>

 
        /*-----------------------------------------*/
        /*---           Estructures             ---*/
        /*-----------------------------------------*/

typedef struct MPSID_SURF_EXTRAS {

char *ifile;  
char *ofile;

char *imgType;

int octaves;
float hess_thrs;

} MPSID_SURF_EXTRAS;
  
   
        /*-----------------------------------------*/
        /*---        def. funciones             ---*/
        /*-----------------------------------------*/

int mpsid_surf(void *data);
        
#endif

