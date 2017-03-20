
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdint.h>

#include <math.h>
#include <stdbool.h>

#include <sys/time.h>

#include <getopt.h>

#include "surf_opts.h"

    /*-----------------------------------------*/
    /*---            functions              ---*/
    /*-----------------------------------------*/


    /*-----------------------------------------*/
    /*---                main               ---*/
    /*-----------------------------------------*/

int main(int argc, char **argv)
{
    
  MPSID_PARSE_UOP parse;
  MPSID_SURF_EXTRAS extras;
  int l, n, k, Ncols, Nrows, Nbands;

  /* set extra parameters (input/output file) */
  parse.op.extras = ((void *) (&extras));
    
  /*---------------------------*/
  /*   read input parameters   */
  /*---------------------------*/

  parse_input(argc, argv, &parse);

  /*-------------------------*/
  /* choose CPU-Tick version */
  /*-------------------------*/

  mpsid_setMethod_TimeMeasurement(&(parse.ppars));

  /*----------------------------------------------------*/
  /*          get memory >>> (HOST & DEV) <<<           */
  /*----------------------------------------------------*/
   
  mpsid_GetMemUOP(&(parse.op));

  /*-----------------------------------*/
  /*      read/generate input data     */
  /*-----------------------------------*/    

  read_ppm2uchar(extras.ifile, ((unsigned char *) parse.op.aux[0]->buf),parse.op.In[0]->nData, ((char *) &extras.imgType));

 
  /* cast from char to float */

  Ncols = parse.op.In[0]->dim[0];
  Nrows = parse.op.In[0]->dim[1];
  Nbands = parse.op.In[0]->dim[2];
  
   
  for(k=0; k<Nrows; ++k) {
    for(n=0; n<Ncols; ++n)
      for(l=0; l<Nbands; ++l) ((float *) parse.op.In[0]->buf)[(k*Ncols+n)*Nbands+l] = ((float) (((unsigned char *) parse.op.aux[0]->buf)[(k*Ncols+n)*Nbands+l]));
      
  }
  
  /*--------------------------*/
  /*    call main routine     */
  /*--------------------------*/

  for(l=0; l<parse.ppars.loops; ++l) {
      
  //  mpsid_startCPUTick(&(parse.ppars));
    mpsid_surf(((void *) (&parse)));
  //  mpsid_endCPUTick(&(parse.ppars));

  }
  
   
    /*----------------------------*/
    /*     save output data       */
    /*----------------------------*/
  //save_data(((float *) parse.op.aux[2]->buf), parse.op.aux[2]->nData, "integral.raw"); 
  //save_data(((float *) parse.op.aux[6]->buf), parse.op.aux[3]->nData, "response_0.raw"); 	
  //save_data(((float *) parse.op.aux[8]->buf), parse.op.aux[4]->nData, "response_1.raw"); 	
  //save_data(((float *) parse.op.aux[9]->buf), parse.op.aux[5]->nData, "response_2.raw"); 	
  save_data(((float *) parse.op.Out[0]->buf), parse.op.Out[0]->nData, "intPoints.raw"); 	  
  save_data(((float *) parse.op.Out[1]->buf), parse.op.Out[1]->nData, "descriptor.raw"); 	
  /* -------------------- */
  /* Free  >>> Device <<< */
  /* -------------------- */
    
  mpsid_CleanupUOP(&(parse.op));
      
  return 0;   
 
}

