
#include <mpsidSurf.h>

int mpsid_surf(void *data)
{
  int Ncols, Nrows, octaves, nscales;
  float hess_thrs;
 
  float *iimg; 
  float *det_hess; 
  
  float *IntPoints;
  float *Descriptors;
  
  int k,j;
	
  MPSID_PARSE_UOP *parse = ((MPSID_PARSE_UOP *) data);
  MPSID_SURF_EXTRAS *extras = ((MPSID_SURF_EXTRAS *) (parse->op.extras));
  
  octaves 	= ((int)	extras->octaves);
  hess_thrs = ((float)	extras->hess_thrs);
  
  Nrows		= ((int)  parse->op.In[0]->dim[1]);
  Ncols		= ((int)  parse->op.In[0]->dim[0]);
      
  iimg 		= ((float *)  parse->op.aux[2]->buf);
  
  IntPoints		= ((float *)  parse->op.Out[0]->buf);
  Descriptors 	= ((float *)  parse->op.Out[1]->buf);
    
  mpsid_IntegralImage(iimg, ((float *) parse->op.In[0]->buf) , Nrows, Ncols);
  
  nscales = 4 + (octaves -1)*2;
  
  /*********************************************************************/
  /**                      SETTING SCALE INFO                         **/
  /*********************************************************************/
  //La memoria para los det_hessian ya ha sido reservada , se encuentra en aux[3] hasta el final
  MPSID_SCALE_INFO scale[nscales]; 
 

  int s; 
  
  s = 2;
    
  int *temp_nscale, *temp_filter_size; 
  int step;
   
  temp_nscale = (int*)calloc( nscales, sizeof(int));  
  temp_filter_size = (int*)calloc( nscales, sizeof(int));  
  
  for (k=0;k<4;k++) temp_nscale[k] = 0;
  for (k=4;k<nscales;k++) temp_nscale[k] = k - 2;	
  for (k=0;k<nscales;k++) temp_nscale[k] = (int)exp2(temp_nscale[k]/2); 
  
  memcpy(temp_filter_size, temp_nscale,nscales*sizeof(int));
  
  int foo = 0;
  for (k=0;k<nscales;k++) {
	  temp_filter_size[k] = foo + temp_filter_size[k];
	  foo = temp_filter_size[k]; 
  }
    
  for (k=0;k<nscales;k++) temp_filter_size[k] = 3 + 6*temp_filter_size[k];	
 
	
   for (k=0; k<nscales; k++ ){
	  step = temp_nscale[k];
	  det_hess 	= ((float *)  parse->op.aux[3+k]->buf);
	  mpsid_SetScale(k,&scale[k],det_hess, s*step, temp_filter_size[k]);
   }

	free(temp_nscale);
	free(temp_filter_size);
	
  /*********************************************************************/
  /**                       COMPUTE RESPONSES                         **/
  /*********************************************************************/
  
	for (k =0; k<nscales; k++) 	mpsid_getResponse (&scale[k], iimg, Nrows, Ncols); 

	
  /*********************************************************************/
  /**                      GET INTEREST POINTS                        **/
  /*********************************************************************/
	
	int *filter_map;
	int t,b,m;
	
	int nIntPoints;

	filter_map = (int*)calloc( octaves*4, sizeof(int)); // For 5 octaves
	
	filter_map[0] = 0;
	filter_map[1] = 1;
	filter_map[2] = 2;
	filter_map[3] = 3;
	
	
	for(k=1; k<octaves; k++) {
		filter_map[0+k*4] = k+(k-1);
		filter_map[1+k*4] = k+(k-1)+2;
		filter_map[2+k*4] = k+(k-1)+3;
		filter_map[3+k*4] = k+(k-1)+4;	
	}

	nIntPoints = 0;
	for (k=0; k<octaves; k++) {
		for (j=0; j<2;j++) {
			b = filter_map[ k*4 + j];
			m = filter_map[ k*4 + j+1];
			t = filter_map[ k*4 + j+2];
			mpsid_findExtremum (&nIntPoints, IntPoints ,&scale[t], &scale[m], &scale[b], hess_thrs, Nrows, Ncols) ;
		}
	}
	
	IntPoints[0] = nIntPoints;
	
	
   /*********************************************************************/
   /**                   DESCRIBE INTEREST POINTS                      **/
   /*********************************************************************/
	
	mpsid_describePoints(  IntPoints, Descriptors, iimg, Nrows, Ncols);
	
	free(filter_map);
	

  return 0;

}
