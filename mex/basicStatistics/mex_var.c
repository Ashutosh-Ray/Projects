#include "mex.h"
#include "matrix.h"

void 	convert_double2float( double *input_double, float *output_float,int Ntot);
void 	convert_float2double( float *input_float, double *output_double,int Ntot);
float 	mean_SIMD(float *data, int N);
float 	var_SIMD(float *data, int N, float mean);


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
	// Q - Neighboor    |   P - Node
	//  col - major
	int   	k, N;
	double 	mean = 0;
	float 	*varf, *dataf;
	double 	*vard, *datad;

	/* Find the dimensions of the data */
    N 			= mxGetM(prhs[0]);  // Nx1 - Number of rows
	
    /* Retrieve the input data */
    datad 		=  mxGetPr(prhs[0]);
    dataf 		= (float*) mxMalloc(sizeof(float)*N);
    convert_double2float(datad, dataf,N);
    
    /* Create an mxArray for the output data */
    plhs[0] 	= mxCreateDoubleMatrix(1, 1, mxREAL); 	
    vard 		= mxGetPr(plhs[0]);
    varf    	= (float*) mxMalloc(sizeof(float));
    
	mean 		= mean_SIMD(dataf,N);	
	varf[0]		= var_SIMD(dataf,N,mean);
	
	vard[0] 	= (double) varf[0];
	
	/*
	mean 		= (double)mean;
	vard[0]  =0;
	for (k=0;k<N;k++){
		vard[0] += (datad[k] - mean)*(datad[k] - mean);
	}
	vard[0] /= (N-1);
	*/
	mxFree(dataf);
	mxFree(varf);
}
  

void convert_double2float( double *input_double, float *output_float,int Ntot)
{
    int i;
    for (i = 0; i < Ntot; i++)
    {
                output_float[i] = (float) input_double[i];
    }
}

void convert_float2double( float *input_float, double *output_double,int Ntot)
{
    int i;
    for (i = 0; i < Ntot; i++)
    {
                output_double[i] = (double) input_float[i];
    }
}

float mean_SIMD(float *data, int N)
{
	int k,l;
	float *mean, M;

	l = N%4; 
  
	posix_memalign( (void**)&mean, 16, sizeof(float));
	
    mean[0]=0;
       
    __asm__ __volatile__(
          "movss  0x0(%1), %%xmm2       \n\t" 			
        :
        : "r" (data),"r" (mean)
        : "memory");
    
    for (k=0 ; k<((N-l)>>2) ; k++){
    __asm__ __volatile__(
          "movups 0x0(%0), %%xmm0       \n\t"            /* xmm0: | p3 | p2 | p1 | p0 | */
          "haddps %%xmm0 , %%xmm0 		\n\t"			 
          "haddps %%xmm0 , %%xmm0 		\n\t"
          "addss  %%xmm0 , %%xmm2       \n\t" 
        :
        : "r" (data), "r" (mean)
        : "memory");
        data += 4;
     }
  
	__asm__ __volatile__(
          "movss  %%xmm2 , 0x0(%1)       \n\t" 
        :
        : "r" (data), "r" (mean)
        : "memory");

	for ( k=0; k<l; k++){
		mean[0] += data[k];
		}
		
	M = mean[0]/N;	
	
	free(mean);
    return(M);
}

float var_SIMD(float *data, int N, float mean)
{
	int k,l;
	float *var, *meanptr, V;

	l = N%4; 
  
	posix_memalign( (void**)&var, 16, sizeof(float));
	posix_memalign( (void**)&meanptr, 16, sizeof(float));
	
    var[0]		= 0;
    meanptr[0] 	= mean;
       
    __asm__ __volatile__(
          "movss  0x0(%1), %%xmm2       \n\t" 	       	/* xmm2: | 0 | 0 | 0 | v | */
          "movss  0x0(%2), %%xmm1 		\n\t"			/* xmm1: | 0 | 0 | 0 | m | */
          "shufps $0x00, %%xmm1, %%xmm1  \n\t" 		    /* xmm3: | m | m | m | m | */
        :
        : "r" (data),"r" (var), "r" (meanptr)
        : "memory");
    
    for (k=0 ; k<((N-l)>>2) ; k++){
    __asm__ __volatile__(
          "movups 0x0(%0), %%xmm0       \n\t"            /* xmm0: | d3 | d2 | d1 | d0 | */
          "subps  %%xmm1 , %%xmm0 		\n\t"			 /* xmm0: | m-d3 | m-d2 | m-d1 | m-d0 | */
          "mulps  %%xmm0 , %%xmm0 		\n\t"			 /* xmm0: | md3^2 | md2^2 | md1^2 | md0^2 | */
          "haddps %%xmm0 , %%xmm0 		\n\t"			 
          "haddps %%xmm0 , %%xmm0 		\n\t"
          "addss  %%xmm0 , %%xmm2       \n\t" 
        :
        : "r" (data), "r" (var), "r" (meanptr)
        : "memory");
        data += 4;
     }
  
	__asm__ __volatile__(
          "movss  %%xmm2 , 0x0(%1)       \n\t" 
        :
        : "r" (data), "r" (var), "r" (meanptr)
        : "memory");

	for ( k=0; k<l; k++){
		var[0] += (data[k]-meanptr[0])*(data[k]-meanptr[0]);
		}
		
	V = var[0]/(N-1);	
	
	free(var);
    return(V);
}
