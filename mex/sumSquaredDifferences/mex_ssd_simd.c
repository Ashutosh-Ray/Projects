#include "mex.h"
#include "matrix.h"

float ssd(float *q, int lq, float *p, int lp, int N);
float ssd_SIMD(float *q, int lq, float *p, int lp, int N);
void convert_double2float( double *input_double, float *output_float,int Ntot);
void convert_float2double( float *input_float, double *output_double,int Ntot);

void mexFunction(int nlhs, mxArray *plhs[ ],
			int nrhs, const mxArray *prhs[ ]) 

{
	// Q - Neighboor    |   P - Node
	//  col - major
	int i, j, pm,qm, N,channels;
	float *e, *p, *q;
	double *ed, *pd ,*qd;
	double *Lp , *Lq , *cptr; 

	/* Find the dimensions of the data */
    pm = mxGetN(prhs[0]);  // 3xN*Lp
	qm = mxGetN(prhs[1]);  // 3xN*Lq
	
    /* Retrieve the input data */
   // p = (float*)mxGetPr(prhs[0]);
   // q = (float*)mxGetPr(prhs[1]);
    
    pd = mxGetPr(prhs[0]);
    qd = mxGetPr(prhs[1]);
    Lp = mxGetPr(prhs[2]); // Number of p-labels  
    Lq = mxGetPr(prhs[3]); // Number of q-labels  
    cptr = mxGetPr(prhs[4]); // Numer of channels
    channels = (int)cptr[0];
    N = qm / (int)Lq[0]; // Number of pixels
    
    /* Create an mxArray for the output data */
    // plhs[0] = mxCreateNumericMatrix(Lp[0], Lq[0], mxSINGLE_CLASS, 0); 
    plhs[0]= mxCreateDoubleMatrix(Lp[0], Lq[0], mxREAL);	
    ed = mxGetPr(plhs[0]);
    
    p = (float*) mxMalloc(sizeof(float)*channels*N*Lp[0]);
	q = (float*) mxMalloc(sizeof(float)*channels*N*Lq[0]);
	e = (float*) mxMalloc(sizeof(float)*Lp[0]*Lq[0]);
	
	convert_double2float(pd,p,  channels*N*Lp[0]);
	convert_double2float(qd,q,  channels*N*Lq[0]);
    //e = (float*)mxGetPr(plhs[0]);
   
	//printf("N = %i  qm = %i Lq =%i \n",N,qm,(int)Lq[0]);
    for (j = 0; j < (int)Lq[0] ; j++){
		for (i =0 ; i < (int)Lp[0] ; i++){
			e[j*(int)Lp[0] + i] =  ssd_SIMD(q,j,p,i,N*channels);
		}
	}    
	
	convert_float2double(e, ed, Lp[0]*Lq[0]);
	
	mxFree(p);
	mxFree(q);
	mxFree(e);
}
  
  
float ssd_SIMD(float *q, int lq, float *p, int lp, int N)
{
  int k, l,p_offset, q_offset;
  float *Ps, *Qs, *e, E;
 
	p_offset = lp*N;
	q_offset = lq*N;
	
	l = N%4; 
  
   Ps = p + p_offset;
   Qs = q + q_offset;
   /*
   printf ("\n\n");
   for (k=0 ; k<N; k++){
		printf("Qs = %f  | Ps = %f \n",Qs[k],Ps[k]);
	}
   
    printf("dir p %p    | dir q %p\n",p,q);
    printf("p_off %i    |  q_off %i\n",p_offset,q_offset);
    printf("dir ps %p   | dir qs %p\n",Ps,Qs);
    */
	posix_memalign( (void**)&e, 16, sizeof(float));
	
    e[0]=0;
       
    __asm__ __volatile__(
          "movss  0x0(%2), %%xmm2       \n\t" 
        :
        : "r" (Qs), "r" (Ps), "r" (e)
        : "memory");
    
    for (k=0 ; k<((N-l)>>2) ; k++){
    __asm__ __volatile__(
          "movups 0x0(%0), %%xmm0       \n\t"            /* xmm0: | p3 | p2 | p1 | p0 | */
          "movups 0x0(%1), %%xmm1  		\n\t"	 		 /* xmm1: | q3 | q2 | q1 | q0 | */
          "subps  %%xmm1 , %%xmm0       \n\t"			 /* xmm0: | q3-p3 | q2-p2 | q1-p1 | q0-p0 | */
          "mulps  %%xmm0 , %%xmm0 		\n\t"			 /* xmm0: | pq3^2 | qp2^2 | qp1^2 | qp0^2 | */
          "haddps %%xmm0 , %%xmm0 		\n\t"
          "haddps %%xmm0 , %%xmm0 		\n\t"
          "addss  %%xmm0 , %%xmm2       \n\t" 
        :
        : "r" (Qs), "r" (Ps), "r" (e)
        : "memory");
        Ps += 4;
        Qs += 4; 
     }
  
	__asm__ __volatile__(
          "movss  %%xmm2 , 0x0(%2)       \n\t" 
        :
        : "r" (Qs), "r" (Ps), "r" (e)
        : "memory");

	for ( k=0; k<l; k++){
		e[0] += (Ps[k]-Qs[k])*(Ps[k]-Qs[k]);
		}
		
	E = e[0];	
	
	free(e);
    return(E);
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
