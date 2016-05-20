#include "mex.h"
#include "matrix.h"

float ssd(float *q, int lq, float *p, int lp, int N);

void mexFunction(int nlhs, mxArray *plhs[ ],
			int nrhs, const mxArray *prhs[ ]) 

{
	// Q - Neighboor    |   P - Node
	//  col - major
	int i, j, pm,qm, N;
	float *e, *p, *q;
	double *Lp , *Lq; 

	/* Find the dimensions of the data */
    pm = mxGetN(prhs[0]);  // 3xN*Lp
	qm = mxGetN(prhs[1]);  // 3xN*Lq
	
    /* Retrieve the input data */
    p = (float*)mxGetPr(prhs[0]);
    q = (float*)mxGetPr(prhs[1]);
    Lp = mxGetPr(prhs[2]); // Number of p-labels  
    Lq = mxGetPr(prhs[3]); // Number of q-labels  
    N = qm / (int)Lq[0]; // Number of pixels
    
    /* Create an mxArray for the output data */
    plhs[0] = mxCreateNumericMatrix(Lp[0], Lq[0], mxSINGLE_CLASS, 0); 	
    e = (float*)mxGetPr(plhs[0]);
   /* 
    printf("Lp = %i  | Lq = %i \n",(int)Lp[0],(int)Lq[0]);
    printf("Pm = %i  | Qm = %i \n",pm,qm);
    printf("PatchDim = %i \n",N);
    */
    for (j = 0; j < (int)Lq[0] ; j++){
		for (i =0 ; i < (int)Lp[0] ; i++){
			e[j*(int)Lp[0] + i] =  ssd(q,j,p,i,N);
		}
	}    
}
  
  
 float ssd(float *q, int lq, float *p, int lp, int N)
 {
	 // lq -> offset = lq*N*3
	 // lp -> offset = lp*N*3
	 float e = 0;
	 int i,j, p_offset, q_offset; 
	 
	 p_offset = lp*N*3;
	 q_offset = lq*N*3;
	 
	 for (j=0; j<N; ++j)
	 {
	   for (i=0; i<3; ++i)
		{
			e += (q[q_offset + i+j*3] - p[p_offset + i+j*3])*(q[q_offset + i+j*3] - p[p_offset + i+j*3]);
		}
	 }
	
	 return e;
 }
