
/*>>> siempre el 1er include <<<*/
#include <mpsid.h>
#include <mpsidParseCodelets.c>

/*>>> otros includes ansiC de MPSID (orden es importante) <<<*/
#include <mpsidCpuTick.c>
#include <mpsidUtils.c>
#include <mpsidMemmgt.c>
#include <mpsidSurfFunctions.c>

/*>>> includes cuda (particular al programa transpose) <<<*/
#include "mpsidSurf.c"

enum {

  CPU_SERIAL = 1,
  LOOPS,
  INFILE,
  OUTFILE,
  OCTAVES,
  HESSIAN_THRESHOLD, 
  VERBOSE,
  CPUTICKS,
  SAVETICKS

};

struct option opts_list[] = {

  {"serial",1,0,CPU_SERIAL},
  {"loops",1,0,LOOPS},
  {"ii",1,0,INFILE},
  {"io",1,0,OUTFILE},
  {"octaves",1,0,OCTAVES},
  {"hthrs",1,0,HESSIAN_THRESHOLD},
  {"verbose",0,0,VERBOSE},
  {"cpuTime",0,0,CPUTICKS},
  {"saveTicks",0,0,SAVETICKS},
  {0,0,0,0}

};

        /*-----------------------------------------*/
        /*---         def. funtions             ---*/
        /*-----------------------------------------*/

int parse_input(int argc, char **argv, MPSID_PARSE_UOP *parse);

        /*-----------------------------------------*/
        /*---            routines               ---*/
        /*-----------------------------------------*/

int parse_input(int argc, char **argv, MPSID_PARSE_UOP *parse)
{

  int c, k;

  char ofile[] = "default";
  MPSID_SURF_EXTRAS *extras;

  extras = ((MPSID_SURF_EXTRAS *) parse->op.extras);


	/* Default values */
  extras->ifile = NULL;
  extras->ofile = NULL;
  extras->octaves = 5;
  extras->hess_thrs = 6.5; // matlab = 0.0001 image [0 1], 0.0001 * 255*255 
  
  extras->imgType  = ((char *) calloc(2, sizeof(char)));  

  parse->ppars.flagClkTicks = MPSID_CPUTICK;
  parse->ppars.loops = 1;
  parse->ppars.stages = 3;
  parse->ppars.serial = true;
  parse->ppars.mpsid_get_tick = NULL;
  parse->ppars.mpsid_get_base = NULL;
  parse->ppars.flagPrintLog = 1;      
  parse->ppars.flagPrintPart = 1;
  parse->ppars.cpu_tick.elapsed  = NULL;
  parse->ppars.cuda_tick.elapsed = NULL;


  while((c = getopt_long(argc, argv, "", opts_list, NULL)) != EOF)

  switch(c) {

    case CPU_SERIAL:
      mpsid_parse_set_cpuidSerial(&(parse->ppars.serial));
      break;

    case LOOPS:
      mpsid_parse_set_PositiveValue(&parse->ppars.loops, "Number of loops must be positive", 1);
      break;

    case CPUTICKS:
      parse->ppars.flagClkTicks = MPSID_CPUTICK;
      break;

    case SAVETICKS:
      parse->ppars.flagPrintLog = 0;
      parse->ppars.flagClkTicks = MPSID_ALLTICK;
      break;

    case VERBOSE:
      /* do something? */
      break;

    case OCTAVES: 
      mpsid_parse_set_PositiveValue((int *) &(extras->octaves), "Number of octaves must be greater than 1 and less or equal than 5", extras->octaves);
      break;

    case HESSIAN_THRESHOLD: 
      mpsid_parse_set_PositiveFloatValue((float *) &(extras->hess_thrs), "Threshold must be positive", extras->hess_thrs);
      break;

    case INFILE:
      mpsid_parse_set_Filename(&(extras->ifile));
      break;

    case OUTFILE:
      mpsid_parse_set_Filename(&(extras->ofile));
      break;

    default:
      printf("\nNot a valid command\n");
      return(-1);      

  }

  if(extras->ifile == NULL) {
    err_info("no input file... exiting\n");
    exit(1);
  }
  
  int nresponse;
  
  //  por escala se halla un det_hessian
  nresponse = 4 + (extras->octaves - 1)*2; 	
	
    /***************************************************************
	  IMPORTANT!!!!! 
	  * 
	* surf: 1 inputs (image), 
	* 2 output ( points, descriptors), 
	* 03 aux vectors 	|   SO FAR 
	* 					| 2 for reading /writting files 
	* 					| 1 for integral image
	* 					| el resto es para las escalas ( det_hessian)
	***************************************************************/
  mpsid_InitUOp(&(parse->op), 1, 2, 3 + nresponse); 
  
  mpsid_InitVector(parse->op.In[0], 3, sizeof(float));   /* INPUT */ 

  for(k=0; k<parse->op.nOutputs; ++k) mpsid_InitVector(parse->op.Out[k], 2, sizeof(float));

  mpsid_InitVector(parse->op.aux[0], 3, sizeof(char)); /* aux char for reading file */
  mpsid_InitVector(parse->op.aux[1], 3, sizeof(char)); /* aux char for writing file */
  mpsid_InitVector(parse->op.aux[2], 2, sizeof(float)); /* aux float for integral image */
   
  for(k=3; k<parse->op.nAux; ++k) mpsid_InitVector(parse->op.aux[k], 2, sizeof(float));

   
  for(k=0; k<parse->op.nInputs; ++k) {
    parse->op.In[k]->flagCUDAMemType = MPSID_CUDA_MEM_NONE;
    parse->op.In[k]->flagHostGetMem = MPSID_MEM_POSIX_ALGN;
    parse->op.In[k]->major = MPSID_ROWMAJOR;
  }
  
  for(k=0; k<parse->op.nOutputs; ++k) {
	parse->op.Out[k]->flagCUDAMemType = MPSID_CUDA_MEM_NONE;
    parse->op.Out[k]->flagHostGetMem = MPSID_MEM_POSIX_ALGN;
    parse->op.Out[k]->major = MPSID_ROWMAJOR;
  }
		
  for(k=0; k<parse->op.nAux; ++k) {
    parse->op.aux[k]->flagCUDAMemType = MPSID_CUDA_MEM_NONE;
    parse->op.aux[k]->flagHostGetMem = MPSID_MEM_POSIX_ALGN;
    parse->op.aux[k]->major = MPSID_ROWMAJOR;
  }
  	
  read_headppm(extras->ifile, &(parse->op.In[0]->dim[0]), &(parse->op.In[0]->dim[1]), &(parse->op.In[0]->dim[2]), ((char *) &extras->imgType));


	/**  Setting Output dimensions  */
  parse->op.Out[0]->dim[0] = parse->op.In[0]->dim[1]*parse->op.In[0]->dim[0];     
  parse->op.Out[0]->dim[1] = 5;    // Rows (x,y,scale,orientation,laplacian)
  
  parse->op.Out[1]->dim[0] = parse->op.In[0]->dim[1]*parse->op.In[0]->dim[0];     
  parse->op.Out[1]->dim[1] = 64;    // Rows

	/**  Setting Aux dimensions  */
  for(k=0; k<2; ++k) {
    parse->op.aux[k]->dim[0] = parse->op.In[0]->dim[0];     
    parse->op.aux[k]->dim[1] = parse->op.In[0]->dim[1]; 
	parse->op.aux[k]->dim[2] = parse->op.In[0]->dim[2];
  }
  

	/* todos nuestras auxiliarias son matrices, hasta el momento*/
  for(k=2; k<parse->op.nAux; ++k) {
    parse->op.aux[k]->dim[0] = parse->op.In[0]->dim[0];     
    parse->op.aux[k]->dim[1] = parse->op.In[0]->dim[1]; 
  }
     

  /* check output file */
  if(extras->ofile == NULL) {
     extras->ofile = calloc( strlen(ofile)+1, sizeof(char));
     strcpy(extras->ofile, ofile);
  }
    
  mpsid_SetAllVectorsnBytes(parse->op);

  mpsid_initTicks(&(parse->ppars));

  return(1);

}

