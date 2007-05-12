#include "copyright.h"
/*==============================================================================
 * FILE: output.c
 *
 * PURPOSE: Controls output of data.  Output is divided into three types:
 *   1. dump_*(): ALL variables are written in * format over WHOLE grid
 *   2. output_*(): ONE variable is written in * format with various options
 *   3. restarts: special form of a dump, includes extra data
 *   The number and types of outputs are all controlled by <ouputN> blocks in
 *   the input files, parsed by the functions in par.c.  
 *
 * TOTAL NUMBER of outputs is controlled by 'maxout' in <job> block in input
 *   file.  Only the first 'maxout' <outputN> blocks are processed, where
 *   N < maxout.  If N > maxout, that <outputN> block is ignored.
 *
 * OPTIONS available in an <outputN> block are:
 *   out         = all,d,M1,M2,M3,E,B1c,B2c,B3c,ME,V1,V2,V3,P,S,cs2
 *   out_fmt     = bin,dx,hst,tab,rst,vtk,fits,pdf,pgm,ppm
 *   dat_fmt     = format string used to write tabular output (e.g. 12.5e)
 *   dt          = problem time between outputs
 *   id          = any string
 *   dmin/dmax   = max/min applied to all outputs
 *   ix1,ix2,ix3 = range of indices to be dumped (see parse_slice())
 *   palette     = rainbow,jh_colors,idl1,idl2,step8,step32,heat
 *   usr_expr_flag = 1 for user-defined expression (defined in problem.c)
 *   
 * EXAMPLE of an <outputN> block for a VTK dump:
 *   <ouput1>
 *   out_fmt = vtk
 *   out_dt  = 0.1
 *
 * EXAMPLE of an <outputN> block for a ppm image of a XY-slice in ppm format:
 *   <output5>
 *   out_fmt = ppm
 *   dt      = 100.0
 *   out     = d
 *   id      = d
 *   ix3     = 64
 *   dmin    = 0.25
 *   dmax    = 2.9
 *   palette = rainbow
 *
 * EXAMPLE of an <outputN> block for restarts:
 *   <ouput3>
 *   out_fmt = rst
 *   out_dt  = 1.0
 *
 * CONTROL of output proceeds as follows:
 *  -init_output(): called by main(), parses the first maxout output blocks.
 *     The info in each block is stored in an element of a global array of
 *     "Output_s" structures, including a pointer to the appropriate output
 *     function.
 *  -data_output(): called in main loop, compares integration time with time 
 *     for output for each element in Output array, and calls output functions.
 *    
 *   To add permanently a new type of dump X, write a new function dump_X, then
 *   modify init_output() to set the output function pointer when out_fmt=X
 *   in the input file (see below for examples of bin, hst, rst, etc.)
 *    
 *   To add permanently a new type of output X, write a new function output_X, 
 *   modify init_output() to set the output function pointer when out_fmt=X
 *   in the input file (see below for examples of pgm, ppm, etc.)
 *
 *   To add a problem-specific user-defined output function in the problem
 *   definition file...
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   init_output() -
 *   data_output() -
 *   add_output()
 *   add_rst_out()
 *   data_output_destruct()
 *   data_output_enroll()
 *   subset1,2,3()   -
 *
 * VARIABLE TYPE AND STRUCTURE DEFINITIONS: none
 *============================================================================*/

#ifndef __POWERPC__
#define HAVE_DLFCN
#endif

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_DLFCN
#include <dlfcn.h>
#endif

#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "palette.h"
#include "prototypes.h"

#define MAXOUT_DEFAULT     10

static int out_count = 0;    /* Number of Outputs initialized in the OutArray */
static size_t out_size  = 0;       /* Size of the array OutArray[] */
static Output *OutArray = NULL;    /* Array of Output modes */
static Output rst_out;             /* Restart Output */
static int rst_flag = 0;           /* (0,1) -> Restart Outputs are (off,on) */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   expr_*
 *   get_expr
 *   free_output
 *   parse_slice
 *   getRGB
 *============================================================================*/

Real expr_d  (const Grid *pG, const int i, const int j, const int k);
Real expr_M1 (const Grid *pG, const int i, const int j, const int k);
Real expr_M2 (const Grid *pG, const int i, const int j, const int k);
Real expr_M3 (const Grid *pG, const int i, const int j, const int k);
Real expr_E  (const Grid *pG, const int i, const int j, const int k);
Real expr_B1c(const Grid *pG, const int i, const int j, const int k);
Real expr_B2c(const Grid *pG, const int i, const int j, const int k);
Real expr_B3c(const Grid *pG, const int i, const int j, const int k);
Real expr_ME (const Grid *pG, const int i, const int j, const int k);
Real expr_V1 (const Grid *pG, const int i, const int j, const int k);
Real expr_V2 (const Grid *pG, const int i, const int j, const int k);
Real expr_V3 (const Grid *pG, const int i, const int j, const int k);
Real expr_P  (const Grid *pG, const int i, const int j, const int k);
Real expr_cs2(const Grid *pG, const int i, const int j, const int k);
Real expr_S  (const Grid *pG, const int i, const int j, const int k);
static Gasfun_t getexpr(const int n, const char *expr);
static void free_output(Output *pout);
static void parse_slice(char *block, char *axname, int nx, Real x, Real dx,
			int *l, int *u, int *new_nx, Real *new_x, Real *new_dx);
float *getRGB(char *name);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* init_output:  */

void init_output(Grid *pGrid)
{
  int i,j,outn,nout=0,maxout;
  char block[80], *fmt, defid[10];
  Output new_out;
  int usr_expr_flag;

  maxout = par_geti_def("job","maxout",MAXOUT_DEFAULT);

/*--------------------- begin loop over maxout output numbers ----------------*/
  for (outn=1; outn<=maxout; outn++) {       /* check all outN up to 'maxout' */

    sprintf(block,"output%d",outn);
/* An output format is required, if not, we move on. */
    if(par_exist(block,"out_fmt") == 0){
      ath_perr(-1,"[init_output]: %s/out_fmt does not exist\n",block);
      continue;
    }

/* Zero (NULL) all members of the Output structure */
    memset(&new_out,0,sizeof(Output));

/* The next output time and number */
    new_out.t   = par_getd_def(block,"time",pGrid->time);
    new_out.num = par_geti_def(block,"num",0);

    new_out.dt  = par_getd(block,"dt");
    new_out.n   = outn;

/* set id in output filename to input string if present, otherwise use "outN"
 * as default, where N is output number */
    sprintf(defid,"out%d",outn);
    new_out.id = par_gets_def(block,"id",defid);

    fmt = new_out.out_fmt = par_gets(block,"out_fmt");

/* out:     controls what variable can be output (all, or any of expr_*)
 * out_fmt: controls format of output (single variable) or dump (all variables)
 *
 * if "out" doesn't exist, we assume 'all' variables are meant to be dumped */

    new_out.out = par_gets_def(block,"out","all");

/* First handle data dumps (ouput of ALL variables) */

    if(strcmp(new_out.out,"all") == 0){
/* check for valid data dump: dump format = {bin, dx, hst, tab, rst, vtk} */
      if (strcmp(fmt,"bin")==0){
	new_out.fun = dump_binary;
	goto add_it;
      }
      else if (strcmp(fmt,"dx")==0){
	new_out.fun = dump_dx;
	if (par_exist(block,"dat_fmt"))
	  new_out.dat_fmt = par_gets(block,"dat_fmt");
	goto add_it;
      }
      else if (strcmp(fmt,"hst")==0){
	new_out.fun = dump_history;
	if (par_exist(block,"dat_fmt"))
	  new_out.dat_fmt = par_gets(block,"dat_fmt");
	goto add_it;
      }
      else if (strcmp(fmt,"tab")==0){
	new_out.fun = dump_tab;
	if (par_exist(block,"dat_fmt"))
	  new_out.dat_fmt = par_gets(block,"dat_fmt");
	goto add_it;
      }
      else if (strcmp(fmt,"rst")==0){
	new_out.fun = dump_restart;
	add_rst_out(&new_out);
	continue;
      }
      else if (strcmp(fmt,"vtk")==0){
	new_out.fun = dump_vtk;
	goto add_it;
      }
      else{    /* Unknown data dump (fatal error) */
	ath_error("Unsupported dump mode for %s/out_fmt=%s for out=all\n",
          block,fmt);
      }
    }

/* Now handle data outputs (ouput of SINGLE variable).  There are lots more
 * options for outputs than dumps.  Need to choose variable, format, size
 * of domain to be output, scaling to min/max (if necessary),...    */

/* Is this a user defined expression? This allows the user to output any
 * problem-specific quantity using the formats and options supported here.
 * new_out.out must point to an expression defined in the user's problem.c */

    if(par_exist(block,"usr_expr_flag"))
      usr_expr_flag = par_geti(block,"usr_expr_flag");
    else
      usr_expr_flag = 0;

/* Get the expression function pointer */
    if(usr_expr_flag)
      new_out.expr = get_usr_expr(new_out.out);
    else
      new_out.expr = getexpr(nout, new_out.out);

    if (new_out.expr == NULL) {
      ath_perr(-1,"Could not parse expression %s, skipping it\n",
	      new_out.out);
      free_output(&new_out);
      continue;
    }

/* ix1, ix2, ix3:  parse and set size and coordinates of output "grid" */
    parse_slice(block,"ix1" ,pGrid->Nx1, pGrid->x1_0, pGrid->dx1,
		&new_out.ix1l, &new_out.ix1u, 
		&new_out.Nx1, &new_out.x1_0, &new_out.dx1);
    parse_slice(block,"ix2", pGrid->Nx2, pGrid->x2_0, pGrid->dx2,
		&new_out.ix2l, &new_out.ix2u, 
		&new_out.Nx2, &new_out.x2_0, &new_out.dx2);
    parse_slice(block,"ix3", pGrid->Nx3, pGrid->x3_0, pGrid->dx3,
		&new_out.ix3l, &new_out.ix3u, 
		&new_out.Nx3, &new_out.x3_0, &new_out.dx3);

/* notice the way how the ndim is increased while filling dim[] */
    new_out.ndim = 0;
    if (new_out.Nx1 > 1) new_out.ndim++;
    if (new_out.Nx2 > 1) new_out.ndim++;
    if (new_out.Nx3 > 1) new_out.ndim++;
    ath_pout(1,"DEBUG %d  ndim=%d  Nx=[%d %d %d]\n",
            new_out.n, new_out.ndim, new_out.Nx1,new_out.Nx2,new_out.Nx3);

/* dmin/dmax & sdmin/sdmax */
    if(par_exist(block,"dmin") != 0){ /* Use a fixed minimum scale? */
      new_out.sdmin = 1;
      new_out.dmin = par_getd(block,"dmin");
    }

    if(par_exist(block,"dmax") != 0){ /* Use a fixed maximum scale? */
      new_out.sdmax = 1;
      new_out.dmax = par_getd(block,"dmax");
    }

/* palette: default is rainbow */
    if (strcmp(fmt,"ppm") == 0) {
      new_out.palette = par_gets_def(block,"palette","rainbow");

      new_out.rgb = getRGB(new_out.palette);
      if ( (new_out.der = (float *) malloc(3*256*sizeof(float))) == NULL) {
	free_output(&new_out);
	ath_error("[init_output]: malloc returned a NULL pointer\n");
      }
      for(j=0; j<3; j++)    /* compute derivates to speed up interpolations */
	for (i=0; i<255; i++)
	  new_out.der[3*i+j] = new_out.rgb[3*(i+1)+j] - new_out.rgb[3*i+j];
    }

/* check for valid data output option (output of single variables)
 *  output format = {fits, pdf, pgm, ppm, tab}.  Note for pdf and tab outputs
 *  we also get the format for the print statements.
 */

    if (strcmp(fmt,"fits")==0)
      new_out.fun = output_fits;
    else if (strcmp(fmt,"pdf")==0){
      new_out.fun = output_pdf;
      if (par_exist(block,"dat_fmt"))
	new_out.dat_fmt = par_gets(block,"dat_fmt");
    }
    else if (strcmp(fmt,"pgm")==0)
      new_out.fun = output_pgm;
    else if (strcmp(fmt,"ppm")==0)
      new_out.fun = output_ppm;
    else if (strcmp(fmt,"tab")==0){
      new_out.fun = output_tab;
      if (par_exist(block,"dat_fmt"))
	new_out.dat_fmt = par_gets(block,"dat_fmt");
    }
    else {
/* unknown output format is fatal */
      free_output(&new_out);
      ath_error("Unsupported %s/out_fmt=%s\n",block,fmt);
    }

  add_it:

/* DEBUG */
    ath_pout(1,"OUTPUT: %d %d %s %s [%g : %g]\n",
             new_out.n, new_out.ndim, new_out.out_fmt,
             new_out.out, new_out.dmin, new_out.dmax);

    if(add_output(&new_out)) free_output(&new_out);
    else{
      nout++;
      ath_pout(0,"Added out%d\n",nout);
    }
  }
/*------------------------- end loop over output numbers ---------------------*/

}

/*----------------------------------------------------------------------------*/
/* data_output:  Called by main(), tests whether time for output, and calls
 *   appropriate output functions    */

void data_output(Grid *pGrid, Domain *pD, const int flag)
{
  int n;
  int dump_flag[MAXOUT_DEFAULT+1];
  char block[80];

/* Loop over all elements in output array
 * set dump flag to input argument, check whether time for output */
  for (n=0; n<out_count; n++) {
    dump_flag[n] = flag;
    if (pGrid->time >= OutArray[n].t) {
      OutArray[n].t += OutArray[n].dt;
      dump_flag[n] = 1;
    }
  }

/* Now check for restart dump, and make restart if dump_flag != 0 */
  if(rst_flag){
    dump_flag[out_count] = flag;
    if(pGrid->time >= rst_out.t){
      rst_out.t += rst_out.dt;
      dump_flag[out_count] = 1;
    }

    if(dump_flag[out_count] != 0){
/* Update the output numbers and times in the output blocks */
      for(n=0; n<out_count; n++){
/* User enrolled outputs have outn < 0 */
	if(OutArray[n].n > 0){
	  sprintf(block,"output%d",OutArray[n].n);
          if (dump_flag[n] != 0) {
/* About to write this output, so increase the output
 * number given in the restart file */
	    par_seti(block,"num","%d",OutArray[n].num+1,"Next Output Number");
          } else {
	    par_seti(block,"num","%d",OutArray[n].num,"Next Output Number");
          }
	  par_setd(block,"time","%.15e",OutArray[n].t,"Next Output Time");
	}
      }
/* Now do the same for the restart output block */
      sprintf(block,"output%d",rst_out.n);
      par_seti(block,"num","%d",rst_out.num+1,"Next Output Number");
      par_setd(block,"time","%.15e",rst_out.t,"Next Output Time");

/* Write the restart file */
      (*(rst_out.fun))(pGrid,pD,&(rst_out));

      rst_out.num++;
    }
  }

/* Loop over all elements in output array
 * If dump_flag != 0, make output */
  for (n=0; n<out_count; n++) {
    if(dump_flag[n] != 0) {
      (*OutArray[n].fun)(pGrid,pD,&(OutArray[n]));
      OutArray[n].num++;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* add_output: adds element initialized by init_output() to array of Outputs  */

int add_output(Output *new_out)
{
  size_t new_size;
  Output *new_array;

  if ((size_t)out_count >= out_size) {
/* Grow the output array by some small value, say 5 at a time */
    new_size = out_size + 5;
    if((new_array = realloc(OutArray,new_size*sizeof(Output))) == NULL){
      ath_error("[add_output]: Error growing output array\n");
    }
    OutArray = new_array;
    out_size = new_size;
  }

  OutArray[out_count] = *new_out;
  out_count++;

  return 0; /* Success */
}

/*----------------------------------------------------------------------------*/
/* add_rst: sets restart Output structure to value set by init_output() */

void add_rst_out(Output *new_out)
{
  rst_flag = 1;
  rst_out = *new_out;

  return;
}

/*----------------------------------------------------------------------------*/
/* data_output_destruct: free all memory associated with Output, called by
 *   main() at end of run */

void data_output_destruct(void)
{
  int i;

  for (i=0; i<out_count; i++) {
/* print the global min/max computed over the calculation */
    if (OutArray[i].out != NULL){
      if(strcmp(OutArray[i].out,"all") != 0)
	ath_pout(0,"Global min/max for %s: %g %g\n",OutArray[i].out,
		 OutArray[i].gmin, OutArray[i].gmax);

      free(OutArray[i].out);
    }
    if (OutArray[i].out_fmt != NULL) free(OutArray[i].out_fmt);
    if (OutArray[i].dat_fmt != NULL) free(OutArray[i].dat_fmt);
    if (OutArray[i].id      != NULL) free(OutArray[i].id);
  }

  if(rst_flag){
    if (rst_out.out     != NULL) free(rst_out.out);
    if (rst_out.out_fmt != NULL) free(rst_out.out_fmt);
    if (rst_out.dat_fmt != NULL) free(rst_out.dat_fmt);
    if (rst_out.id      != NULL) free(rst_out.id);
  }

  if (OutArray != NULL) {
    free(OutArray);
    OutArray = NULL;
    out_count = 0;
    out_size = 0;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* data_output_enroll: Enroll a user-defined data output function */

void data_output_enroll(Real time, Real dt, int num, const VGFunout_t fun,
			const char *fmt, const Gasfun_t expr, int n,
			const Real dmin, const Real dmax, int sdmin, int sdmax)
{
  Output new_out;

/* Zero (NULL) all members of the Output */
  memset(&new_out,0,sizeof(Output));

/* Set the input values and strdup the fmt string */
  new_out.n     = n;
  new_out.t     = time;
  new_out.dt    = dt;
  new_out.num   = num;
  new_out.dmin  = dmin;
  new_out.dmax  = dmax;
  new_out.sdmin = sdmin;
  new_out.sdmax = sdmax;
  new_out.fun   = fun;
  new_out.expr  = expr;

  if(fmt != NULL){
    if((new_out.out_fmt = ath_strdup(fmt)) == NULL)
      ath_perr(-1,"[dump_user_enroll]: Warning out_fmt strdup failed\n");
  }

  if(add_output(&new_out) && fmt != NULL) free(&(new_out.out_fmt));

  return;
}

/*----------------------------------------------------------------------------*/
/* subset3:  there is only one way to copy a cube into a cube */

float ***subset3(Grid *pgrid, Output *pout)
{
  float ***data;
  int Nx1, Nx2, Nx3;
  int i,j,k, il, jl, kl, iu, ju, ku;

  if (pout->ndim != 3)
    ath_error("[subset3] <output%d> %s has dimension %d, not 3\n",
              pout->n,pout->out, pout->ndim);

  Nx1 = pout->Nx1;
  Nx2 = pout->Nx2;
  Nx3 = pout->Nx3;

#ifdef WRITE_GHOST_CELLS
  if(pgrid->Nx1 > 1){
    il = pgrid->is - nghost;
    iu = pgrid->ie + nghost;
  }
  else{
    il = pgrid->is;
    iu = pgrid->ie;
  }

  if(pgrid->Nx2 > 1){
    jl = pgrid->js - nghost;
    ju = pgrid->je + nghost;
  }
 else{
    jl = pgrid->js;
    ju = pgrid->je;
  }
  if(pgrid->Nx3 > 1){
    kl = pgrid->ks - nghost;
    ku = pgrid->ke + nghost;
  }
  else{
    kl = pgrid->ks;
    ku = pgrid->ke;
  }
#else
  il = pgrid->is;
  iu = pgrid->ie;
  jl = pgrid->js;
  ju = pgrid->je;
  kl = pgrid->ks;
  ku = pgrid->ke;
#endif
  ath_pout(1,"subset2: lu's:  %d %d    %d %d     %d %d\n",
          il,  iu,  jl,  ju,  kl,  ku);

  data = (float ***) calloc_3d_array(Nx3,Nx2,Nx1,sizeof(float));
  for (k=0; k<Nx3; k++)
    for (j=0; j<Nx2; j++)
      for (i=0; i<Nx1; i++)
        data[k][j][i] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
  return data;
}

/*----------------------------------------------------------------------------*/
/* subset2:  this handles all 3 cases of subsetting a cube to a plane */

float **subset2(Grid *pgrid, Output *pout)
{
  float **data, factor;
  int Nx1, Nx2, Nx3;
  int i,j,k, il, jl, kl;

  if (pout->ndim != 2)
    ath_error("[subset2] <output%d> %s has dimension %d, not 2\n",
	      pout->n,pout->out, pout->ndim);

  Nx1 = pout->Nx1;
  Nx2 = pout->Nx2;
  Nx3 = pout->Nx3;

#ifdef WRITE_GHOST_CELLS
  if(pgrid->Nx1 > 1){
    il = pgrid->is - nghost;
  }
  else{
    il = pgrid->is;
  }

  if(pgrid->Nx2 > 1){
    jl = pgrid->js - nghost;
  }
  else{
    jl = pgrid->js;
  }
  if(pgrid->Nx3 > 1){
    kl = pgrid->ks - nghost;
  }
  else{
    kl = pgrid->ks;
  }
#else
  il = pgrid->is;
  jl = pgrid->js;
  kl = pgrid->ks;
#endif

/* handle TWO_DIM simulations */

  if (pgrid->Nx3 == 1) {
    data = (float **) calloc_2d_array(Nx2,Nx1,sizeof(float));
    for (j=0; j<Nx2; j++) {
      for (i=0; i<Nx1; i++) {
	data[j][i] = (*pout->expr)(pgrid,i+il,j+jl,0);
      }
    }
    return data;
  } else if (pgrid->Nx2 == 1  || pgrid->Nx2 == 1) {
    ath_error("subset2: cannot handle twodim with nx1 or nx2 being 1 yet");
  }
	  
  if (pout->Nx3 == 1) {
/* Nx3,Nx2,Nx1 -> Nx2,Nx1 */
    data = (float **) calloc_2d_array(Nx2,Nx1,sizeof(float));
    factor = 1.0/(pout->ix3u - pout->ix3l+1);
    for (j=0; j<Nx2; j++) {
      for (i=0; i<Nx1; i++) {
	data[j][i] = 0.0;
	for (k=pout->ix3l; k<=pout->ix3u; k++)
	  data[j][i] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
	data[j][i] *= factor;
      }
    }
  } else if (pout->Nx2 == 1) {
/* Nx3,Nx2,Nx1 -> Nx3,Nx1 */
    data = (float **) calloc_2d_array(Nx3,Nx1,sizeof(float));
    factor = 1.0/(pout->ix2u - pout->ix2l+1);
    for (k=0; k<Nx3; k++) {
      for (i=0; i<Nx1; i++) {
	data[k][i] = 0.0;
	for (j=pout->ix2l; j<=pout->ix2u; j++)
	  data[k][i] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
	data[k][i] *= factor;
      }
    }
  } else if (pout->Nx1 == 1) {
/* Nx3,Nx2,Nx1 -> Nx3,Nx2 */
    data = (float **) calloc_2d_array(Nx3,Nx2,sizeof(float));
    factor = 1.0/(pout->ix1u - pout->ix1l+1);
    for (k=0; k<Nx3; k++) {
      for (j=0; j<Nx2; j++) {
	data[k][j] = 0.0;
	for (i=pout->ix1l; i<=pout->ix1u; i++)
	  data[k][j] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
	data[k][j] *= factor;
      }
    }
  } else
    ath_perr(-1,"Should not reach here\n");
  return data;
}

/*----------------------------------------------------------------------------*/
/* subset1:  this handles all 3 cases of subsetting a cube to a vector */

float *subset1(Grid *pgrid, Output *pout)
{
  float *data, factor;
  int Nx1, Nx2, Nx3;
  int i,j,k, il, jl, kl;

  if (pout->ndim != 1)
    ath_error("[subset1] <output%d> %s has dimension %d, not 1\n",
	      pout->n,pout->out, pout->ndim);

  Nx1 = pout->Nx1;
  Nx2 = pout->Nx2;
  Nx3 = pout->Nx3;

#ifdef WRITE_GHOST_CELLS
  if(pgrid->Nx1 > 1){
    il = pgrid->is - nghost;
  }
  else{
    il = pgrid->is;
  }

  if(pgrid->Nx2 > 1){
    jl = pgrid->js - nghost;
  }
  else{
    jl = pgrid->js;
  }
  if(pgrid->Nx3 > 1){
    kl = pgrid->ks - nghost;
  }
  else{
    kl = pgrid->ks;
  }
#else
  il = pgrid->is;
  jl = pgrid->js;
  kl = pgrid->ks;
#endif

  if (pout->Nx3 > 1) {
/* Nx3,Nx2,Nx1 -> Nx3 */
    data = (float *) calloc_1d_array(Nx3,sizeof(float));
    factor = 1.0/(pout->ix2u - pout->ix2l+1)/(pout->ix1u - pout->ix1l+1);
    for (k=0; k<Nx3; k++) {
      data[k] = 0.0;
      for (j=0; j<Nx2; j++)
	for (i=0; i<Nx1; i++)
	  data[k] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
      data[k] *= factor;
    }
  } else if (pout->Nx2 > 1) {
/* Nx3,Nx2,Nx1 -> Nx2 */
    data = (float *) calloc_1d_array(Nx2,sizeof(float));
    factor = 1.0/(pout->ix3u - pout->ix3l+1)/(pout->ix1u - pout->ix1l+1);
    for (j=0; j<Nx2; j++) {
      data[j] = 0.0;
      for (k=0; k<Nx3; k++)
	for (i=0; i<Nx1; i++)
	  data[j] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
      data[j] *= factor;
    }
  } else if (pout->Nx1 > 1) {
/* Nx3,Nx2,Nx1 -> Nx1 */
    data = (float *) calloc_1d_array(Nx1,sizeof(float));
    factor = 1.0/(pout->ix3u - pout->ix3l+1)/(pout->ix2u - pout->ix2l+1);
    for (i=0; i<Nx1; i++) {
      data[i] = 0.0;
      for (k=0; k<Nx3; k++)
	for (j=0; j<Nx2; j++)
	  data[i] += (*pout->expr)(pgrid,i+il,j+jl,k+kl);
      data[i] *= factor;
    }
  } else {
    ath_perr(-1,"Should not reach here\n");
  }

  return data;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/* expr_*: where * are the conserved variables d,M1,M2,M3,E */

Real expr_d(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].d;
}
Real expr_M1(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M1;
}
Real expr_M2(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M2;
}
Real expr_M3(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M3;
}
#ifndef ISOTHERMAL
Real expr_E(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].E;
}
#endif

/*--------------------------------------------------------------------------- */
/* expr_*: where * are magnetic field variables: B1c, B2c, B3c, B^2 */

#ifdef MHD
Real expr_B1c(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].B1c;
}
Real expr_B2c(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].B2c;
}
Real expr_B3c(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].B3c;
}
Real expr_ME(const Grid *pG, const int i, const int j, const int k) {
  return 0.5*(pG->U[k][j][i].B1c*pG->U[k][j][i].B1c + 
	      pG->U[k][j][i].B2c*pG->U[k][j][i].B2c + 
	      pG->U[k][j][i].B3c*pG->U[k][j][i].B3c);
}
#endif

/*--------------------------------------------------------------------------- */
/* expr_*: where * are the primitive variables */

Real expr_V1(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M1/pG->U[k][j][i].d;
}
Real expr_V2(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M2/pG->U[k][j][i].d;
}
Real expr_V3(const Grid *pG, const int i, const int j, const int k) {
  return pG->U[k][j][i].M3/pG->U[k][j][i].d;
}

Real expr_P(const Grid *pG, const int i, const int j, const int k) {
#ifdef ISOTHERMAL
  return  pG->U[k][j][i].d*Iso_csound2;
#else
  Gas *gp = &(pG->U[k][j][i]);
  return Gamma_1*(gp->E 
#ifdef MHD
		  - 0.5*(gp->B1c*gp->B1c + gp->B2c*gp->B2c + gp->B3c*gp->B3c)
#endif /* MHD */
		  - 0.5*(gp->M1*gp->M1 + gp->M2*gp->M2 + gp->M3*gp->M3)/gp->d);
#endif /* ISOTHERMAL */
}

/*--------------------------------------------------------------------------- */
/* expr_cs2: sound speed squared  */

#ifndef ISOTHERMAL
Real expr_cs2(const Grid *pG, const int i, const int j, const int k)
{
  Gas *gp = &(pG->U[k][j][i]);
  return (Gamma*Gamma_1*(gp->E 
#ifdef MHD
	  - 0.5*(gp->B1c*gp->B1c + gp->B2c*gp->B2c + gp->B3c*gp->B3c)
#endif /* MHD */
	  - 0.5*(gp->M1*gp->M1 + gp->M2*gp->M2 + gp->M3*gp->M3)/gp->d)/gp->d);
}
#endif /* ISOTHERMAL */


/*--------------------------------------------------------------------------- */
/* expr_S: entropy = P/d^{Gamma}  */

#ifdef ADIABATIC
Real expr_S(const Grid *pG, const int i, const int j, const int k)
{
  Gas *gp = &(pG->U[k][j][i]);
  Real P = Gamma_1*(gp->E 
#ifdef MHD
		   - 0.5*(gp->B1c*gp->B1c + gp->B2c*gp->B2c + gp->B3c*gp->B3c)
#endif /* MHD */
		   - 0.5*(gp->M1*gp->M1 + gp->M2*gp->M2 + gp->M3*gp->M3)/gp->d);
  return P/pow((double)gp->d, (double)Gamma);
}
#endif /* ADIABATIC */

/*--------------------------------------------------------------------------- */
/* getexpr: return a function pointer for a simple expression - no parsing.
 *   For a user defined expression, get_usr_expr() in problem.c is used.  */

static Gasfun_t getexpr(const int n, const char *expr)
{
  char ename[32];

  sprintf(ename,"expr_out%d",n);

  if (strcmp(expr,"d")==0)
    return expr_d;
  else if (strcmp(expr,"M1")==0)
    return expr_M1;
  else if (strcmp(expr,"M2")==0)
    return expr_M2;
  else if (strcmp(expr,"M3")==0)
    return expr_M3;
#ifndef ISOTHERMAL
  else if (strcmp(expr,"E")==0)
    return expr_E;
#endif /* ISOTHERMAL */
#ifdef MHD
  else if (strcmp(expr,"B1c")==0)
    return expr_B1c;
  else if (strcmp(expr,"B2c")==0)
    return expr_B2c;
  else if (strcmp(expr,"B3c")==0)
    return expr_B3c;
  else if (strcmp(expr,"ME")==0)
    return expr_ME;
#endif
  else if (strcmp(expr,"V1")==0)
    return expr_V1;
  else if (strcmp(expr,"V2")==0)
    return expr_V2;
  else if (strcmp(expr,"V3")==0)
    return expr_V3;
  else if (strcmp(expr,"P")==0)
    return expr_P;
#ifndef ISOTHERMAL
  else if (strcmp(expr,"cs2")==0)
    return  expr_cs2;
#endif /* ISOTHERMAL */
#ifdef ADIABATIC
  else if (strcmp(expr,"S")==0)
    return  expr_S;
#endif /* ADIABATIC */
  else {
    FILE *fp = fopen("tmp.expr.c","w");
    fprintf(fp,"#include \"athena.h\"\n");
    fprintf(fp,"#include \"prototypes.h\"\n");
    fprintf(fp,"#include \"gasexpr.h\"\n");
    fprintf(fp,"Real %s(const Grid *pG, const int i, const int j, const int k) {\n",ename);
    fprintf(fp,"  return %s;\n",expr);
    fprintf(fp,"}\n");
    fclose(fp);
/* compile and link it to a shared object */
/* but for now just return NULL */
    ath_perr(-1,"Warning: dynamics expressions not supported yet\n");
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/* free_output: free memory associated with Output structure.  Only used when
 *   error occurs in adding a new output; this function frees memory and returns
 *   control to calling function */

static void free_output(Output *pOut)
{
  if(pOut->out     != NULL) free(pOut->out);
  if(pOut->out_fmt != NULL) free(pOut->out_fmt);
  if(pOut->dat_fmt != NULL) free(pOut->dat_fmt);
  if(pOut->id      != NULL) free(pOut->id);
  return;
}

/*----------------------------------------------------------------------------*/
/* parse_slice: sets the lower and upper bounds of a slice along an axis, 
 *   using values of ix1, ix2 or ix3 in the <output> block. The user can select
 *   bounds in the range 1..n (0 is not used for anything).  Valid formats are:
 *       ix1 = 5                a single value along x1 will be selected
 *       ix1 = 5:10             values 5 through 10 along x1 will be averaged
 *       ix1 = :                all planes along x1 (1:n) will be averaged
 *   and also
 *       ix1 = 5:               5:n actually
 *       ix1 = :10              1:10 actually
 *   If values for ix1,ix2,ix3 are not set in the <output> block, that
 *   dimension is not reduced, and the upper/lower bounds are set to -1.
 *   This function only parses the input text to extract the integer indices.
 *   the actual slicing and averaging is done by subset1,2().
 *
 *   The integer values of ix1,ix2,ix3 refer to the integer index of grid cells.
 *   If output of ghost cells is DISABLED (default), then ix1=5 refers to the
 *   fifth cell beyond the ghost zones (with index i=nghost+ix1).  If output
 *   of ghost cells is ENABLED (code configured with --enable-ghost), then
 *   ix1=5 refers to the fifth cell in the array (index i=ix1)
 *
 *   Also note that the user select an axis 1 based (i.e. 1..N) but they
 *   are stored 0 based (i.e. 0..N-1)
 */

static void parse_slice(char *block, char *axname, int nx, Real x, Real dx,
			int *l, int *u, int *new_nx, Real *new_x, Real *new_dx)
{
  char *expr, *cp;

  if (par_exist(block,axname)) {
    expr = par_gets(block,axname);
    cp = strchr(expr, ':');
    if (cp) {             /* either ':'  or 'lower:upper'  */
      *cp++ = 0;
      while (*cp && isspace(*cp)) cp++;
      if (*cp)
	*u = atoi(cp)-1;
      else
	*u = nx-1;
      cp = expr;
      while (*cp && isspace(*cp)) cp++;
      if (*cp)
	*l = atoi(cp)-1;
      else
	*l = 0;
    } else {               /* single slice  */
      *l = *u = atoi(expr)-1;
    }
    if (*l < 0 || *l >= nx) {
      ath_error("[parse_slice]: %s has illegal slice %d, not in range 1..%d\n",
	        expr,*l,nx);
    }
    free(expr);
    ath_pout(1,"DEBUG: parse_slice: %s/%s = %s =>  %d .. %d\n",
              block,axname,expr,*l,*u);
  } else {             /* no slicing in this axis */
    ath_pout(1,"DEBUG: parse_slice: %s/%s => no slicing \n",block,axname);
    *l = *u = -1;
  }
  if (*l == -1) {
    *new_nx = nx;
    *new_x  = x;
    *new_dx = dx;

/* notice that we don't use the FITS pixel centered definition of a WCS, but the
 * "left" edge of the cell. Thus the computational domain goes from x .. x+dx*nx
 */ 

  } else {
    *new_nx = (*u - *l + 1);
    *new_dx = dx;
    *new_x  = x + (*l+1)*dx ;
  }
}

/*----------------------------------------------------------------------------*/
/* getRGB: function for accessing palettes stored stored in structure RGB.
 *   Compares argument with strings (names) of palettes in RGB, and returns 
 *   pointer to first element of matching palette.  */

float *getRGB(char *name)
{
  int i;

  for (i=0; rgb[i].name && rgb[i].rgb; i++) {
    if (strcmp(name,rgb[i].name) == 0)
      return rgb[i].rgb;
  }

/* failed to find a matching palette: print them all and exit...  */

  ath_perr(-1,"Fatal error: could not find palette=%s, valid names are:\n",
    name);
  for (i=0; rgb[i].name && rgb[i].rgb; i++)
    ath_perr(-1,"%s ",rgb[i].name);
  ath_perr(-1,"\n");
  exit(EXIT_FAILURE);

  return NULL; /* Never reached, avoids compiler warnings */
}
