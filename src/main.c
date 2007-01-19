#include "copyright.h"
#define MAIN_C
/*==============================================================================
 * //////////////////////////// ATHENA Main Program \\\\\\\\\\\\\\\\\\\\\\\\\\\
 *
 *  Athena - Developed by JM Stone, TA Gardiner, PJ Teuben, & JF Hawley
 *  See the GNU General Public License for usage restrictions. 
 *
 *============================================================================*/
static char *athena_version = "version 3.0 - 01-JAN-2007";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* If there is a wall-time limit for MPI calculations run under a batch queue
 * system like PBS, then Athena will stop itself if there is less than
 * WTIME_WINDOW wall-seconds remaining to allow for time to write out the last
 * data files.
 */
#define WTIME_WINDOW 600.0

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   change_rundir - creates and outputs data to new directory
 *   usage         - outputs help message and terminates execution
 *============================================================================*/

void change_rundir(char *name);
void usage(char *prog);

/*----------------------------------------------------------------------------*/
/* main: Athena main program  */

int main(int argc, char *argv[])
{
  VGFun_t integrate;               /* pointer to integrator, set at runtime */
  Grid level0_Grid;                /* only level0 grids in this version */
  Domain level0_Domain;
  int ires=0;          /* restart flag, is set to 1 if -r argument on cmdline */
  int i,nlim,done=0,zones;
  Real tlim;
  double cpu_time, zcs;
  char *definput = "athinput", *rundir = NULL, *res_file = NULL;
  char *athinput = definput;
  long clk_tck = sysconf(_SC_CLK_TCK);
  struct tms tbuf;
  clock_t time0,time1, have_times;
  struct timeval tvs, tve;
#ifdef MPI_PARALLEL
  char *pc, *suffix, *name, new_name[MAXLEN];
  int len, h, m, s, err, use_wtlim=0;
  double wtstart, wtime, wtend;
#endif /* MPI_PARALLEL */

/* MPICH modifies argc and argv when calling MPI_Init() so that
 * after calling MPI_Init() only athena's command line arguments are
 * present and they are the same for the parent and child processes
 * alike. Note that this is not required by the MPI standard. */
#ifdef MPI_PARALLEL
  if(MPI_SUCCESS != MPI_Init(&argc, &argv))
    ath_error("[main]: Error on calling MPI_Init\n");
#endif /* MPI_PARALLEL */

/*--- Step 1. ----------------------------------------------------------------*/
/* Check for command line options and respond.  See comments in usage()
 * for description of options.  */

  for (i=1; i<argc; i++) {
/* If argv[i] is a 2 character string of the form "-?" then: */
    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'i':                                /* -i <file>   */
	athinput = argv[++i];
	break;
      case 'r':                                /* -r <file>   */
	ires = 1;
	res_file = argv[++i];
/* If the input file has not yet been set, use the restart file */
	if(athinput == definput) athinput = res_file;
	break;
      case 'd':                                /* -d <directory>   */
	rundir = argv[++i];
	break;
      case 'n':                                /* -n */
	done = 1;
	break;
      case 'h':                                /* -h */
	usage(argv[0]);
	break;
      case 'c':                                /* -c */
	show_config();
	exit(0);
	break;
#ifdef MPI_PARALLEL
      case 't':                                /* -t hh:mm:ss */
	use_wtlim = 1; /* Logical to use a wall time limit */
	sscanf(argv[++i],"%d:%d:%d",&h,&m,&s);
	wtend = s + 60*(m + 60*h);
	printf("Wall time limit: %d hrs, %d min, %d sec\n",h,m,s);
	wtstart = MPI_Wtime();
	break;
#endif /* MPI_PARALLEL */
#ifndef MPI_PARALLEL
      default:
	usage(argv[0]);
	break;
#endif /* MPI_PARALLEL */
      }
    }
  }

/*--- Step 2. (MPI_PARALLEL) -------------------------------------------------*/
/* For MPI_PARALLEL jobs, have parent read input file, parse command line, and
 * distribute information to children  */

#ifdef MPI_PARALLEL
/* Get my task id (rank in MPI) */
  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD,&(level0_Grid.my_id)))
    ath_error("Error on calling MPI_Comm_rank\n");

/* Get the number of processes */
  if(MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD,&(level0_Grid.nproc)))
    ath_error("Error on calling MPI_Comm_size\n");

/* Only parent (my_id==0) reads input parameter file, parses command line */
  if(level0_Grid.my_id == 0){
    par_open(athinput); 
    par_cmdline(argc,argv);
  }

/* Distribute the contents of the (updated) parameter file to the children. */
  par_dist_mpi(level0_Grid.my_id,MPI_COMM_WORLD);

/* Each child uses my_id in all output filenames to distinguish its output.
 * Output from the parent process does not have my_id in the filename */
  if(level0_Grid.my_id != 0){
    name = par_gets("job","problem_id");
    sprintf(new_name,"%s-id%d",name,level0_Grid.my_id);
    free(name);
    par_sets("job","problem_id",new_name,NULL);
  }

  show_config_par(); /* Add the configure block to the parameter database */

/* Share the restart flag with the children */
  if(MPI_SUCCESS != MPI_Bcast(&ires, 1, MPI_INT, 0, MPI_COMM_WORLD))
    ath_error("[main]: Error on calling MPI_Bcast\n");

/* Parent needs to send the restart file name to the children.  This requires 
 * sending the length of the restart filename string, the string, and then
 * having each child add my_id to the name so it opens the appropriate file */

/* Find length of restart filename */
  if(ires){ 
    if(level0_Grid.my_id == 0)
      len = 1 + (int)strlen(res_file);

/* Share this length with the children */
    if(MPI_SUCCESS != MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD))
      ath_error("[main]: Error on calling MPI_Bcast\n");

    if(len + 10 > MAXLEN)
      ath_error("[main]: Restart filename length = %d is too large\n",len);
    
/* Share the restart filename with the children */
    if(level0_Grid.my_id == 0) strcpy(new_name, res_file);
    if(MPI_SUCCESS != MPI_Bcast(new_name, len, MPI_CHAR, 0, MPI_COMM_WORLD))
      ath_error("[main]: Error on calling MPI_Bcast\n");

/* Assume the restart file name is of the form
       [/some/dir/elsewhere/]basename.0000.rst and search for the
       periods in the name. */
/* DOES THIS REQUIRE NAME HAVE 4 INTEGERS? */
    pc = &(new_name[len - 5]);
    if(*pc != '.'){
      ath_error("[main]: Bad Restart filename: %s\n",new_name);
    }

    do{ /* Position the char pointer at the first period */
      pc--;
      if(pc == new_name)
	ath_error("[main]: Bad Restart filename: %s\n",new_name);
    }while(*pc != '.');

/* Add my_id to the filename */
    if(level0_Grid.my_id == 0) {     /* I'm the parent */
      strcpy(new_name, res_file);
    } else {                         /* I'm a child */
      suffix = ath_strdup(pc);
      sprintf(pc,"-id%d%s",level0_Grid.my_id,suffix);
      free(suffix);
      res_file = new_name;
    }
  }

  if(done){          /* Quit MPI_PARALLEL job if code was run with -n option. */
    par_dump(0,stdout);   
    par_close();
    MPI_Finalize();
    return 0;
  }

/*--- Step 2. (MPI_SERIAL) ---------------------------------------------------*/
/* If this is *not* an MPI_PARALLEL job, there is only one process to open and
 * read input file  */

#else
  level0_Grid.my_id = 0;
  level0_Grid.nproc = 1;
  par_open(athinput);   /* opens AND reads */
  par_cmdline(argc,argv);
  show_config_par();   /* Add the configure block to the parameter database */

  if(done){      /* Quit non-MPI_PARALLEL job if code was run with -n option. */
    par_dump(0,stdout);
    par_close();
    return 0;
  }

#endif /* MPI_PARALLEL */

/*--- Step 3. ----------------------------------------------------------------*/
/* set variables in <time> block (these control execution time) */

  CourNo = par_getd("time","cour_no");
  nlim = par_geti_def("time","nlim",-1);
  tlim = par_getd("time","tlim");

/* Set variables in <job> block, EOS parameters from <problem> block.  */

  level0_Grid.outfilename = par_gets("job","problem_id");
#ifdef ISOTHERMAL
  Iso_csound = par_getd("problem","iso_csound");
  Iso_csound2 = Iso_csound*Iso_csound;
#else
  Gamma = par_getd("problem","gamma");
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;
#endif

/*--- Step 4. ----------------------------------------------------------------*/
/* Initialize and partition the computational Domain across processors.  Then
 * initialize each individual Grid block within Domain  */

  init_domain(&level0_Grid, &level0_Domain);
  init_grid  (&level0_Grid, &level0_Domain);

  if (level0_Grid.Nx1 > 1 && level0_Grid.Nx2 > 1 && level0_Grid.Nx3 > 1
    && CourNo >= 0.5) 
    ath_error("CourNo=%e , must be < 0.5 with 3D integrator\n",CourNo);

/*--- Step 5. ----------------------------------------------------------------*/
/* Set function pointers for integrate (based on dimensionality) */

  integrate = integrate_init(level0_Grid.Nx1,level0_Grid.Nx2,level0_Grid.Nx3);

/*--- Step 6. ----------------------------------------------------------------*/
/* Set initial conditions, either by reading from restart or calling
 * problem generator */

  if(ires) 
    restart_grid_block(res_file, &level0_Grid, &level0_Domain);  /*  Restart */
  else     
    problem(&level0_Grid, &level0_Domain);      /* New problem */

/*--- Step 7. ----------------------------------------------------------------*/
/* set boundary value function pointers using BC flags in <grid> blocks, then
 * set boundary conditions for initial conditions  */

  set_bvals_init(&level0_Grid, &level0_Domain);
  set_bvals(&level0_Grid);                            
  if(ires == 0) new_dt(&level0_Grid);

/*--- Step 8. ----------------------------------------------------------------*/
/* Set output modes (based on <ouput> blocks in input file).
 * Allocate temporary arrays needed by solver */

  init_output(&level0_Grid); 
  lr_states_init(level0_Grid.Nx1,level0_Grid.Nx2,level0_Grid.Nx3);

/*--- Step 9. ----------------------------------------------------------------*/
/* Setup complete, force dump of initial conditions with flag=1 */

  par_dump(0,stdout);      /* Dump a copy of the parsed information to stdout */
  change_rundir(rundir);
  if (ires==0) data_output(&level0_Grid, &level0_Domain, 1);

  ath_sig_init(); /* Install a signal handler */

  gettimeofday(&tvs,NULL);
  if((have_times = times(&tbuf)) > 0)
    time0 = tbuf.tms_utime + tbuf.tms_stime;
  else
    time0 = clock();

  if (level0_Grid.my_id == 0) {
    printf("\nSetup complete, entering main loop...\n\n");
    printf("cycle=%i time=%e dt=%e\n",level0_Grid.nstep,level0_Grid.time,
	 level0_Grid.dt);
  }

/*--- Step 10. ---------------------------------------------------------------*/
/* START OF MAIN INTEGRATION LOOP ==============================================
 * Steps are: (1) Check for data ouput
 *            (2) Integrate level0 grid
 *            (3) Set boundary values
 *            (4) Set new timestep
 */

  while (level0_Grid.time < tlim && (nlim < 0 || level0_Grid.nstep < nlim)) {

    data_output(&level0_Grid, &level0_Domain, 0);

/* modify timestep so loop finishes at t=tlim exactly */
    if ((tlim-level0_Grid.time) < level0_Grid.dt) {
      level0_Grid.dt = (tlim-level0_Grid.time);
    }

    (*integrate)(&level0_Grid);

    Userwork_in_loop(&level0_Grid, &level0_Domain);
    set_bvals(&level0_Grid);

    level0_Grid.nstep++;
    level0_Grid.time += level0_Grid.dt;
    new_dt(&level0_Grid);

#ifdef MPI_PARALLEL
    if(use_wtlim){
      if(level0_Grid.my_id == 0) /* Calculate the time remaining */
	wtime = wtend - (MPI_Wtime() - wtstart + WTIME_WINDOW);

      err = MPI_Bcast(&wtime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(wtime < 0.0 || err != MPI_SUCCESS) break;
    }
#endif /* MPI_PARALLEL */

    if(ath_sig_act(&level0_Grid) != 0) break;

    if (level0_Grid.my_id == 0) {
      printf("cycle=%i time=%e dt=%e\n",level0_Grid.nstep,level0_Grid.time,
              level0_Grid.dt);
    }
    fflush(stdout);
    fflush(stderr);

  }
/* END OF MAIN INTEGRATION LOOP ==============================================*/

/*--- Step 11. ---------------------------------------------------------------*/
/* Finish up by computing zc/sec, dumping data, and deallocate memory */

  if (level0_Grid.nstep == nlim)
    printf("\nterminating on cycle limit\n");
#ifdef MPI_PARALLEL
  else if(use_wtlim && (wtime < 0.0))
    printf("\nterminating on wall-time limit\n");
#endif /* MPI_PARALLEL */
  else
    printf("\nterminating on time limit\n");

  gettimeofday(&tve,NULL);
  if(have_times > 0) {
    times(&tbuf);
    time1 = tbuf.tms_utime + tbuf.tms_stime;
    cpu_time = (time1 > time0 ? (double)(time1 - time0) : 1.0)/
      (double)clk_tck;
  } else {
    time1 = clock();
    cpu_time = (time1 > time0 ? (double)(time1 - time0) : 1.0)/
      (double)CLOCKS_PER_SEC;
  }

  zones = level0_Grid.Nx1*level0_Grid.Nx2*level0_Grid.Nx3;
  zcs = (double)zones*(double)(level0_Grid.nstep)/cpu_time;

/* Calculate and print the zone-cycles / cpu-second */
  if (level0_Grid.my_id == 0) {
    printf("  tlim= %e   nlim= %i\n",tlim,nlim);
    printf("  time= %e  cycle= %i\n",level0_Grid.time,level0_Grid.nstep);
    printf("\nzone-cycles/cpu-second = %e\n",zcs);
  }

/* Calculate and print the zone-cycles / wall-second */
  cpu_time = (double)(tve.tv_sec - tvs.tv_sec) +
    1.0e-6*(double)(tve.tv_usec - tvs.tv_usec);
  zcs = (double)zones*(double)(level0_Grid.nstep)/cpu_time;
  if (level0_Grid.my_id == 0) {
    printf("\nelapsed wall time = %e sec.\n",cpu_time);
    printf("\nzone-cycles/wall-second = %e\n",zcs);
  }

#ifdef MPI_PARALLEL
  zones = (level0_Domain.ixe - level0_Domain.ixs + 1)
         *(level0_Domain.jxe - level0_Domain.jxs + 1)
         *(level0_Domain.kxe - level0_Domain.kxs + 1);
  zcs = (double)zones*(double)(level0_Grid.nstep)/cpu_time;
  if (level0_Grid.my_id == 0) {
    printf("\ntotal zone-cycles/wall-second = %e\n",zcs);
  }
#endif /* MPI_PARALLEL */

/* complete any final User work, and make last dump */

  Userwork_after_loop(&level0_Grid, &level0_Domain);
  data_output(&level0_Grid, &level0_Domain, 1);     /* always write last data */

/* Free all temporary arrays */

  lr_states_destruct();
  integrate_destruct();
  data_output_destruct();
  par_close();       

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif /* MPI_PARALLEL */

  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*  change_rundir: change run directory;  create it if it does not exist yet
 */

void change_rundir(char *name)
{
  int err=0;
#ifdef MPI_PARALLEL
  int rerr, gerr;
#endif /* MPI_PARALLEL */

  if (name == NULL || *name == 0) return;
  fprintf(stderr,"Changing run directory to \"%s\"\n",name);

  if (chdir(name)) {
    if (mkdir(name,0775)) {
      fprintf(stderr,"Failed to create run directory %s\n",name);
      err = 1;
    }
    else if (chdir(name)) {
      fprintf(stderr,"Cannot change directory to %s\n",name);
      err = 1;
    }
  }

#ifdef MPI_PARALLEL

  rerr = MPI_Allreduce(&err, &gerr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(rerr) fprintf(stderr,"[change_rundir]: MPI_Allreduce error = %d\n",rerr);

  if(rerr || gerr){
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

#else

  if(err) exit(EXIT_FAILURE);

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*  usage: outputs help
 *    athena_version is hardwired at beginning of this file
 *    CONFIGURE_DATE is macro set when configure script runs */

void usage(char *prog)
{
  fprintf(stderr,"Athena %s\n",athena_version);
  fprintf(stderr,"  Last configure: %s\n",CONFIGURE_DATE);
  fprintf(stderr,"\nUsage: %s [options] [block/par=value ...]\n",prog);
  fprintf(stderr,"\nOptions:\n");
  fprintf(stderr,"  -i <file>       Alternate input file [athinput]\n");
  fprintf(stderr,"  -d <directory>  Alternate run dir [current dir]\n");
  fprintf(stderr,"  -h              This Help, and configuration settings\n");
  fprintf(stderr,"  -n              Parse input, but don't run program\n"); 
  fprintf(stderr,"  -c              Show Configuration details and quit\n"); 
  fprintf(stderr,"  -r <file>       Restart a simulation with this file\n");
  show_config();
  exit(0);
}
