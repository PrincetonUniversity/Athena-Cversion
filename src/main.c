#include "copyright.h"
/*==============================================================================
 * //////////////////////////// ATHENA Main Program \\\\\\\\\\\\\\\\\\\\\\\\\\\
 *
 *  Athena - Developed by JM Stone, TA Gardiner, PJ Teuben, & JF Hawley
 *  See the GNU General Public License for usage restrictions. 
 *
 *============================================================================*/
static char *athena_version = "version 3.0 - XX-XXX-2006";

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
#include "prototypes.h"

/* Maximum length for a basename, or path + restart file name in the
 * case of an MPI Parallel calculation. */
#define MAXLEN 256

/* If there is a wall-time limit for MPI calculations run under a batch queue
 * system like PBS, then Athena will stop itself if there is less than
 * WTIME_WINDOW wall-seconds remaining to allow for time to write out the last
 * data files.
 */
#define WTIME_WINDOW 60.0

#ifdef MPI_PARALLEL
/* Min and Max bounding coordinates of the Complete Computational Grid. */
extern int cg_ixs, cg_jxs, cg_kxs;
extern int cg_ixe, cg_jxe, cg_kxe;
#endif /* MPI_PARALLEL */

Real CourNo; /* The Courant, Friedrichs, & Lewy (CFL) Number */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   change_rundir - creates and outputs data to new directory
 *   usage         - outputs help message and terminates execution
 *============================================================================*/

void change_rundir(char *name);
void usage(char *prog);

/*----------------------------------------------------------------------------*/
/* main:  main program  */

int main(int argc, char *argv[])
{
  VGFun_t integrate;               /* pointer to integrator, set at runtime */
  Grid grid_level0;                /* only level0 grids in this version */
  int i,ires=0,nlim,done=0,zones;
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

/*----------------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------------*/
/* For MPI_PARALLEL jobs, have parent read input file and distribute
 * information to children  */

#ifdef MPI_PARALLEL
/* Get my task id, or rank as it is called in MPI */
  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD,&(grid_level0.my_id)))
    ath_error("Error on calling MPI_Comm_rank\n");

/* Get the number of processes */
  if(MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD,&(grid_level0.nproc)))
    ath_error("Error on calling MPI_Comm_size\n");

/* Only parent (my_id==0) reads input parameter file, parses command line */
  if(grid_level0.my_id == 0){
    par_open(athinput); 
    par_cmdline(argc,argv);
  }

/* Distribute the contents of the (updated) parameter file to the children. */
  par_dist_mpi(grid_level0.my_id,MPI_COMM_WORLD);

/* Each child uses my_id in all output filenames to distinguish its output.
 * Output from the parent process does not have my_id in the filename */
  if(grid_level0.my_id != 0){
    name = par_gets("job","problem_id");
    sprintf(new_name,"%s-id%d",name,grid_level0.my_id);
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
    if(grid_level0.my_id == 0)
      len = 1 + (int)strlen(res_file);

/* Share this length with the children */
    if(MPI_SUCCESS != MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD))
      ath_error("[main]: Error on calling MPI_Bcast\n");

    if(len + 10 > MAXLEN)
      ath_error("[main]: Restart filename length = %d is too large\n",len);

/* Share the restart filename with the children */
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
    if(grid_level0.my_id == 0) {     /* I'm the parent */
      strcpy(new_name, res_file);
    } else {                         /* I'm a child */
      suffix = ath_strdup(pc);
      sprintf(pc,"-id%d%s",grid_level0.my_id,suffix);
      free(suffix);
      res_file = new_name;
    }
  }

  if(done){           /* Quit MPI_PARALLEL job if code was run with -n option. */
    par_dump(0,stdout);   
    par_close();
    MPI_Finalize();
    return 0;
  }

#else

/*----------------------------------------------------------------------------*/
/* If this is *not* an MPI_PARALLEL job, there is only one process to open and
 * read input file  */

  par_open(athinput);   /* opens AND reads */
  par_cmdline(argc,argv);
  show_config_par();   /* Add the configure block to the parameter database */

  if(done){       /* Quit non-MPI_PARALLEL job if code was run with -n option. */
    par_dump(0,stdout);
    par_close();
    return 0;
  }

#endif /* MPI_PARALLEL */

/*----------------------------------------------------------------------------*/
/* set variables in <time> block (these control execution time) */

  CourNo = par_getd("time","cour_no");
  nlim = par_geti_def("time","nlim",-1);
  tlim = par_getd("time","tlim");

/* Set variables in <job> block.  */

  grid_level0.outfilename = par_gets("job","problem_id");
  set_eos_param();

/* For both new and restart runs, the grid is initialized from parameters that
 * have already been read from the input file */

  init_grid_block(&grid_level0);

/* Set function pointers for integrate (based on dimensionality) */

  integrate = integrate_init(grid_level0.Nx1,grid_level0.Nx2,grid_level0.Nx3);

/* Set initial conditions, either by reading from restart or calling
 * problem generator */

  if(ires) 
    restart_grid_block(res_file,&grid_level0);  /*  Restart */
  else     
    problem(&grid_level0);                      /* New problem */

/* set boundary value function pointers using BC flags, read from <grid> blocks */
  set_bvals_init(&grid_level0);
  set_bvals(&grid_level0);                                 /* in input file */
  if(ires == 0) init_dt(&grid_level0);

/* Set output modes (based on <ouput> blocks in input file).
 * Allocate temporary arrays needed by solver */

  init_output(&grid_level0); 
  lr_states_init(grid_level0.Nx1,grid_level0.Nx2,grid_level0.Nx3);

/* Setup complete, dump initial conditions */

  par_dump(0,stdout); /* Dump a copy of the parsed information to stdout */
  change_rundir(rundir);

  ath_sig_init(); /* Install a signal handler */

  Userwork_before_loop(&grid_level0);
  gettimeofday(&tvs,NULL);
  if((have_times = times(&tbuf)) > 0)
    time0 = tbuf.tms_utime + tbuf.tms_stime;
  else
    time0 = clock();

  printf("\nSetup complete, entering main loop...\n\n");
  printf("cycle=%i time=%e dt=%e\n",grid_level0.nstep,grid_level0.time,
	 grid_level0.dt);

/*----------------------------------------------------------------------------*/
/* START OF MAIN INTEGRATION LOOP ==============================================
 * Steps are: (1) Check for data ouput
 *            (2) Integrate level0 grid
 *            (3) Set boundary values
 */

  while (grid_level0.time < tlim && (nlim < 0 || grid_level0.nstep < nlim)) {

    data_output(&grid_level0, 0);

    if ((tlim-grid_level0.time) < grid_level0.dt) {
      grid_level0.dt = (tlim-grid_level0.time);
    }

    printf("cycle=%i time=%e dt=%e\n",grid_level0.nstep,grid_level0.time,
	   grid_level0.dt);

    (*integrate)(&grid_level0);

#ifdef MPI_PARALLEL
    sync_dt(&grid_level0);
#endif /* MPI_PARALLEL */

    set_bvals(&grid_level0);

#ifdef MPI_PARALLEL
    if(use_wtlim){
      if(grid_level0.my_id == 0) /* Calculate the time remaining */
	wtime = wtend - (MPI_Wtime() - wtstart + WTIME_WINDOW);

      err = MPI_Bcast(&wtime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(wtime < 0.0 || err != MPI_SUCCESS) break;
    }
#endif /* MPI_PARALLEL */

    if(ath_sig_act(&grid_level0) != 0) break;
  }
/* END OF MAIN INTEGRATION LOOP ===============================================*/
/*------------------------------------------------------------------------------
 * Finish up by computing zc/sec, dumping data, and deallocate memory */

  if (grid_level0.nstep == nlim)
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

  zones = grid_level0.Nx1*grid_level0.Nx2*grid_level0.Nx3;
  zcs = (double)zones*(double)(grid_level0.nstep)/cpu_time;

/* Calculate and print the zone-cycles / cpu-second */
  printf("  tlim= %e   nlim= %i\n",tlim,nlim);
  printf("  time= %e  cycle= %i\n",grid_level0.time,grid_level0.nstep);
  printf("\nzone-cycles/cpu-second = %e\n",zcs);

/* Calculate and print the zone-cycles / wall-second */
  cpu_time = (double)(tve.tv_sec - tvs.tv_sec) +
    1.0e-6*(double)(tve.tv_usec - tvs.tv_usec);
  printf("\nelapsed wall time = %e sec.\n",cpu_time);
  zcs = (double)zones*(double)(grid_level0.nstep)/cpu_time;
  printf("\nzone-cycles/wall-second = %e\n",zcs);

#ifdef MPI_PARALLEL
  zones = (cg_ixe - cg_ixs + 1)*(cg_jxe - cg_jxs + 1)*(cg_kxe - cg_kxs + 1);
  zcs = (double)zones*(double)(grid_level0.nstep)/cpu_time;
  printf("\ntotal zone-cycles/wall-second = %e\n",zcs);
#endif /* MPI_PARALLEL */

/* complete any final User work, and make last dump */

  Userwork_after_loop(&grid_level0);
  data_output(&grid_level0, 1);   /* always write the last data */

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
