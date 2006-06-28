#include "copyright.h"
/*==============================================================================
 * FILE: ath_signal.c
 *
 * PURPOSE: Implements very simple signal handling.  Since signals can
 *   come in at any time, these functions set a static global
 *   variable which will be checked and reset if necessary -- TAG 8/19/2004
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   ath_sig_init()
 *   ath_sig_act()
 *
 * VARIABLE TYPE AND STRUCTURE DEFINITIONS: none.
 *============================================================================*/

#include <signal.h>
#include <stdio.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

static volatile int sig_caught = 0;   /* caught signal */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   handler()
 *============================================================================*/

static void handler(int s);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ath_sig_init:   */

void ath_sig_init(void)
{
  signal(SIGTERM, handler); /* Define the signal handler function */
  return;
}

/*----------------------------------------------------------------------------*/
/* ath_sig_act: Handles response to any received signals.  At the moment,
 *   only response to SIGTERM is implemented.   */

#ifdef MPI_PARALLEL
int ath_sig_act(Grid *pG)
{
  MPI_Status stat;
  int i, err, chsig, gsig;

  if(sig_caught == SIGTERM)
    fprintf(stderr,"[ath_sig_act]: Caught SIGTERM (my_id = %d)\n",pG->my_id);

  if(pG->my_id == 0){   /* I'm the parent */
    gsig = sig_caught;   /* Initialize global signal */

    for(i=1; i<pG->nproc; i++){   /* Collect the children's signal */
/* Receive a message from any child process */
      err = MPI_Recv(&chsig, 1, MPI_INT, MPI_ANY_SOURCE,
		     sync_step_tag, MPI_COMM_WORLD, &stat);
      if(err != MPI_SUCCESS)
	ath_error("[ath_sig_act]: MPI_Recv error = %d\n",err);

/* This logic can be extended in the future for different signals */
      gsig = (chsig == SIGTERM ? SIGTERM : gsig);
    }
  }
  else{   /* I'm a child */
    chsig = sig_caught;

/* Send the sig_caught value to the parent */
    err = MPI_Send(&chsig, 1, MPI_INT, 0, sync_step_tag, MPI_COMM_WORLD);
    if(err != MPI_SUCCESS)
      ath_error("[ath_sig_act]: MPI_Send error = %d\n",err);
  }

/* Now share the signal, chosen by the parent, with all of the children */
  err = MPI_Bcast(&gsig, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(err != MPI_SUCCESS)
    ath_error("[ath_sig_act]: Error on calling MPI_Bcast\n");

  sig_caught = 0; /* Reset the signal */

  if(gsig == SIGTERM){
    puts("Caught SIGTERM: Terminating program execution");
    return 1;
  }

  return 0;
}

#else

int ath_sig_act(Grid *pG)
{
  if(sig_caught == SIGTERM){
    puts("Caught SIGTERM: Terminating program execution");
    return 1;
  }
  return 0;
}

#endif /* MPI_PARALLEL */

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* handler:  */

static void handler(int s){
  sig_caught = s;
  signal(s, handler);   /* Reinstall the signal handler function */
  return;
}
