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

int ath_sig_act(int *piquit)
{

#ifdef MPI_PARALLEL
  int err, sig = sig_caught > *piquit ? sig_caught : *piquit;

  err = MPI_Allreduce(&sig, piquit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(err) sim_error("[sim_sig_act]: MPI_Allreduce returned error = %d\n",err);

#else /* SERIAL */

  *piquit = sig_caught > *piquit ? sig_caught : *piquit;

#endif /* MPI_PARALLEL */

  if(sig_caught == SIGTERM)
    ath_pout(0,"Caught SIGTERM: Terminating program execution\n");

  sig_caught = 0; /* Reset the signal */

  return *piquit;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* handler:  */

static void handler(int s){
  sig_caught = s;
  signal(s, handler);   /* Reinstall the signal handler function */
  return;
}
