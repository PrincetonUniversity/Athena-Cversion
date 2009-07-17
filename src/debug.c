#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "debug.h"
#include "cyl.h"

// #ifdef CYLINDRICAL

/*============================================================================
 * DISPLAY  FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*  FUNCTION debug_header
 *
 *    THIS FUNCTION DISPLAYS A STANDARDIZED HEADER FOR VISUALLY SEPARATING OUT
 *    SECTIONS OF THE CODE.
 */
void debug_header(const int level, char *section_name)
{
  int k;
  char bar[BAR_LENGTH];

  for (k=0; k<BAR_LENGTH; k++) bar[k] = BAR_SYMBOL;
  ath_pout(level, "\n%s\n%s\n%s\n", bar,section_name,bar);
}


/*----------------------------------------------------------------------------*/
/*  FUNCTION print_cons1d
 *
 *    THIS FUNCTION PRINTS OUT THE CONTENTS OF A SINGLE Cons1D VECTOR AFTER
 *    ROTATING THE VECTOR COMPONENTS SO THAT THEY ARE LABELED CORRECTLY.
 */
int print_cons1d(const int level, char *name, const Cons1D *U,
                 const int k, const int j, const int i, const int dir)
{
  int n; 

  ath_pout(level, "\n%s[%d][%d][%d]\n", name,k,j,i);

  switch(dir) {
    case(1):  /* (X,Y,Z) */
              ath_pout(level, "\t | d  |   | " FMT " |\n", U->d);
              ath_pout(level, "\t | Mx |   | " FMT " |\n", U->Mx);
              ath_pout(level, "\t | My |   | " FMT " |\n", U->My);
              ath_pout(level, "\t | Mz |   | " FMT " |\n", U->Mz);
#ifndef ISOTHERMAL
              ath_pout(level, "\t | E  |   | " FMT " |\n", U->E);
#endif /* ISOTHERMAL */
#ifdef MHD
              ath_pout(level, "\t | By |   | " FMT " |\n", U->By);
              ath_pout(level, "\t | Bz |   | " FMT " |\n", U->Bz);
#endif /* MHD */
#if (NSCALARS > 0)
              for (n=0; n<NSCALARS; n++) {
              ath_pout(level, "\t | s%d  |   | " FMT " |\n", n,U->s[n]);
              }
#endif /* NSCALARS */
#if !defined(ISOTHERMAL) && defined(CYLINDRICAL)
              ath_pout(level, "\t | Pf |   | " FMT " |\n", U->Pflux);
#endif /* ISOTHERMAL */
              break;    

    case(2):  /* (Z,X,Y) */
              ath_pout(level, "\t | d  |   | " FMT " |\n", U->d);
              ath_pout(level, "\t | Mx |   | " FMT " |\n", U->Mz);
              ath_pout(level, "\t | My |   | " FMT " |\n", U->Mx);
              ath_pout(level, "\t | Mz |   | " FMT " |\n", U->My);
#ifndef ISOTHERMAL
              ath_pout(level, "\t | E  |   | " FMT " |\n", U->E);
#endif /* ISOTHERMAL */
#ifdef MHD
              ath_pout(level, "\t | Bx |   | " FMT " |\n", U->Bz);
              ath_pout(level, "\t | Bz |   | " FMT " |\n", U->By);
#endif /* MHD */
#if (NSCALARS > 0)
              for (n=0; n<NSCALARS; n++) {
              ath_pout(level, "\t | s%d  |   | " FMT " |\n", n,U->s[n]);
              }
#endif /* NSCALARS */
#if !defined(ISOTHERMAL) && defined(CYLINDRICAL)
              ath_pout(level, "\t | Pf |   | " FMT " |\n", U->Pflux);
#endif /* ISOTHERMAL */
              break;

    case(3):  /* (Y,Z,X) */
              ath_pout(level, "\t | d  |   | " FMT " |\n", U->d);
              ath_pout(level, "\t | Mx |   | " FMT " |\n", U->My);
              ath_pout(level, "\t | My |   | " FMT " |\n", U->Mz);
              ath_pout(level, "\t | Mz |   | " FMT " |\n", U->Mx);
#ifndef ISOTHERMAL
              ath_pout(level, "\t | E  |   | " FMT " |\n", U->E);
#endif /* ISOTHERMAL */
#ifdef MHD
              ath_pout(level, "\t | Bx |   | " FMT " |\n", U->By);
              ath_pout(level, "\t | By |   | " FMT " |\n", U->Bz);
#endif /* MHD */
#if (NSCALARS > 0)
              for (n=0; n<NSCALARS; n++) {
              ath_pout(level, "\t | s%d  |   | " FMT " |\n", n,U->s[n]);
              }
#endif /* NSCALARS */
#if !defined(ISOTHERMAL) && defined(CYLINDRICAL)
              ath_pout(level, "\t | Pf |   | " FMT " |\n", U->Pflux);
#endif /* ISOTHERMAL */
              break;

    default:  ath_error("[print_cons1d]:  Invalid direction!\n");
  }

  ath_pout(level, "\n");

  return 1;
}


/*----------------------------------------------------------------------------*/
/*  FUNCTION print_prim1d
 *
 *    THIS FUNCTION PRINTS OUT THE CONTENTS OF A SINGLE Prim1D VECTOR AFTER
 *    ROTATING THE VECTOR COMPONENTS SO THAT THEY ARE LABELED CORRECTLY.
 */
int print_prim1d(const int level, char *name, const Prim1D *W,
                 const int k, const int j, const int i, const int dir)
{
  int n;

  ath_pout(level, "\n%s[%d][%d][%d]\n", name,k,j,i);

  switch(dir) {
    case(1):  /* (X,Y,Z) */
              ath_pout(level, "\t | d  |   | " FMT " |\n", W->d);
              ath_pout(level, "\t | Vx |   | " FMT " |\n", W->Vx);
              ath_pout(level, "\t | Vy |   | " FMT " |\n", W->Vy);
              ath_pout(level, "\t | Vz |   | " FMT " |\n", W->Vz);
#ifndef ISOTHERMAL
              ath_pout(level, "\t | P  |   | " FMT " |\n", W->P);
#endif /* ISOTHERMAL */
#ifdef MHD
              ath_pout(level, "\t | By |   | " FMT " |\n", W->By);
              ath_pout(level, "\t | Bz |   | " FMT " |\n", W->Bz);
#endif /* MHD */
#if (NSCALARS > 0)
              for (n=0; n<NSCALARS; n++) {
              ath_pout(level, "\t | r%d  |   | " FMT " |\n", n,W->r[n]);
              }
#endif /* NSCALARS */
              break;    

    case(2):  /* (Z,X,Y) */
              ath_pout(level, "\t | d  |   | " FMT " |\n", W->d);
              ath_pout(level, "\t | Vx |   | " FMT " |\n", W->Vz);
              ath_pout(level, "\t | Vy |   | " FMT " |\n", W->Vx);
              ath_pout(level, "\t | Vz |   | " FMT " |\n", W->Vy);
#ifndef ISOTHERMAL
              ath_pout(level, "\t | P  |   | " FMT " |\n", W->P);
#endif /* ISOTHERMAL */
#ifdef MHD
              ath_pout(level, "\t | Bx |   | " FMT " |\n", W->Bz);
              ath_pout(level, "\t | Bz |   | " FMT " |\n", W->By);
#endif /* MHD */
#if (NSCALARS > 0)
              for (n=0; n<NSCALARS; n++) {
              ath_pout(level, "\t | r%d  |   | " FMT " |\n", n,W->r[n]);
              }
#endif /* NSCALARS */
              break;

    case(3):  /* (Y,Z,X) */
              ath_pout(level, "\t | d  |   | " FMT " |\n", W->d);
              ath_pout(level, "\t | Vx |   | " FMT " |\n", W->Vy);
              ath_pout(level, "\t | Vy |   | " FMT " |\n", W->Vz);
              ath_pout(level, "\t | Vz |   | " FMT " |\n", W->Vx);
#ifndef ISOTHERMAL
              ath_pout(level, "\t | P  |   | " FMT " |\n", W->P);
#endif /* ISOTHERMAL */
#ifdef MHD
              ath_pout(level, "\t | Bx |   | " FMT " |\n", W->By);
              ath_pout(level, "\t | By |   | " FMT " |\n", W->Bz);
#endif /* MHD */
#if (NSCALARS > 0)
              for (n=0; n<NSCALARS; n++) {
              ath_pout(level, "\t | r%d  |   | " FMT " |\n", n,W->r[n]);
              }
#endif /* NSCALARS */
              break;

    default:  ath_error("[print_prim1d]:  Invalid direction!\n");
  }

  ath_pout(level, "\n");

  return 1;
}


/*----------------------------------------------------------------------------*/
/*  FUNCTION print_gas
 *
 *    THIS FUNCTION PRINTS OUT THE CONTENTS OF A SINGLE Gas VECTOR.
 */
int print_gas(const int level, char *name, const Gas *U, const int k, const int j, const int i)
{
  int n;

  ath_pout(level, "\n%s[%d][%d][%d]\n", name,k,j,i);

  ath_pout(level, "\t | d  |   | " FMT " |\n", U->d);
  ath_pout(level, "\t | M1 |   | " FMT " |\n", U->M1);
  ath_pout(level, "\t | M2 |   | " FMT " |\n", U->M2);
  ath_pout(level, "\t | M3 |   | " FMT " |\n", U->M3);

#ifdef MHD
  ath_pout(level, "\t | B1 | = | " FMT " |\n", U->B1c);
  ath_pout(level, "\t | B2 |   | " FMT " |\n", U->B2c);
  ath_pout(level, "\t | B3 |   | " FMT " |\n", U->B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
  ath_pout(level, "\t | E  |   | " FMT " |\n", U->E);
#endif 
#if (NSCALARS > 0)
  for(n=0; n< NSCALARS; n++) {
  ath_pout(level, "\t | s%d |   | " FMT " |\n", n,U->s[n]);
  }
#endif

  ath_pout(level, "\n");

  return 1;
}

// #endif /* CYLINDRICAL */
