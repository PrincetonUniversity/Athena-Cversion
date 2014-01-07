#include "copyright.h"
/*=============================================================================
 * FILE: units.c
 * PURPOSE: Contains a function to initialize a structure to hold the code
 *   units as well as several physical constants in those units.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *
 *   init_units()
 *   init_consts()
 *
 * History:
 *   Written by Aaron Skinner, Dec 2012
==============================================================================*/
#include <math.h>
#include "athena.h"
#include "prototypes.h"
#include "globals.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/

/*==============================================================================
 * PUBLIC FUNCTIONS:
 *============================================================================*/

/*-----------------------------CHANGE UNITS-----------------------------------*/
void init_units(UnitS *unit)
{
  Real cm,g,s,K,dyne,erg;
  ConstS consts;

  
  /* Set new base units */
  unit->Vcode = unit->Lcode/unit->Tcode;
  unit->Dcode = unit->Mcode/CUBE(unit->Lcode);
  cm          = 1.0/unit->Lcode;
  g           = 1.0/unit->Mcode;
  s           = 1.0/unit->Tcode;
  dyne        = g*cm/(s*s);
  erg         = g*cm*cm/(s*s);
  K           = 1.0;
  unit->cm    = cm;
  unit->g     = g;
  unit->s     = s;
  unit->dyne  = dyne;
  unit->erg   = erg;
  unit->K     = K;
  
  init_consts(&consts);

  /* Convert common physical constants/units from cgs to new units */
  unit->G     = consts.G * cm*cm*cm/(g*s*s);
  unit->Msun  = consts.Msun * g;
  unit->Lsun  = consts.Lsun * erg/s;
  unit->Myr   = consts.Myr * s;
  unit->pc    = consts.pc * cm;
  unit->kpc   = consts.kpc * cm;
  unit->kms   = consts.kms * cm/s;
  unit->mH    = consts.mH * g;
  unit->aR    = consts.aR * erg/(cm*cm*cm*K*K*K*K);
  unit->kB    = consts.kB * erg/K;
  unit->c     = consts.c * cm/s;
}

void init_consts(ConstS *consts)
{
  /* common physical constants in cgs */
  consts->G     = 6.67259e-8;
  consts->Msun  = 1.9891e+33;
  consts->Lsun  = 3.8268e+33;
  consts->Myr   = 3.155815e+13;
  consts->pc    = 3.085678e+18;
  consts->kpc   = 3.085678e+21;
  consts->kms   = 1.0e+5;
  consts->mH    = 1.6733e-24;
  consts->aR    = 7.5646e-15;
  consts->kB    = 1.380658e-16;
  consts->c     = 2.99792458e+10;
}
