#include "copyright.h"
/*==============================================================================
 * FILE: convert_var_sr.c
 *
 * PURPOSE: Functions to convert conservative to primitive vars, and vice versa,
 *   for special relativistic hydrodynamics and MHD..
 *
 * REFERENCE: S. Noble et al., "Primitive Variable Solvers for Conservative
 *   General Relativistic Magnetohydrodynamics", ApJ. 641, 626 (2006)
 *
 * HISTORY: First version written by Kevin Tian, summer 2007.
 *          Rewritten and revised by Jonathan Fulton, Senior Thesis 2009.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   Cons1D_to_Prim1D() - converts 1D vector (Bx passed through arguments)
 *   Prim1D_to_Cons1D() - converts 1D vector (Bx passed through arguments)
 * There are two versions of each function, one for hydro and MHD respectively.
 *============================================================================*/

#include <math.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef SPECIAL_RELATIVITY /* This version for special relativity only */

#ifdef MHD
/* The ratio (Gamma - 1) / Gamma where Gamma is the adiabatic constant */
Real Gamma_1overGamma;

/* prototypes for private functions needed with MHD */
static void calc_Function(Real Bsq, Real Msq, Real MdotB, Real MdotBsq, Real E,
   Real d, Real W, Real vsq, Real *deltaW, Real *deltavsq, Real *F1, Real *F2);
static Real calc_Vsq(Real Bsq, Real Msq, Real MdotB, Real MdotBsq, Real xW);
static int solve2D(Real Bsq, Real Msq, Real MdotB, Real MdotBsq, Real E, Real d,
   Real *xW, Real *xvsq);
#endif /* MHD */


#ifdef HYDRO
/*----------------------------------------------------------------------------*/
/* Cons1D_to_Prim1D: HYDRODYNAMICS VERSION
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 * Bx is passed in through the argument list.
 */

void Cons1D_to_Prim1D(const Cons1D *U, Prim1D *W MHDARG( , const Real *Bx))
{
  Real Msq, M, ME, Dsq, Gamma_1sq, denom;
  Real a3, a2, a1, a0;
  Real i1, i2, i3;
  Real iR, iS, iT;
  Real ix1=0.0, iB, iC, v, vOverM;

  /* Step One: Reduce equations to a quartic polynomial in v */

  Msq = SQR(U->Mx) + SQR(U->My) + SQR(U->Mz);
  M = sqrt(Msq);

  if (fabs(M) < TINY_NUMBER) {

    v = 0.0;   /* if M is zero, then v must be identically zero */

  } else {

    ME = M * U->E;
    Dsq = SQR(U->d);

    Gamma_1sq = SQR(Gamma_1);

    denom = 1.0 / (Gamma_1sq * (Msq + Dsq));

    a3 = (-2.0 * Gamma * Gamma_1 * ME) * denom;
    a2 = (SQR(Gamma) * SQR(U->E) + 2.0*Gamma_1*Msq - Gamma_1sq*Dsq)*denom;
    a1 = (-2.0 * Gamma * ME) * denom;
    a0 = Msq * denom;

    /* Step Two: Solve the polynomial analytically */

    i1 = -a2;
    i2 = a3 * a1 - 4.0 * a0;
    i3 = 4.0 * a2 * a0 - SQR(a1) - SQR(a3) * a0;

    iR = (9.0 * i1 * i2 - 27.0 * i3 - 2.0 * SQR(i1) * i1) / 54.0;
    iS = (3.0 * i2 - SQR(a2)) / 9.0;
    iT = SQR(iR) + SQR(iS) * iS;

    /* iT may be negative, but ix1 then adds a complex number and its conjugate,
     * giving us a real value, so ix1 is real no matter what */
    if (iT < 0) {
      ix1 = 2.0*pow(sqrt(iR*iR + iT),(ONE_3RD))*cos(atan2(sqrt(-iT),iR)/3.0)
             - i1/3.0;
    } else {
      ix1 = pow((iR + sqrt(iT)),(ONE_3RD)) + pow((iR - sqrt(iT)),(ONE_3RD)) 
             - i1/3.0;
    }

    iB = 0.5*(a3 + sqrt(SQR(a3) - 4.0*a2 + 4.0*ix1));
    iC = 0.5*(ix1 - sqrt(SQR(ix1) - 4.0*a0));
    v = (-iB + sqrt(SQR(iB) - 4.0*iC))/2.0;

  }

  /* Step Three: Resolve primitives with the value of v */

  /* ensure that v is physical */
  v = MAX(v, 0.0);
  v = MIN(v, 1.0 - 1.0e-15);

  vOverM = v / M;

  if (fabs(M) < TINY_NUMBER) {
    vOverM = 0.0;
  }

  W->d = sqrt(1.0 - SQR(v)) * U->d;

  W->Vx = U->Mx * vOverM;
  W->Vy = U->My * vOverM;
  W->Vz = U->Mz * vOverM;

  W->P = Gamma_1*((U->E - U->Mx*W->Vx - U->My*W->Vy - U->Mz*W->Vz) - W->d);

  return;
}
#endif /* HYDRO */

#ifdef HYDRO
/*----------------------------------------------------------------------------*/
/* Prim1D_to_Cons1D: HYDRODYNAMICS VERSION
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 * Bx is passed in through the argument list. 
 */

void Prim1D_to_Cons1D(Cons1D *U, const Prim1D *W MHDARG( , const Real *Bx))
{
  Real vsq = SQR(W->Vx) + SQR(W->Vy) + SQR(W->Vz);
  Real V0 = 1.0 / (1.0 - vsq);

  Real wV0sq = (W->d + W->P + W->P/Gamma_1) * V0;

  V0 = sqrt(V0);

  U->E = wV0sq - W->P;
  U->Mx = wV0sq * W->Vx;
  U->My = wV0sq * W->Vy;
  U->Mz = wV0sq * W->Vz;

  U->d = V0 * W->d;

  return;
}
#endif /* HYDRO */

#ifdef MHD
/*----------------------------------------------------------------------------*/
/* Cons1D_to_Prim1D: MHD VERSION 
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 * Bx is passed in through the argument list.
 *
 * IMPORTANT: This algorithm uses an iterative (Newton-Raphson) root-finding
 * step, which requires an initial guess for W.  This guess should be passed in
 * through W in the argument list, to be replaced with the calculated values.
 */

void Cons1D_to_Prim1D(const Cons1D *U, Prim1D *W, const Real *Bx)
{
  Real xW, xvsq, u0;
  Real xWplusBsq, MdotBoverxW;
  Real Bsq = 0.0, Msq = 0.0, MdotB = 0.0, MdotBsq = 0.0, E = 0.0, d = 0.0;

  int nr_success;
  int w_incs = 0;

  Gamma_1overGamma = Gamma_1/Gamma;
  Bsq = SQR((*Bx)) + SQR(U->By) + SQR(U->Bz);
  Msq = SQR(U->Mx) + SQR(U->My) + SQR(U->Mz);

  MdotB = U->Mx * (*Bx) + U->My * U->By + U->Mz * U->Bz;
  MdotBsq = SQR(MdotB);

  E = U->E;
  d = U->d;

  /* making a guess based on prior values for W */
  /* using a vsq that solves F1 is from Noble et al. */
  xvsq = SQR(W->Vx) + SQR(W->Vy) + SQR(W->Vz);
  xW = fabs((W->d + Gamma/Gamma_1 * W->P) / (1 - xvsq));
  xvsq = calc_Vsq(Bsq, Msq, MdotB, MdotBsq, xW);


/* Here, we ensure that the calculated value of xvsq in terms of xW is physical.
 *   In particular, we want xvsq < 1. Since xvsq is approximately proportional
 *   to 1/W^2, we can decrease xvsq by increasing xW.
 * However, we bound the number of decreases allowed because we may lose
 *   information otherwise.
 */
  while (xvsq >= 1.0 && w_incs < 15) {
    xW = xW * 10.0;
    xvsq = calc_Vsq(Bsq, Msq, MdotB, MdotBsq, xW);
    w_incs ++;
  }

/* In that offchance that we did not permit ourselves to inrease xW
 *   sufficiently, we still need to ensure the physicality of vsq.
 * The value (1 - 1.e-15) is from Noble et al. and should not cause overflow
 */
  xvsq = MIN(xvsq, 1.0 - 1.0e-15);

  nr_success = solve2D(Bsq, Msq, MdotB, MdotBsq, E, d, &xW, &xvsq);

/* If solve2D failed to converge, we make a second attempt. The only reparable
 *   failure I found was when vsq is close to 1, and the NR method attempts to
 *   increase vsq above 1. The attempt below tries to remedy this by decreasing
 *   vsq significantly and then allowing the NR method to raise it to the
 *   desired value.
 */
  if (nr_success != 0) {

    xvsq = SQR(W->Vx) + SQR(W->Vy) + SQR(W->Vz);
    xW = fabs((W->d + Gamma/Gamma_1 * W->P) / (1.0 - xvsq));
    xvsq = calc_Vsq(Bsq, Msq, MdotB, MdotBsq, xW);

    w_incs = 0;

/* The .001 threshold is fairly arbitrary */
    while (xvsq > 0.001 && w_incs < 30) {
      xW = xW * 10.0;
      xvsq = calc_Vsq(Bsq, Msq, MdotB, MdotBsq, xW);
      w_incs ++;
    }

    xvsq = MIN(xvsq, 1.0 - 1.0e-15); /* should not become useful */

    nr_success = solve2D(Bsq, Msq, MdotB, MdotBsq, E, d, &xW, &xvsq);

/* If we failed again, then we're out of ideas */
    if (nr_success != 0) {
      ath_error("[Cons1D_to_Prim1D]: NR iteration failed to converge\n");
    }

  }

  /* now hopefully we have answers */

  u0 = 1.0 / sqrt(1.0 - xvsq);
  xWplusBsq = 1.0 / (xW + Bsq);
  MdotBoverxW = MdotB / xW;

  W->d = U->d / u0;

  W->Vx = U->Mx;
  W->Vy = U->My;
  W->Vz = U->Mz;

  W->Vx += MdotBoverxW * (*Bx);
  W->Vy += MdotBoverxW * U->By;

  W->Vx *= xWplusBsq;
  W->Vy *= xWplusBsq;
  W->Vz *= xWplusBsq;

  W->P = (1.0 - xvsq) * Gamma_1 / Gamma * (xW - u0 * U->d);

  W->P = MAX(W->P, TINY_NUMBER); /* To ensure P is non-zero */

  W->By = U->By;
  W->Bz = U->Bz;

  return;
}
#endif /* MHD */

#ifdef MHD
/*----------------------------------------------------------------------------*/
/* Prim1D_to_Cons1D:  MHD VERSION
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 * Bx is passed in through the argument list.
 */

void Prim1D_to_Cons1D(Cons1D *U, const Prim1D *W, const Real *Bx)
{
  Real vsq, u0, Bsq = 0.0, vDotB = 0.0, wu0sq;

  vsq = SQR(W->Vx) + SQR(W->Vy) + SQR(W->Vz);
  u0 = 1.0 / sqrt(1.0 - vsq);

  Bsq = SQR(*Bx) + SQR(W->By) + SQR(W->Bz);

  vDotB = (*Bx)*W->Vx + W->By*W->Vy + W->Bz*W->Vz;

  wu0sq = (W->d + Gamma/Gamma_1 * W->P) * SQR(u0);

  U->d = u0 * W->d;

  U->Mx = (wu0sq  + Bsq) * W->Vx;
  U->My = (wu0sq  + Bsq) * W->Vy;
  U->Mz = (wu0sq  + Bsq) * W->Vz;

  U->E = wu0sq - W->P + (1.0 + vsq) * Bsq/2.0 - SQR(vDotB) / 2.0;

  U->Mx -= vDotB * (*Bx);
  U->My -= vDotB * W->By;
  U->Mz -= vDotB * W->Bz;

  U->By = W->By;
  U->Bz = W->Bz;

  return;
}
#endif /* MHD */

#ifdef MHD
/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*//* calcFunction:
 *   Calculate the values of F1, F2, deltaW, deltavsq given the other values.
 */
static void calc_Function(Real Bsq, Real Msq, Real MdotB, Real MdotBsq, Real E,
  Real d, Real W, Real vsq, Real *deltaW, Real *deltavsq, Real *F1, Real *F2) 
{
  Real J11,J12,J21,J22, detinv;
  Real BsqplusW = Bsq + W;
  Real BsqplusWqsq = SQR(BsqplusW);
  Real vsq_1 = 1.0 - vsq;
  Real doversqrtvsq_1 = d / sqrt(vsq_1);

/* Change these for whatever equation of state is desired */
/* We have an ideal equation of state below */
  Real P = Gamma_1overGamma * vsq_1 * (W - doversqrtvsq_1);
  Real dPdvsq = Gamma_1overGamma * (0.5 * doversqrtvsq_1);
  Real dPdW = Gamma_1overGamma * vsq_1;

  Real Wsq = SQR(W);
  Real Wcu = Wsq * W;

  *F1 = Msq - vsq * BsqplusWqsq + MdotBsq * (BsqplusW + W) / Wsq;
  *F2 = E - 0.5 * Bsq * (1.0 + vsq) + MdotBsq / (2.0 * Wsq) - W + P;

  J11 = -2.0 * BsqplusW * (vsq + MdotBsq / Wcu);
  J12 = - BsqplusWqsq; 
  J21 = -1.0 + dPdW - MdotBsq / Wcu; 
  J22 = -0.5 * Bsq - dPdvsq;

  detinv = 1.0 / (J11 * J22 - J12 * J21);

/* The NR step is
   [ deltaW   ] = -[ J11  J12 ]-1 [ *F1 ]
   [ deltavsq ]    [ J21  J22 ]   [ *F2 ]

   deltaF = -J^{-1} F
*/
  *deltaW = -(J22 * (*F1) - J12 * (*F2)) * detinv;
  *deltavsq = -(-J21 * (*F1) + J11 * (*F2)) * detinv;
}

/*----------------------------------------------------------------------------*/
/* calc_Vsq:
 *   Given Bsq, Msq, MdotB, MdotBsq, W, calculate the value of vsq
 */

static Real calc_Vsq(Real Bsq, Real Msq, Real MdotB, Real MdotBsq, Real xW) 
{
  Real BsqplusW = Bsq + xW;
  Real BsqplusWqsq = SQR(BsqplusW);
  Real Wsq = xW * xW;

  return (Msq * Wsq + MdotBsq * (Bsq + 2.0 * xW)) / (BsqplusWqsq * Wsq);
}

/*----------------------------------------------------------------------------*/
/* solve2D:
 *   Applies the Newton-Raphson method to numerically calculate the values of W
 *     and vsq given the values of Bsq, Msq, MdotB, MdotBsq, E, d, and guesses
 *     for W and vsq (in *xW and *xvsq)
 *   Places the values (if we find them) in *xW and *xvsq.
 *   Returns 0 if convergence was deemed succesful, non-zero otherwise.
 *     Returns 1 if W did not seem to be converging after MAX_CYCLES
 *     Returns 2 if vsq did not seem to be converging after MAX_CYCLES**
 *     Returns 3 if the guesses for W and/or vsq were unphysical (W < 0, vsq < 0
 *       or vsq >= 1)
 *   Convergence is measured by deltaW/W going to zero..
 */

/* The maximum number of NR cycles to perform per attempt, and the maximum
     number of extra cycles to run after deltaW/W is considered sufficiently
     small, respectively */
#define MAX_CYCLES 200
#define MAX_EXTRA 3

/* The permitted errors. We seek deltaW/W < ERROR_EXP, but consider it as
     converged if deltaW/W < ERROR_MIN.
   We also expect deltavsq/vsq < ERROR_VSQ */
#define ERROR_MIN 1.e-10
#define ERROR_EXP 1.e-16
#define ERROR_VSQ 1.e-10

static int solve2D(Real Bsq, Real Msq, Real MdotB, Real MdotBsq, Real E, Real d,
  Real *xW, Real *xvsq) 
{
  Real W = *xW, vsq = *xvsq;
  Real W_old, vsq_old;
  Real deltaW = 0.0, deltavsq = 0.0;
  Real F1,F2;
  Real errW, errvsq;
  int iterate = 1;
  int extra = 0;
  int cycles = 0;

  if (W < 0.0 || vsq < 0.0 || vsq >= 1.0) {
    return 3;
  }

  while (iterate == 1) {
    calc_Function(Bsq,Msq,MdotB,MdotBsq,E,d,W,vsq,&deltaW,&deltavsq,&F1,&F2);

    W_old = W;
    vsq_old = vsq;

    W += deltaW;
    vsq += deltavsq;

    errW = (fabs(W) < TINY_NUMBER) ? fabs(deltaW) : fabs(deltaW/W);

    W = fabs(W);  /* ensure that W is physical (W > 0) */
/* ensure that vsq is physical (0 <= vsq < 1) */
    vsq = MAX(vsq, 0.0);
    vsq = MIN(vsq, 1.0 - 1.0e-15); /* 1 - 1.e-15 is from Noble et al. */

    cycles ++;

    if (errW < ERROR_EXP) {
      extra ++;
    }

    if (cycles >= MAX_CYCLES) {
      iterate = 0;
    } else if (extra >= MAX_EXTRA) {
      iterate = 0;
    }

  }

  *xW = W;
  *xvsq = vsq;

/* check vsq error to catch traps near vsq == 1 */
  errvsq = (fabs(vsq) < TINY_NUMBER) ? fabs(deltavsq) : fabs(deltavsq/vsq);

/* so if the convergence is not sufficient, we return non-zero */
/* errvsq is a sign that the vsq value has problems */
  if (errvsq > ERROR_VSQ) {
    return 2;
  } else if (errW > ERROR_MIN) {
    return 1;
  } else if (errW <= ERROR_MIN && errW > ERROR_EXP) {
    return 0; /* in case we care to treat it differently */
  } else {
    return 0;
  }
}
#undef MAX_CYCLES
#undef MAX_EXTRA
#undef ERROR_MIN
#undef ERROR_EXP
#undef ERROR_VSQ

#endif /* MHD */
#endif /* SPECIAL_RELATIVITY */
