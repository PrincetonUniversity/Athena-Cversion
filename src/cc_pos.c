#include "copyright.h"
/*==============================================================================
 * FILE: cc_pos.c
 *
 * PURPOSE: Functions to compute (x1,x2,x3) positions of cells i,j,k.  
 *   In a nested grid, each Grid structure is a patch in a larger computational
 *   domain (with the exception of the level0 grid).  The displacement of the
 *   origin of the Grid from the origin of the computational domain (level0 
 *   grid) is x1_{disp} = idisp*dx1.  Furthermore, the origin of the level0
 *   grid can be displaced by a distance x1_{0} from the origin of x1.  Thus,
 *   the x1 position of the center of cell i (x1_{cc,i}) in any level Grid is
 *            x1_{cc,i} = x1_{0} + ((i + idisp) + 0.5)*dx1
 *   Similarly for x2 and x3.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   cc_pos() - given i,j,k returns cell-centered x1,x2,x3
 *   fc_pos() - given i,j,k returns face-centered x1,x2,x3
 *   ccr_pos() - given i,j,k returns cell-centered x1,x2,x3 for RadGridS
 *   x1cc - given i, returns cell-centered x1. Same for x2cc and x3cc
 *   celli - given x, returns containing cell first index. Same holds for cellj  *and cellk
 *============================================================================*/

#include "athena.h"
#include "defs.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* cc_pos:  */

void cc_pos(const GridS *pG, const int i, const int j,const int k,
	    Real *px1, Real *px2, Real *px3)
{
  *px1 = pG->MinX[0] + ((Real)(i - pG->is) + 0.5)*pG->dx1;
  *px2 = pG->MinX[1] + ((Real)(j - pG->js) + 0.5)*pG->dx2;
  *px3 = pG->MinX[2] + ((Real)(k - pG->ks) + 0.5)*pG->dx3;
  return;
}

/*----------------------------------------------------------------------------*/
/* fc_pos:  */

void fc_pos(const GridS *pG, const int i, const int j,const int k,
	    Real *px1, Real *px2, Real *px3)
{
  *px1 = pG->MinX[0] + ((Real)(i - pG->is))*pG->dx1;
  *px2 = pG->MinX[1] + ((Real)(j - pG->js))*pG->dx2;
  *px3 = pG->MinX[2] + ((Real)(k - pG->ks))*pG->dx3;
  return;
}

#ifdef CYLINDRICAL
Real x1vc(const GridS* pG, const int i)
{
  Real x1cc = pG->MinX[0] + ((Real)(i - pG->is) + 0.5)*pG->dx1;
  return x1cc + SQR(pG->dx1)/(12.0*x1cc);
}
#endif

#ifdef PARTICLES
/*============================================================================
cell-location functions 
Created: Emmanuel Jacquet, Mar. 2008
Modified: Xuening Bai, Dec. 2008
============================================================================*/

/*----------------------------------------------------------------------------*/
/* Input: pGrid: grid; x: global x coordinate;
 *        dx1_1: 1/dx1 (to improve performance)
 * Output: i: i-index containing x; a: grid index coordinate of x;
 * Return: 0: x is on the left of the ith cell;
 *         1: x is on the right of the ith cell;
 */

int celli(const GridS* pGrid, const Real x, const Real dx1_1, int *i, Real *a)
{
  *a = (x - pGrid->x1_0) * dx1_1 - pGrid->idisp;
  *i = (int)(*a);
  if (((*a)-(*i)) < 0.5) return 0;	/* in the left half of the cell*/
  else return 1;			/* in the right half of the cell*/
}

Real x1cc(const Grid* pGrid, const int i)
{
  return (pGrid->x1_0 + (i + pGrid->idisp + 0.5) * pGrid->dx1);
}

int cellj(const Grid* pGrid, const Real y, const Real dx2_1, int *j, Real *b)
{
  *b = (y - pGrid->x2_0) * dx2_1 - pGrid->jdisp;
  *j = (int)(*b);
  if (((*b)-(*j)) < 0.5) return 0;	/* in the left half of the cell*/
  else return 1;			/* in the right half of the cell*/
}

Real x2cc(const Grid* pGrid, const int j)
{
  return (pGrid->x2_0 + (j + pGrid->jdisp + 0.5) * pGrid->dx2);
}

int cellk(const Grid* pGrid, const Real z, const Real dx3_1, int *k, Real *c)
{
  *c = (z - pGrid->x3_0) * dx3_1 - pGrid->kdisp;
  *k = (int)(*c);
  if (((*c)-(*k)) < 0.5) return 0;	/* in the left half of the cell*/
  else return 1;			/* in the right half of the cell*/
}

Real x3cc(const Grid* pGrid, const int k)
{
  return (pGrid->x3_0 + (k + pGrid->kdisp + 0.5) * pGrid->dx3);
}

#endif /* PARTICLES */

#ifdef RADIATION_TRANSFER
/*----------------------------------------------------------------------------*/
/* ccr_pos:  */

void ccr_pos(const RadGridS *pRG, const int i, const int j,const int k,
	    Real *px1, Real *px2, Real *px3)
{
  *px1 = pRG->MinX[0] + ((Real)(i - pRG->is) + 0.5)*pRG->dx1;
  *px2 = pRG->MinX[1] + ((Real)(j - pRG->js) + 0.5)*pRG->dx2;
  *px3 = pRG->MinX[2] + ((Real)(k - pRG->ks) + 0.5)*pRG->dx3;
  return;
}
#endif /* RADIATION_TRANSFER */
