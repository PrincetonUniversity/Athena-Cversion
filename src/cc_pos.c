#include "copyright.h"
/*============================================================================*/
/*! \file cc_pos.c
 *  \brief Functions to compute (x1,x2,x3) positions of cells i,j,k.  
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
 * - cc_pos() - given i,j,k returns cell-centered x1,x2,x3
 * - fc_pos() - given i,j,k returns face-centered x1,x2,x3
 * - x1cc() - given i, returns cell-centered x1.
 * - x2cc() - given j, returns cell-centered x2.
 * - x3cc() - given k, returns cell-centered x3.
 * - celli() - given x, returns containing cell first index. 
 * - cellj() - given y, returns containing cell first index.  
 * - cellk() - given y, returns containing cell first index.		      */
/*============================================================================*/

#include "athena.h"
#include "defs.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void cc_pos(const GridS *pG, const int i, const int j,const int k,
 *                  Real *px1, Real *px2, Real *px3)
 *  \brief given i,j,k returns cell-centered x1,x2,x3
 */
void cc_pos(const GridS *pG, const int i, const int j,const int k,
	    Real *px1, Real *px2, Real *px3)
{
  *px1 = pG->MinX[0] + ((Real)(i - pG->is) + 0.5)*pG->dx1;
  *px2 = pG->MinX[1] + ((Real)(j - pG->js) + 0.5)*pG->dx2;
  *px3 = pG->MinX[2] + ((Real)(k - pG->ks) + 0.5)*pG->dx3;
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void fc_pos(const GridS *pG, const int i, const int j,const int k,
 *                  Real *px1, Real *px2, Real *px3)
 *  \brief given i,j,k returns face-centered x1,x2,x3
 */

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
============================================================================*/

/*----------------------------------------------------------------------------*/
/* Input: pGrid: grid; x: global x coordinate;
 *        dx1_1: 1/dx1 (to improve performance)
 * Output: i: i-index containing x; a: grid index coordinate of x;
 * Return: 0: x is on the left of the ith cell;
 *         1: x is on the right of the ith cell;
 */

/*! \fn int celli(const GridS* pG, const Real x, const Real dx1_1, 
 *		  int *i, Real *a)
 *  \brief given x, returns containing cell first index.  */
int celli(const GridS* pG, const Real x, const Real dx1_1, int *i, Real *a)
{
  *a = (x - pG->MinX[0]) * dx1_1 + pG->is;
  *i = (int)(*a);
  if (((*a)-(*i)) < 0.5) return 0;	/* in the left half of the cell*/
  else return 1;			/* in the right half of the cell*/
}

/*! \fn Real x1cc(const GridS* pG, const int i)
 *  \brief given i, returns cell-centered x1. */
Real x1cc(const GridS* pG, const int i)
{
  return (pG->MinX[0] + ((Real)(i - pG->is) + 0.5) * pG->dx1);
}

/*! \fn cellj(const GridS* pG, const Real y, const Real dx2_1, 
 *	      int *j, Real *b)
 *  \brief given y, returns containing cell first index.  */
int cellj(const GridS* pG, const Real y, const Real dx2_1, int *j, Real *b)
{
  *b = (y - pG->MinX[1]) * dx2_1 + pG->js;
  *j = (int)(*b);
  if (((*b)-(*j)) < 0.5) return 0;	/* in the left half of the cell*/
  else return 1;			/* in the right half of the cell*/
}

/*! \fn Real x2cc(const GridS* pG, const int j)
 *  \brief given j, returns cell-centered x2. */
Real x2cc(const GridS* pG, const int j)
{
  return (pG->MinX[1] + ((Real)(j - pG->js) + 0.5) * pG->dx2);
}

/*! \fn int cellk(const GridS* pG, const Real z, const Real dx3_1, 
 *                int *k, Real *c)
 *  \brief given z, returns containing cell first index.  */ 
int cellk(const GridS* pG, const Real z, const Real dx3_1, int *k, Real *c)
{
  *c = (z - pG->MinX[2]) * dx3_1 + pG->ks;
  *k = (int)(*c);
  if (((*c)-(*k)) < 0.5) return 0;	/* in the left half of the cell*/
  else return 1;			/* in the right half of the cell*/
}

/*! \fn Real x3cc(const GridS* pG, const int k)
 *  \brief given k, returns cell-centered x3. */
Real x3cc(const GridS* pG, const int k)
{
  return (pG->MinX[2] + ((Real)(k - pG->ks) + 0.5) * pG->dx3);
}

#endif /* PARTICLES */
