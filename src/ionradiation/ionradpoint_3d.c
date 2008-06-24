#include "copyright.h"
/*==============================================================================
 * FILE: ionrad_3d.c
 *
 * PURPOSE: Contains functions to compute the photoionization rate
 *   from point sources, using the algorithm described in Krumholz,
 *   Stone, & Gardiner (2007).
 *
 *   Use of these routines requires that --enable-ion-radpoint be set
 *   at compile time.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   add_radpoint_3d             - adds a new radiation source
 *   restore_radpoint_3d         - restores a radiation source from
 *                                    a checkpoint
 *   ion_radpoint_get_ranstate   - returns state of random number generator
 *                                    so it can be stored in checkpoints
 *   ion_radpoint_set_ranstate   - sets the state of the random number
 *                                    generator; used on restarts
 *   ion_radpoint_init_domain_3d - handles internal initialization
 *   get_ph_rate_point           - computed photoionzation rate from
 *                                    a point source
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"
#include "chealpix.h"
#include "ionrad.h"

#ifdef ION_RADPOINT

/* Point source only stuff */

#define NSIDE(l)       ( 1<<(l) )  /* HealPix nside parameter */
static int rebuild_interval;       /* Number of time steps between 
				      rebuilding trees */
static int is_periodic1=-1;        /* Periodic or not */
static int is_periodic2=-1;
static int is_periodic3=-1;
static Real x1lo, x2lo, x3lo;      /* Physical lower corner of grid */
static Real x1hi, x2hi, x3hi;      /* Physical upper corner of grid */
static int ray_number;             /* Ray refinement factor */
static int min_tree_level;         /* Minimum level for ray tree */

/* Variable describing the state of the random number generator */
static Rad_Ran2_State ranstate;


/* ------------------------------------------------------------
 * Utility routines
 * ------------------------------------------------------------
 *
 * This is a somewhat catchall category of routines.
 *
 */

/* Routine to find the intersection of a ray originating inside a 
 * box with the walls of the box. The
 * routine returns the coordinates of the intersection (x1_int,
 * x2_int, x3_int), the distance to the intersection, and the face of
 * the box where the intersection takes place (0/1 = low/high x, 2/3 =
 * low/high y, 4/5 = low/high z).  To handle the pathological case
 * where the ray exactly hits a box edge or corner, we also return
 * n_int, the number of intersections (up to 3). The additional
 * intersections are returned in the 2nd and 3rd elements of face. If
 * the ray and box don't intersect, n_int is set to 0. Elements of
 * face that do not have corresponding intersections are unaltered.
 */
void find_ray_box_intersect(Real x1_0, Real x2_0, Real x3_0,
			    Real n1, Real n2, Real n3,
			    Real x1blo, Real x2blo, Real x3blo,
			    Real x1bhi, Real x2bhi, Real x3bhi,
			    Real *x1_int, Real *x2_int, Real *x3_int,
			    Real *dist_int, int *nface, int *face,
			    int nexclude, int *exclude)
{
  Real dist[6];
  int i;

  /* Distance along the ray to all 6 cube faces. Note that we check if
   * n's are 0 to avoid getting inf's or NaN's.     
   */
  dist[0] = fabs(n1) < EPS ? LARGE : (x1blo - x1_0) / n1;
  dist[1] = fabs(n1) < EPS ? LARGE : (x1bhi - x1_0) / n1;
  dist[2] = fabs(n2) < EPS ? LARGE : (x2blo - x2_0) / n2;
  dist[3] = fabs(n2) < EPS ? LARGE : (x2bhi - x2_0) / n2;
  dist[4] = fabs(n3) < EPS ? LARGE : (x3blo - x3_0) / n3;
  dist[5] = fabs(n3) < EPS ? LARGE : (x3bhi - x3_0) / n3;

  /* For faces we shouldn't consider, set distance = LARGE */
  for (i=0; i<nexclude; i++)
    dist[exclude[i]] = LARGE;

  /* Pick minimum distance */
  *dist_int = LARGE;
  *nface = 0;
  for (i=0; i<6; i++) {
    if ((*dist_int > dist[i]) && (dist[i] > 0.0)) {
      *nface = 1;
      *dist_int = dist[i];
      face[0] = i;
    } else if (fabs(*dist_int-dist[i]) < 0.0) {
      face[*nface] = i;
      (*nface)++;
    }
  }
  if (*nface == 0) {
    return;
  }

  /* Set position of intersection. Note that we use this somewhat more
   * cumbersome piece of code than might be necessary to guarantee
   * that the position in the direction of intersection gets set to
   * exactly the box wall coordinate. This is useful when calling this
   * routine multiple times to walk along a ray, since it guarantees
   * that we won't wind up hitting the same wall twice in a row.
   */
  *x1_int = x1_0 + (*dist_int)*n1;
  *x2_int = x2_0 + (*dist_int)*n2;
  *x3_int = x3_0 + (*dist_int)*n3;
  for (i=0; i<*nface; i++) {
    switch (face[i]) {
    case 0: { *x1_int = x1blo; break; }
    case 1: { *x1_int = x1bhi; break; }
    case 2: { *x2_int = x2blo; break; }
    case 3: { *x2_int = x2bhi; break; }
    case 4: { *x3_int = x3blo; break; }
    case 5: { *x3_int = x3bhi; break; }
    }
  }
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1+IMM1/NTAB_RAN2)
#define RNMX (1.0-DBL_EPSILON)

/* The routine ran2() is extracted from the Numerical Recipes in C
   (version 2) code.  I've modified it to use doubles instead of
   floats. -- T. A. Gardiner -- Aug. 12, 2003

   MRK: a special note about this random number generator: we keep
   this random number generator here, and with its own random seed and
   internal state, to make sure that the same state is maintained on
   all processors for the purpose of generating ray trees, and to
   ensure that we can save and restore that state. The reproducibility
   is desirable for debugging.
*/

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
   with Bays-Durham shuffle and added safeguards.  Returns a uniform
   random deviate between 0.0 and 1.0 (exclusive of the endpoint
   values).  Call with idum = a negative integer to initialize;
   thereafter, do not alter idum between successive deviates in a
   sequence.  RNMX should appriximate the largest floating point value
   that is less than 1. */
double ion_radtransfer_ran2(long int *idum){
  int j;
  long int k;
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    ranstate.idum2=(*idum);
    for (j=NTAB_RAN2+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB_RAN2) ranstate.iv[j] = *idum;
    }
    ranstate.iy=ranstate.iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=ranstate.idum2/IQ2;
  ranstate.idum2=IA2*(ranstate.idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (ranstate.idum2 < 0) ranstate.idum2 += IM2;
  j=(int)(ranstate.iy/NDIV);              /* Will be in the range 0...NTAB_RAN2-1 */
  ranstate.iy=ranstate.iv[j]-ranstate.idum2;                /* Here idum is shuffled, idum and idum2 */
  ranstate.iv[j] = *idum;                 /* are combined to generate output */
  if (ranstate.iy < 1) ranstate.iy += IMM1;
  if ((temp=AM*ranstate.iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NDIV
#undef RNMX


/* ------------------------------------------------------------
 * Tree routines
 * ------------------------------------------------------------
 */

/* Routine to return the level for a given index */
int ray_level(int rayptr) {
  int lev=-2;
  rayptr+=4;
  while (rayptr > 0) { lev++; rayptr >>= 2; }
  return(lev);
}

/* Routine to find the next ray that is not a child */
Ray *next_non_child_ray(Ray *this_ray) {
  Ray *sibling;
  int i;

  /* Handle case of base rays differently */
  if (this_ray->parent==NULL) {
    /* Return next ray of 12 base ones, or null if this is the 12th */
    if (this_ray->mynum==11) {
      return(NULL);
    } else {
      return(this_ray+1);
    }
  } else {
    /* Return the closest sibling of this ray to its right. If there
       is no sibling, return the next non-child ray of this ray's
       parent. */
    for (i=this_ray->mynum+1; i<4; i++) {
      sibling = this_ray->parent->child[i];
      if (sibling) break;
    }
    if (i==4) return(next_non_child_ray(this_ray->parent));
    else return(sibling);
  }
}

/* Routine to find the next ray for a given ray. */
Ray *next_ray(Ray *this_ray) {
  int i;

  /* If this ray has a child, the next is its child */
  for (i=0; i<4; i++) {
    if (this_ray->child[i]) return(this_ray->child[i]);
  }
  return(next_non_child_ray(this_ray));
}


/* Free a ray */
void free_ray(Ray *ray) {
  int n;

#ifdef MPI_PARALLEL
  for (n=0; n<ray->nproc; n++) {
    if (ray->i_list) {
      if (ray->i_list[n]) { 
	free(ray->i_list[n]); 
	ray->i_list[n]=NULL;
      }
    }
    if (ray->j_list) {
      if (ray->j_list[n]) {
	free(ray->j_list[n]);
	ray->j_list[n]=NULL;
      }
    }
    if (ray->k_list) {
      if (ray->k_list[n]) {
	free(ray->k_list[n]);
	ray->k_list[n]=NULL;
      }
    }
    if (ray->cell_length_list) {
      if (ray->cell_length_list[n]) { 
	free(ray->cell_length_list[n]); 
	ray->cell_length_list[n]=NULL;
      }
    }
  }
  if (ray->proc_list) {
    free(ray->proc_list);
    ray->proc_list=NULL;
  }
  if (ray->x1_dom) {
    free(ray->x1_dom);
    ray->x1_dom=NULL;
  }
  if (ray->x2_dom) {
    free(ray->x2_dom);
    ray->x2_dom=NULL;
  }
  if (ray->x3_dom) {
    free(ray->x3_dom);
    ray->x3_dom=NULL;
  }
  if (ray->proc_length_list) {
    free(ray->proc_length_list);
    ray->proc_length_list=NULL;
  }
  if (ray->ncell) {
    free(ray->ncell);
    ray->ncell=NULL;
  }
#endif
  if (ray->i_list) { 
    free(ray->i_list); 
    ray->i_list=NULL;
  }
  if (ray->j_list) {
    free(ray->j_list);
    ray->j_list=NULL;
  }
  if (ray->k_list) {
    free(ray->k_list);
    ray->k_list=NULL;
  }
  if (ray->cell_length_list) { 
    free(ray->cell_length_list); 
    ray->cell_length_list=NULL;
  }
  for (n=0; n<4; n++) {
    if (ray->child[n]) {
      free_ray(ray->child[n]);
      free(ray->child[n]);
      ray->child[n]=NULL;
    }
  }
}


/* Free a ray tree */
void free_tree(Ray_Tree *tree) {
  int n;

  /* Free rays associated with this tree */
  for (n=0; n<12; n++) free_ray(tree->rays+n);
  free(tree->rays);
  tree->rays=NULL;
  tree->max_level=-1;
}


/* Routine to check whether a ray terminates on this processor's
   domain. */
#ifdef MPI_PARALLEL
int check_ray_end_contained(Ray *ray, Grid *pGrid) {
  Real xplim[2][3];
  Real tol[3];

  /* Return false if this is a null ray */
  if (!ray) return(0);

  /* Figure out positions of processor domain min and max relative to
     endpoint of ray */
  cc_pos(pGrid, pGrid->is, pGrid->js, pGrid->ks, 
	 &(xplim[0][0]), &(xplim[0][1]), &(xplim[0][2]));
  cc_pos(pGrid, pGrid->ie, pGrid->je, pGrid->ke,
         &(xplim[1][0]), &(xplim[1][1]), &(xplim[1][2]));
  xplim[0][0] -= ray->x1_1 + 0.5*pGrid->dx1;
  xplim[0][1] -= ray->x2_1 + 0.5*pGrid->dx2;
  xplim[0][2] -= ray->x3_1 + 0.5*pGrid->dx3;
  xplim[1][0] -= ray->x1_1 - 0.5*pGrid->dx1;
  xplim[1][1] -= ray->x2_1 - 0.5*pGrid->dx2;
  xplim[1][2] -= ray->x3_1 - 0.5*pGrid->dx3;

  /* Compute tolerances. Use a 1 cell tolerance to be safe. */
  tol[0] = pGrid->dx1;
  tol[1] = pGrid->dx2;
  tol[2] = pGrid->dx3;

  /* Check if the processor domain includes the ray end. If so, return
     true. Otherwise, return false. */
  if ((xplim[0][0] <= tol[0]) && (xplim[1][0] >= -tol[0]) &&
      (xplim[0][1] <= tol[1]) && (xplim[1][1] >= -tol[1]) &&
      (xplim[0][2] <= tol[2]) && (xplim[1][2] >= -tol[2])) {
    return(1);
  } else {
    return(0);
  }
}
#endif

/* Routine to figure out whether a given healpix ray or any of its
   children could possibly intersect a given processor domain. This is
   very difficult to do in general, but we only need an approximate
   solution. We estimate that each healpix ray covers at most a certain
   angular diameter. This means that the ray defines a cone. We then ask
   if this cone intersects the rectangular prism that is the processor
   domain. This is a fairly easy problem to do.
*/
#ifdef MPI_PARALLEL
int check_ray_intersect(Real x1, Real x2, Real x3, Real n1, Real n2, 
			Real n3, int lev, Grid *pGrid) {
  Real maxdist;
  Real x, y, z, s, a1, a2, a3;
  Real xplim[2][3], dist, xint1, xint2, xint3, theta;
  int i, j, k, nray;

  /* Figure out positions of cube max and min relative to origin of ray */
  cc_pos(pGrid, pGrid->is, pGrid->js, pGrid->ks, 
	 &(xplim[0][0]), &(xplim[0][1]), &(xplim[0][2]));
  cc_pos(pGrid, pGrid->ie, pGrid->je, pGrid->ke,
         &(xplim[1][0]), &(xplim[1][1]), &(xplim[1][2]));
  xplim[0][0] -= x1 + 0.5*pGrid->dx1;
  xplim[0][1] -= x2 + 0.5*pGrid->dx2;
  xplim[0][2] -= x3 + 0.5*pGrid->dx3;
  xplim[1][0] -= x1 - 0.5*pGrid->dx1;
  xplim[1][1] -= x2 - 0.5*pGrid->dx2;
  xplim[1][2] -= x3 - 0.5*pGrid->dx3;

  /* First check if the ray intersects the processor domain */

  /* Get distance to the plane along which each face lies. If
     the distance is positive, check if the intersection point
     lies within the square on that plane that defines the box
     face. */
  if (fabs(n1) > EPS) {
    dist = xplim[0][0] / n1;
    if (dist >= 0) {
      xint2 = dist*n2;
      xint3 = dist*n3;
      if ((xplim[0][1] <= xint2) && (xint2 <= xplim[1][1]) &&
          (xplim[0][2] <= xint3) && (xint3 <= xplim[1][2]))
        return(1);
    }
    dist = xplim[1][0] / n1;
    if (dist >= 0) {
      xint2 = dist*n2;
      xint3 = dist*n3;
      if ((xplim[0][1] <= xint2) && (xint2 <= xplim[1][1]) &&
          (xplim[0][2] <= xint3) && (xint3 <= xplim[1][2]))
        return(1);
    }
  }
  if (fabs(n2) > EPS) {
    dist = xplim[0][1] / n2;
    if (dist >= 0) {
      xint1 = dist*n1;
      xint3 = dist*n3;
      if ((xplim[0][0] <= xint1) && (xint1 <= xplim[1][0]) &&
          (xplim[0][2] <= xint3) && (xint3 <= xplim[1][2]))
        return(1);
    }
    dist = xplim[1][1] / n2;
    if (dist >= 0) {
      xint1 = dist*n1;
      xint3 = dist*n3;
      if ((xplim[0][0] <= xint1) && (xint1 <= xplim[1][0]) &&
          (xplim[0][2] <= xint3) && (xint3 <= xplim[1][2]))
        return(1);
    }
  }
  if (fabs(n3) > EPS) {
    dist = xplim[0][2] / n3;
    if (dist >= 0) {
      xint1 = dist*n1;
      xint2 = dist*n2;
      if ((xplim[0][0] <= xint1) && (xint1 <= xplim[1][0]) &&
          (xplim[0][1] <= xint2) && (xint2 <= xplim[1][1]))
        return(1);
    }
    dist = xplim[1][2] / n3;
    if (dist >= 0) {
      xint1 = dist*n1;
      xint2 = dist*n2;
      if ((xplim[0][0] <= xint1) && (xint1 <= xplim[1][0]) &&
          (xplim[0][1] <= xint2) && (xint2 <= xplim[1][1]))
        return(1);
    }
  }

  /* Now check if any of the vertices lie within the cone. */

  /* Get max angular distance from pixel center to edge on this level. The
     coefficient in front of maxdist includes the safety factor. */
  nray = 12*(1<<(2*lev));
  maxdist = 4.0/sqrt(nray);

  /* If angular distance is > pi, return true */
  if (maxdist > PI) return(1);

  /* Loop over vertices */
  for (i=0; i<2; i++) {
    x = xplim[i][0];
    for (j=0; j<2; j++) {
      y = xplim[j][1];
      for (k=0; k<2; k++) {
        z = xplim[k][2];

        /* Compute distance from vertex to the closest point on the ray that
           defines the axis of the cone, and the distance along that ray to
           the point of closest approach */
	s = x*n1 + y*n2 + z*n3;		/* Distance along ray */
	a1 = x - n1*s;			/* Vector from ray to vertex */
        a2 = y - n2*s;
        a3 = z - n3*s;

        /* Handle cases where maxdist > PI/2 and < PI/2 separately */
        if (maxdist < PI/2) {

	  /* Intersection is only possible for maxdist < PI/2 if s > 0 */
	  if (s < 0) continue;

	  /* Compute angle from ray to vertex */
	  theta = atan2(sqrt(a1*a1+a2*a2+a3*a3), s);

	  /* Check for intersection */
	  if (theta < maxdist) return(1);

	} else {

	  /* In this case, the opening angle of the cone is > PI/2, so
	     intersection is automatic if s > 0 */
	  if (s >= 0) return(1);

	  /* If s is negative, we want the angle to be larger than
	     PI - maxdist, and to compute the angle we use -s */
          theta = atan2(sqrt(a1*a1+a2*a2+a3*a3), -s);
	  if (theta > PI-maxdist) return(1);

	}
      }
    }
  }

  /* Now check each of the edges of the processor domain. Each edge is a
     segment of a line of fixed x, y, or z. We find the point of minimum 
     angular distance to the cone's center along that line. If that point
     falls within the segment that defines the domain edge, check if, at
     that point, we are within the cone. */
  /* Check edges at constant (x,y) */
  for (i=0; i<2; i++) {
    x = xplim[i][0];
    for (j=0; j<2; j++) {
      y = xplim[j][1];

      /* For this line, find the point of minimum angular distance */
      z = n3*(x*x+y*y)/(n1*x+n2*y);

      /* Does it lie in this processor domain edge? If not, continue. */
      if ((z < xplim[0][2]) || (z > xplim[1][2])) continue;

      /* Check angular distance to this point */
      s = x*n1 + y*n2 + z*n3;         /* Distance along ray */
      a1 = x - n1*s;                  /* Vector from ray to vertex */
      a2 = y - n2*s;
      a3 = z - n3*s;

      /* Handle cases where maxdist > PI/2 and < PI/2 separately */
      if (maxdist < PI/2) {

        /* Intersection is only possible for maxdist < PI/2 if s > 0 */
        if (s < 0) continue;

        /* Compute angle from ray to vertex */
        theta = atan2(sqrt(a1*a1+a2*a2+a3*a3), s);

        /* Check for intersection */
        if (theta < maxdist) return(1);

      } else {

        /* In this case, the opening angle of the cone is > PI/2, so
           intersection is automatic if s > 0 */
        if (s >= 0) return(1);

        /* If s is negative, we want the angle to be larger than
           PI - maxdist, and to compute the angle we use -s */
        theta = atan2(sqrt(a1*a1+a2*a2+a3*a3), -s);
        if (theta > PI-maxdist) return(1);

      }
    }
  }

  /* Repeat this exercise for edges with constant (x,z) and (y,z) */
  for (i=0; i<2; i++) {
    x = xplim[i][0];
    for (k=0; k<2; k++) {
      z = xplim[k][2];
      y = n2*(x*x+z*z)/(n1*x+n3*z);
      if ((y < xplim[0][1]) || (y > xplim[1][1])) continue;
      s = x*n1 + y*n2 + z*n3;         /* Distance along ray */
      a1 = x - n1*s;                  /* Vector from ray to vertex */
      a2 = y - n2*s;
      a3 = z - n3*s;
      if (maxdist < PI/2) {
        if (s < 0) continue;
        theta = atan(sqrt(a1*a1+a2*a2+a3*a3)/s);
        if (theta < maxdist) return(1);
      } else {
        if (s >= 0) return(1);
        theta = atan(sqrt(a1*a1+a2*a2+a3*a3)/(-s));
        if (theta > PI-maxdist) return(1);
      }
    }
  }
  for (j=0; j<2; j++) {
    y = xplim[j][1];
    for (k=0; k<2; k++) {
      z = xplim[k][2];
      x = n1*(y*y+z*z)/(n2*y+n3*z);
      if ((x < xplim[0][0]) || (x > xplim[1][0])) continue;
      s = x*n1 + y*n2 + z*n3;         /* Distance along ray */
      a1 = x - n1*s;                  /* Vector from ray to vertex */
      a2 = y - n2*s;
      a3 = z - n3*s;
      if (maxdist < PI/2) {
        if (s < 0) continue;
        theta = atan(sqrt(a1*a1+a2*a2+a3*a3)/s);
        if (theta < maxdist) return(1);
      } else {
        if (s >= 0) return(1);
        theta = atan(sqrt(a1*a1+a2*a2+a3*a3)/(-s));
        if (theta > PI-maxdist) return(1);
      }
    }
  }

  /* If we're here, there is no intersection */
  return(0);
}
#endif

/* Do basic initialization for a new tree */
void init_tree(Ray_Tree *tree, Real x1, Real x2, Real x3) {

  float alpha, gamma;
  Real domsize, domsize1, domsize2, domsize3;
  int i, j;

  /* Initialize rebuild counter */
  tree->rebuild_ctr = 1 % rebuild_interval;

  /* Compute rotation matrix, 
     { cos(GAMMA),             sin(GAMMA),            0.0       }
     {-cos(ALPHA)*sin(GAMMA)   cos(ALPHA)*cos(GAMMA)  sin(ALPHA)}
     { sin(ALPHA)*sin(GAMMA)  -sin(ALPHA)*cos(GAMMA)  cos(ALPHA)}
     used to rotate the tree to avoid problems coming
     from the fact that some healpix vectors always parallel the
     coordinate axes and that always having the rays point in the same
     direction can cause errors to build up over time. We pick two
     random rotation angles, using the same random seed on each
     processor (and a local random number generator that is only
     accessible to this module) to ensure that the same random angles
     are generated on each processor.
  */
  alpha = 2 * PI * ion_radtransfer_ran2(&(ranstate.seed));
  gamma = 2 * PI * ion_radtransfer_ran2(&(ranstate.seed));

  tree->rotation[0][0] = cos(gamma);
  tree->rotation[1][0] = sin(gamma);
  tree->rotation[2][0] = 0.0;
  tree->rotation[0][1] = -cos(alpha)*sin(gamma);
  tree->rotation[1][1] = cos(alpha)*cos(gamma);
  tree->rotation[2][1] = sin(alpha);
  tree->rotation[0][2] = sin(alpha)*sin(gamma);
  tree->rotation[1][2] = -sin(alpha)*cos(gamma);
  tree->rotation[2][2] = cos(alpha);

  /* Set tree center location. */
  tree->x1 = x1;
  tree->x2 = x2;
  tree->x3 = x3;

  /* Find extent of tree -- this is just the extent of the
   * computational box, except in periodic geometry, when it is the
   * position of the tree center +- half the domain size. */
  if (is_periodic1) {
    domsize1 = x1hi - x1lo;
    tree->x1blo = tree->x1 - domsize1/2;
    tree->x1bhi = tree->x1 + domsize1/2;
  } else {
    tree->x1blo = x1lo;
    tree->x1bhi = x1hi;
  }
  if (is_periodic2) {
    domsize2 = x2hi - x2lo;
    tree->x2blo = tree->x2 - domsize2/2;
    tree->x2bhi = tree->x2 + domsize2/2;
  } else {
    tree->x2blo = x2lo;
    tree->x2bhi = x2hi;
  }
  if (is_periodic3) {
    domsize3 = x3hi - x3lo;
    tree->x3blo = tree->x3 - domsize3/2;
    tree->x3bhi = tree->x3 + domsize3/2;
  } else {
    tree->x3blo = x3lo;
    tree->x3bhi = x3hi;
  }
  domsize = domsize1 < domsize2 ? domsize1 : domsize2;
  domsize = domsize < domsize3 ? domsize : domsize3;

  /* Start by setting up the 12 base rays. In parallel, these
     will be stored on every processor, unlike the rest of the ray
     tree. */
  tree->max_level = 0;
  if (!(tree->rays = 
	calloc(12, sizeof(Ray))))
    goto on_error;
  for (i=0; i<12; i++) {
    tree->rays[i].parent = NULL;
    tree->rays[i].lev = 0;
    tree->rays[i].mynum = i;
    tree->rays[i].levnum = i;
#ifdef MPI_PARALLEL
    tree->rays[i].proc_list = NULL;
    tree->rays[i].x1_dom = NULL;
    tree->rays[i].x2_dom = NULL;
    tree->rays[i].x3_dom = NULL;
    tree->rays[i].ncell = NULL;
#endif
    tree->rays[i].i_list = NULL;
    tree->rays[i].j_list = NULL;
    tree->rays[i].k_list = NULL;
    tree->rays[i].cell_length_list = NULL;
    for (j=0; j<4; j++) tree->rays[i].child[j] = NULL;
  }

  return;

 on_error:
  ath_error("[init_tree]: malloc returned a NULL pointer\n");
}

/* Same as init_tree, but rather than generating a new rotation matrix
   on the fly, we use the one we're given. This is used for restarts,
   to produce a tree with the same orientation as before the restart. */
void init_tree_rotation(Ray_Tree *tree, Real x1, Real x2, Real x3, 
			float rotation[3][3]) {

  Real domsize, domsize1, domsize2, domsize3;
  int i, j;

  /* Set rotation matrix */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++) tree->rotation[i][j] = rotation[i][j];

  /* Set tree center location. */
  tree->x1 = x1;
  tree->x2 = x2;
  tree->x3 = x3;

  /* Find extent of tree -- this is just the extent of the
   * computational box, except in periodic geometry, when it is the
   * position of the tree center +- half the domain size. */
  if (is_periodic1) {
    domsize1 = x1hi - x1lo;
    tree->x1blo = tree->x1 - domsize1/2;
    tree->x1bhi = tree->x1 + domsize1/2;
  } else {
    tree->x1blo = x1lo;
    tree->x1bhi = x1hi;
  }
  if (is_periodic2) {
    domsize2 = x2hi - x2lo;
    tree->x2blo = tree->x2 - domsize2/2;
    tree->x2bhi = tree->x2 + domsize2/2;
  } else {
    tree->x2blo = x2lo;
    tree->x2bhi = x2hi;
  }
  if (is_periodic3) {
    domsize3 = x3hi - x3lo;
    tree->x3blo = tree->x3 - domsize3/2;
    tree->x3bhi = tree->x3 + domsize3/2;
  } else {
    tree->x3blo = x3lo;
    tree->x3bhi = x3hi;
  }
  domsize = domsize1 < domsize2 ? domsize1 : domsize2;
  domsize = domsize < domsize3 ? domsize : domsize3;

  /* Start by setting up the 12 base rays. In parallel, these
     will be stored on every processor, unlike the rest of the ray
     tree. */
  tree->max_level = 0;
  if (!(tree->rays = 
	calloc(12, sizeof(Ray))))
    goto on_error;
  for (i=0; i<12; i++) {
    tree->rays[i].parent = NULL;
    tree->rays[i].lev = 0;
    tree->rays[i].mynum = i;
    tree->rays[i].levnum = i;
#ifdef MPI_PARALLEL
    tree->rays[i].proc_list = NULL;
    tree->rays[i].x1_dom = NULL;
    tree->rays[i].x2_dom = NULL;
    tree->rays[i].x3_dom = NULL;
    tree->rays[i].ncell = NULL;
#endif
    tree->rays[i].i_list = NULL;
    tree->rays[i].j_list = NULL;
    tree->rays[i].k_list = NULL;
    tree->rays[i].cell_length_list = NULL;
    for (j=0; j<4; j++) tree->rays[i].child[j] = NULL;
  }

  return;

 on_error:
  ath_error("[init_tree]: malloc returned a NULL pointer\n");
}

/* Find all the rays for tree n */
void fill_tree_ray_list(Ray_Tree *tree, Grid *pGrid, int new_max_lev) {
  Real vec[3];
  Real ray_length;
  Real x1_int, x2_int, x3_int, dist_int;
  int nface, face[3];
  int nexclude, exclude[3];
  int cur_max_lev;
  int nray;
  int i, j;
  Ray *this_ray;
#ifdef MPI_PARALLEL
  int max_level_glob, err;
#endif

  /* Store current max level */
  cur_max_lev = tree->max_level;

  /* Make sure we have something to do. If not, return. */
  if (cur_max_lev >= new_max_lev) return;

  /* Walk the ray tree */
  this_ray = &(tree->rays[0]);
  while (this_ray) {

    /* Is this ray below our current max level? If so, go to next
       ray and continue traversing tree */
    if (this_ray->lev < cur_max_lev) {
      this_ray = next_ray(this_ray);
      continue;
    }

    /* Is this ray at our current max level? */
    if (this_ray->lev == cur_max_lev) {

      /* Is this ray a leaf? If so, it has no children, so go to the
	 next non-child ray. */
      if (this_ray->leaf == 1) {
	this_ray=next_non_child_ray(this_ray);
	continue;
      }

      /* If we're here, this ray is at our current max level, and is
	 not a leaf. Initialize its children, then go to next ray. */
      for (i=0; i<4; i++) {
	if (!(this_ray->child[i] = malloc(sizeof(Ray))))
	  goto on_error;
	this_ray->child[i]->parent = this_ray;
	this_ray->child[i]->lev = this_ray->lev + 1;
	this_ray->child[i]->mynum = i;
	this_ray->child[i]->levnum = 4*this_ray->levnum + i;
#ifdef MPI_PARALLEL
	this_ray->child[i]->x1_dom = NULL;
	this_ray->child[i]->x2_dom = NULL;
	this_ray->child[i]->x3_dom = NULL;
	this_ray->child[i]->proc_list = NULL;
	this_ray->child[i]->proc_length_list = NULL;
	this_ray->child[i]->ncell = NULL;
#endif
	this_ray->child[i]->i_list = NULL;
	this_ray->child[i]->j_list = NULL;
	this_ray->child[i]->k_list = NULL;
	this_ray->child[i]->cell_length_list = NULL;
	for (j=0; j<4; j++) this_ray->child[i]->child[j] = NULL;
      }
      if (tree->max_level < this_ray->lev+1) 
	tree->max_level = this_ray->lev + 1;
      this_ray = next_ray(this_ray);
      continue;
    }

    /* If we're here, we're on a newly created ray at greater than the
       current max level, so compute the ray's properties. */

    /* Initialize this ray */
#ifdef MPI_PARALLEL
    this_ray->nproc=-1;
    this_ray->ncell=NULL;
    this_ray->proc_list=NULL;
    this_ray->x1_dom=NULL;
    this_ray->x2_dom=NULL;
    this_ray->x3_dom=NULL;
    this_ray->proc_length_list=NULL;
#else
    this_ray->ncell=-1;
#endif
    this_ray->i_list=NULL;
    this_ray->j_list=NULL;
    this_ray->k_list=NULL;
    this_ray->cell_length_list=NULL;
    for (i=0; i<4; i++) this_ray->child[i]=NULL;

    /* Figure out number of rays on current level */
    nray = 12*(1<<(2*this_ray->lev));

    /* Find the start position and direction for this ray, and its
     * maximum length. For lengths, set length of rays on levels below
     * MINLEVEL to small values to guarantee that we always trace out
     * enough rays to avoid introducing artificial asymmetries. Also
     * rotate all coordinates by the rotation matrix.
     */
    pix2vec_nest(NSIDE(this_ray->lev), this_ray->levnum, vec);
    this_ray->n1 = tree->rotation[0][0]*vec[0] +
      tree->rotation[0][1]*vec[1] +
      tree->rotation[0][2]*vec[2];
    this_ray->n2 = tree->rotation[1][0]*vec[0] +
      tree->rotation[1][1]*vec[1] +
      tree->rotation[1][2]*vec[2];
    this_ray->n3 = tree->rotation[2][0]*vec[0] +
      tree->rotation[2][1]*vec[1] +
      tree->rotation[2][2]*vec[2];

    /* Set up ray properties */
    if (this_ray->lev >= min_tree_level) 
      ray_length = sqrt(nray*min_area/(4*PI*ray_number));
    else
      ray_length = 0.0;
    if (this_ray->lev != 0) {
      this_ray->x1_0 = tree->x1 + 
	this_ray->parent->cum_length * this_ray->n1;
      this_ray->x2_0 = tree->x2 + 
	this_ray->parent->cum_length * this_ray->n2;
      this_ray->x3_0 = tree->x3 + 
	this_ray->parent->cum_length * this_ray->n3;
      ray_length -= this_ray->parent->cum_length;
    } else {
      this_ray->x1_0 = tree->x1;
      this_ray->x2_0 = tree->x2;
      this_ray->x3_0 = tree->x3;
    }

    /* Is start of new ray outside of box? If so, set this ray to have
     * length -1 and be a leaf, and continue to next ray.
     */
    if ((this_ray->x1_0 < tree->x1blo) ||
	(this_ray->x2_0 < tree->x2blo) ||
	(this_ray->x3_0 < tree->x3blo) ||
	(this_ray->x1_0 > tree->x1bhi) ||
	(this_ray->x2_0 > tree->x2bhi) ||
	(this_ray->x3_0 > tree->x3bhi)) {
      this_ray->x1_1 = this_ray->x1_0;
      this_ray->x2_1 = this_ray->x2_0;
      this_ray->x3_1 = this_ray->x3_0;
      this_ray->length = -1.0;
      this_ray->leaf = 1;
      for (i=0; i<4; i++) this_ray->child[i] = NULL;
      if (this_ray->lev < new_max_lev) this_ray = next_ray(this_ray);
      else this_ray = next_non_child_ray(this_ray);
      continue;
    }

#ifdef MPI_PARALLEL
    /* If we're working in parallel, then there is no need to populate
       the ray tree with rays in directions such that they are their       
       children cannot possibly intersect this processor. Check if any
       intersection is possible, and if not, go to next ray. However,
       there is one exception to this: it is possible that a ray may
       start completely off this processor, but its parent may be on
       this processor, since the end of a parent and the start of a
       child are not spatially coincident. In that case, the child
       needs to exist on this processor too, so the processor knows to
       send the flux to the child during the photoionization
       update. We therefore check this possibility before checking for
       an intersection.
    */
    if (!check_ray_end_contained(this_ray->parent, pGrid)) {
      if (!check_ray_intersect(this_ray->x1_0, this_ray->x2_0, this_ray->x3_0,
			       this_ray->n1, this_ray->n2, this_ray->n3, 
			       this_ray->lev, pGrid)) {
	this_ray->x1_1 = this_ray->x1_0;
	this_ray->x2_1 = this_ray->x2_0;
	this_ray->x3_1 = this_ray->x3_0;
	this_ray->length = -1.0;
	this_ray->leaf = 1;
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
	if (this_ray->lev < new_max_lev) this_ray = next_ray(this_ray);
	else this_ray = next_non_child_ray(this_ray);
	continue;
      }
    }
#endif

    /* Initialize the list of faces on which not to check for
       intersections with the ray. We do this to prevent finding the
       same intersection twice in a row when we start on a box face,
       but since we are generally going to start in the middle of a
       cell, for now exclude nothing. */
    nexclude = 0;

    /* Is the start of the ray exactly on a domain wall? If so, and
     *  the ray is headed out of the domain, skip. If the ray is
     *  headed into the domain, exclude the face the ray starts on
     *  from consideration as an intersecting face.
     */
    if (fabs(this_ray->x1_0-tree->x1blo) < 
	EPS*(tree->x1bhi-tree->x1blo)) {
      if (this_ray->n1 < 0) {
	this_ray->x1_1 = this_ray->x1_0;
	this_ray->x2_1 = this_ray->x2_0;
	this_ray->x3_1 = this_ray->x3_0;
	this_ray->length = 0.0;
	this_ray->leaf = 1;
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
	this_ray = next_ray(this_ray);
	continue;
      } else {
	exclude[nexclude++] = 0;
      }
    }
    if (fabs(this_ray->x1_0-tree->x1bhi) < 
	EPS*(tree->x1bhi-tree->x1blo)) {
      if (this_ray->n1 > 0) {
	this_ray->x1_1 = this_ray->x1_0;
	this_ray->x2_1 = this_ray->x2_0;
	this_ray->x3_1 = this_ray->x3_0;
	this_ray->length = 0.0;
	this_ray->leaf = 1;
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
	this_ray = next_ray(this_ray);
	continue;
      } else {
	exclude[nexclude++] = 1;
      }
    }
    if (fabs(this_ray->x2_0-tree->x2blo) < 
	EPS*(tree->x2bhi-tree->x2blo)) {
      if (this_ray->n2 < 0) {
	this_ray->x1_1 = this_ray->x1_0;
	this_ray->x2_1 = this_ray->x2_0;
	this_ray->x3_1 = this_ray->x3_0;
	this_ray->length = 0.0;
	this_ray->leaf = 1;
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
	this_ray = next_ray(this_ray);
	continue;
      } else {
	exclude[nexclude++] = 2;
      }
    }
    if (fabs(this_ray->x2_0-tree->x2bhi) < 
	EPS*(tree->x2bhi-tree->x2blo)) {
      if (this_ray->n2 > 0) {
	this_ray->x1_1 = this_ray->x1_0;
	this_ray->x2_1 = this_ray->x2_0;
	this_ray->x3_1 = this_ray->x3_0;
	this_ray->length = 0.0;
	this_ray->leaf = 1;
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
	this_ray = next_ray(this_ray);
	continue;
      } else {
	exclude[nexclude++] = 3;
      }
    }
    if (fabs(this_ray->x3_0-tree->x3blo) < 
	EPS*(tree->x3bhi-tree->x3blo)) {
      if (this_ray->n3 < 0) {
	this_ray->x1_1 = this_ray->x1_0;
	this_ray->x2_1 = this_ray->x2_0;
	this_ray->x3_1 = this_ray->x3_0;
	this_ray->length = 0.0;
	this_ray->leaf = 1;
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
	this_ray = next_ray(this_ray);
	continue;
      } else {
	exclude[nexclude++] = 4;
      }
    }
    if (fabs(this_ray->x3_0-tree->x3bhi) < 
	EPS*(tree->x3bhi-tree->x3blo)) {
      if (this_ray->n3 > 0) {
	this_ray->x1_1 = this_ray->x1_0;
	this_ray->x2_1 = this_ray->x2_0;
	this_ray->x3_1 = this_ray->x3_0;
	this_ray->length = 0.0;
	this_ray->leaf = 1;
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
	this_ray = next_ray(this_ray);
	continue;
      } else {
	exclude[nexclude++] = 5;
      }
    }

    /* Find the intersection of this ray with the edge of the
     * computational domain */
    find_ray_box_intersect(this_ray->x1_0,
			   this_ray->x2_0,
			   this_ray->x3_0,
			   this_ray->n1,
			   this_ray->n2,
			   this_ray->n3,
			   tree->x1blo, tree->x2blo, tree->x3blo,
			   tree->x1bhi, tree->x2bhi, tree->x3bhi,
			   &x1_int, &x2_int, &x3_int, &dist_int, 
			   &nface, face, nexclude, exclude);
    if (nface==0) {
      ath_error("[fill_tree_ray_list]: invalid box-ray intersection\n");
    }     

    /* Is ray length > length to edge of domain? */
    if (ray_length > dist_int) {

      /* Yes, so set end position and call this ray a leaf */
      this_ray->x1_1 = x1_int;
      this_ray->x2_1 = x2_int;
      this_ray->x3_1 = x3_int;
      this_ray->length = dist_int;
      if (this_ray->parent) 
	this_ray->cum_length = this_ray->parent->cum_length + dist_int;
      else this_ray->cum_length = dist_int;
      this_ray->leaf = 1;
      for (i=0; i<4; i++) this_ray->child[i] = NULL;

    } else {

      /* No, this ray isn't long enough */
      this_ray->x1_1 = this_ray->x1_0 + this_ray->n1 * ray_length;
      this_ray->x2_1 = this_ray->x2_0 + this_ray->n2 * ray_length;
      this_ray->x3_1 = this_ray->x3_0 + this_ray->n3 * ray_length;
      this_ray->length = ray_length;
      if (this_ray->parent) 
	this_ray->cum_length = this_ray->parent->cum_length + ray_length;
      else this_ray->cum_length = ray_length;
      this_ray->leaf = 0;

      /* Set up this ray's children if we're going up to the level
	 they're on. Otherwise set them to NULL for now. */
      if (this_ray->lev+1 <= new_max_lev) {
	for (i=0; i<4; i++) {
	  if (!(this_ray->child[i] = malloc(sizeof(Ray))))
	    goto on_error;
	  this_ray->child[i]->parent = this_ray;
	  this_ray->child[i]->lev = this_ray->lev + 1;
	  this_ray->child[i]->mynum = i;
	  this_ray->child[i]->levnum = 4*this_ray->levnum + i;
	}
	if (tree->max_level < this_ray->lev+1) 
	  tree->max_level = this_ray->lev + 1;
      } else {
	for (i=0; i<4; i++) this_ray->child[i] = NULL;
      }
    }

    /* Continue to next ray. */
    if (this_ray->lev < new_max_lev) this_ray = next_ray(this_ray);
    else this_ray = next_non_child_ray(this_ray);
  }

  /* If we have periodic bc's, we now need to go through the tree and
   * wrap around any points that are outside the computational domain.
   */
  if (is_periodic1 || is_periodic2 || is_periodic3) {

    /* Walk tree */
    this_ray = &(tree->rays[0]);
    while (this_ray) {
      if (is_periodic1) {
	if (this_ray->x1_0 < x1lo) this_ray->x1_0 += (x1hi - x1lo);
	else if (this_ray->x1_0 > x1hi) this_ray->x1_0 -= (x1hi - x1lo);
	if (this_ray->x1_1 < x1lo) this_ray->x1_1 += (x1hi - x1lo);
	else if (this_ray->x1_1 > x1hi) this_ray->x1_1 -= (x1hi - x1lo);
      }
      if (is_periodic2) {
	if (this_ray->x2_0 < x2lo) this_ray->x2_0 += (x2hi - x2lo);
	else if (this_ray->x2_0 > x2hi) this_ray->x2_0 -= (x2hi - x2lo);
	if (this_ray->x2_1 < x2lo) this_ray->x2_1 += (x2hi - x2lo);
	else if (this_ray->x2_1 > x2hi) this_ray->x2_1 -= (x2hi - x2lo);
      }
      if (is_periodic3) {
	if (this_ray->x3_0 < x3lo) this_ray->x3_0 += (x3hi - x3lo);
	else if (this_ray->x3_0 > x3hi) this_ray->x3_0 -= (x3hi - x3lo);
	if (this_ray->x3_1 < x3lo) this_ray->x3_1 += (x3hi - x3lo);
	else if (this_ray->x3_1 > x3hi) this_ray->x3_1 -= (x3hi - x3lo);
      }
      if (this_ray->lev < new_max_lev) this_ray = next_ray(this_ray);
      else this_ray = next_non_child_ray(this_ray);
    }
  }

#ifdef MPI_PARALLEL
  /* Set maximum level of tree to maximum level on any processor. This
     is necessary to handle communications correctly when traversing
     the tree later. */
  err = MPI_Allreduce(&(tree->max_level), &max_level_glob, 1, MPI_INT,
		      MPI_MAX, MPI_COMM_WORLD);
  if (err) ath_error("[fill_tree_ray_list]: MPI_Allreduce error = %d\n", err);
  tree->max_level = max_level_glob;
#endif

  return;

 on_error:
  ath_error("[fill_tree_ray_list]: malloc returned a NULL pointer\n");
}


#ifdef MPI_PARALLEL
void fill_tree_proc_list(Ray_Tree *tree, Grid *pGrid, int min_lev) {

  int i, j, k, i_0, j_0, k_0;
  int faceptr;
  Ray *this_ray;
  Real x1_cur, x2_cur, x3_cur;
  Real x1plo, x2plo, x3plo, x1phi, x2phi, x3phi;
  Real x1_int, x2_int, x3_int, dist_int;
  Real ray_length2, distance2;
  int nface, face[3];
  int nexclude, exclude[3];

  /* Traverse the tree, starting from the first base ray */
  this_ray = &(tree->rays[0]);
  while (this_ray) {

    /* If this is below the minimum level we are to start at, skip
       it. */
    if (this_ray->lev < min_lev) {
      if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
      else this_ray = next_non_child_ray(this_ray);
      continue;
    }

    /* Set current position to ray base and length to ray length */
    x1_cur = this_ray->x1_0;
    x2_cur = this_ray->x2_0;
    x3_cur = this_ray->x3_0;
    ray_length2 = (this_ray->length)*(this_ray->length);

    /* Does this ray start outside the problem domain our outside this
     * processor's domain? If so, stop dealing with it and continue to
     * next ray. */
    if (this_ray->length==-1.0) {
      if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
      else this_ray = next_non_child_ray(this_ray);
      continue;
    }

    /* Initialize position and processor count */
    this_ray->nproc = 0;
    if (!(this_ray->x1_dom = malloc(sizeof(Real))))
      goto on_error;
    this_ray->x1_dom[0] = x1_cur;
    if (!(this_ray->x2_dom = malloc(sizeof(Real))))
      goto on_error;
    this_ray->x2_dom[0] = x2_cur;
    if (!(this_ray->x3_dom = malloc(sizeof(Real))))
      goto on_error;
    this_ray->x3_dom[0] = x3_cur;
    this_ray->proc_list = NULL;
    this_ray->proc_length_list = NULL;

    /* Initialize the list of excluded faces */
    nexclude = 0;

    /* Find index of cell containing this position. Handle the special
     * case where the position exactly coincides with a cell edge by
     * using the direction of the ray to pick whether to go left or
     * right. Also exclude that face from consideration as the
     * possible next intersecting face
     */
    i_0 = floor((x1_cur - pGrid->x1_0) / pGrid->dx1);
    j_0 = floor((x2_cur - pGrid->x2_0) / pGrid->dx2);
    k_0 = floor((x3_cur - pGrid->x3_0) / pGrid->dx3);
    if ((x1_cur - pGrid->x1_0) / pGrid->dx1 -
	floor((x1_cur - pGrid->x1_0) / pGrid->dx1) 
	< EPS) {
      if (this_ray->n1 < 0) {
	i_0--;
	exclude[nexclude++] = 1;
      } else {
	exclude[nexclude++] = 0;
      }
    }
    if ((x2_cur - pGrid->x2_0) / pGrid->dx2 -
	floor((x2_cur - pGrid->x2_0) / pGrid->dx2) 
	< EPS) {
      if (this_ray->n2 < 0) {
	j_0--;
	exclude[nexclude++] = 3;
      } else {
	exclude[nexclude++] = 2;
      }
    }
    if ((x3_cur - pGrid->x3_0) / pGrid->dx3 -
	floor((x3_cur - pGrid->x3_0) / pGrid->dx3) 
	< EPS) {
      if (this_ray->n3 < 0) {
	k_0--;
	exclude[nexclude++] = 5;
      } else {
	exclude[nexclude++] = 4;
      }
    }

    /* Find processor on which this index resides */
    for (k=0; pD->GridArray[k][0][0].kge<k_0; k++);
    for (j=0; pD->GridArray[0][j][0].jge<j_0; j++);
    for (i=0; pD->GridArray[0][0][i].ige<i_0; i++);

    /* Walk along ray */
    distance2 = (x1_cur-this_ray->x1_0)*(x1_cur-this_ray->x1_0) +
      (x2_cur-this_ray->x2_0)*(x2_cur-this_ray->x2_0) +
      (x3_cur-this_ray->x3_0)*(x3_cur-this_ray->x3_0);
    while (ray_length2 - distance2 > sqrt(sqrt(EPS))*ray_length2) {

      /* Increment number of processors encountered */
      this_ray->nproc++;

      /* Allocate memory to store information associate with this
	 processor */
      if (this_ray->proc_list) {
	if (!(this_ray->proc_list = 
	      realloc(this_ray->proc_list, 
		      this_ray->nproc*sizeof(int))))
	  goto on_error;
      } else {
	if (!(this_ray->proc_list = malloc(sizeof(int))))
	  goto on_error;
      }
      if (!(this_ray->x1_dom = 
	    realloc(this_ray->x1_dom, (this_ray->nproc+1)*sizeof(Real))))
	goto on_error;
      if (!(this_ray->x2_dom = 
	    realloc(this_ray->x2_dom, (this_ray->nproc+1)*sizeof(Real))))
	goto on_error;
      if (!(this_ray->x3_dom = 
	    realloc(this_ray->x3_dom, (this_ray->nproc+1)*sizeof(Real))))
	goto on_error;
      if (this_ray->proc_length_list) {
	if (!(this_ray->proc_length_list = 
	      realloc(this_ray->proc_length_list, 
		      this_ray->nproc*sizeof(Real))))
	  goto on_error;
      } else {
	if (!(this_ray->proc_length_list = malloc(sizeof(Real))))
	  goto on_error;
      }

      /* Add current processor to list */
      this_ray->proc_list[this_ray->nproc-1] = pD->GridArray[k][j][i].id;

      /* Get location of the walls of this processor's domain in
       * physical space, adjusting for ghost cells
       */
      x1plo = x1lo + pD->GridArray[k][j][i].igs*pGrid->dx1;
      x2plo = x2lo + pD->GridArray[k][j][i].jgs*pGrid->dx2;
      x3plo = x3lo + pD->GridArray[k][j][i].kgs*pGrid->dx3;
      x1phi = x1lo + (pD->GridArray[k][j][i].ige+1)*pGrid->dx1;
      x2phi = x2lo + (pD->GridArray[k][j][i].jge+1)*pGrid->dx2;
      x3phi = x3lo + (pD->GridArray[k][j][i].kge+1)*pGrid->dx3;

      /* Find the intersection of this ray with the walls of this
       * domain. 
       */
      find_ray_box_intersect(x1_cur, x2_cur, x3_cur,
			     this_ray->n1,
			     this_ray->n2,
			     this_ray->n3,
			     x1plo, x2plo, x3plo,
			     x1phi, x2phi, x3phi,
			     &x1_int, &x2_int, &x3_int, &dist_int, 
			     &nface, face, nexclude, exclude);

      if (nface==0) {
	ath_error("[fill_tree_proc_list]: invalid box-ray intersection\n");
      }

      /* Are we done? If so, record info and finish with ray. */
      distance2 = (x1_int-this_ray->x1_0)*(x1_int-this_ray->x1_0) +
	(x2_int-this_ray->x2_0)*(x2_int-this_ray->x2_0) +
	(x3_int-this_ray->x3_0)*(x3_int-this_ray->x3_0);
      if (ray_length2 - distance2 <= sqrt(sqrt(EPS))*ray_length2) {
	this_ray->x1_dom[this_ray->nproc] = this_ray->x1_1;
	this_ray->x2_dom[this_ray->nproc] = this_ray->x2_1;
	this_ray->x3_dom[this_ray->nproc] = this_ray->x3_1;
	this_ray->proc_length_list[this_ray->nproc-1] = 
	  this_ray->length -
	    sqrt((x1_cur-this_ray->x1_0)*(x1_cur-this_ray->x1_0) +
		 (x2_cur-this_ray->x2_0)*(x2_cur-this_ray->x2_0) +
		 (x3_cur-this_ray->x3_0)*(x3_cur-this_ray->x3_0));
	break; /* End loop over ray length */
      }

      /* If we're here, we got to the processor edge without getting to
	 the end of the ray. */

      /* Get index of next processor, and set list of faces not to
	 consider in looking for the next intersection. */
      nexclude = 0;
      for (faceptr=0; faceptr<nface; faceptr++) {
	switch (face[faceptr]) {
	case 0: { 
	  i--; 
	  exclude[nexclude++] = 1;
	  break;
	}
	case 1: {
	  i++;
	  exclude[nexclude++] = 0;
	  break;
	}
	case 2: {
	  j--;
	  exclude[nexclude++] = 3;
	  break;
	}
	case 3: {
	  j++;
	  exclude[nexclude++] = 2;
	  break;
	}
	case 4: {
	  k--;
	  exclude[nexclude++] = 5;
	  break;
	}
	case 5: {
	  k++;
	  exclude[nexclude++] = 4;
	  break;
	}
	}
      }

      /* Move position */
      x1_cur = x1_int;
      x2_cur = x2_int;
      x3_cur = x3_int;

      /* Record the domain entry/exit and length of ray intersecting
	 domain */
      this_ray->x1_dom[this_ray->nproc] = x1_cur;
      this_ray->x2_dom[this_ray->nproc] = x2_cur;
      this_ray->x3_dom[this_ray->nproc] = x3_cur;
      this_ray->proc_length_list[this_ray->nproc-1] = dist_int;

      /* Deal with periodic bc's -- if we've gone off the edge of the
       * grid, loop back on the other side. 
       */
      if (is_periodic1) {
	if (i < 0) {
	  i = NGrid_x1 - 1;
	  x1_cur = x1hi;
	} else if (i >= NGrid_x1) {
	  i = 0;
	  x1_cur = x1lo;
	}
      }
      if (is_periodic2) {
	if (j < 0) {
	  j = NGrid_x2 - 1;
	  x2_cur = x2hi;
	} else if (j >= NGrid_x2) {
	  j = 0;
	  x2_cur = x2lo;
	}
      }
      if (is_periodic3) {
	if (k < 0) {
	  k = NGrid_x3 - 1;
	  x3_cur = x3hi;
	} else if (k >= NGrid_x3) {
	  k = 0;
	  x3_cur = x3lo;
	}
      }
    } /* End loop over ray length */

    /* Proceed to next ray */
    if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
    else this_ray = next_non_child_ray(this_ray);
  }

  return;

 on_error:
  ath_error("[fill_tree_proc_list]: malloc returned a NULL pointer\n");
}


/* Routine to go through a tree and, on each processor, remove rays
   for which are not needed. A processor needs to know about all of
   its parent rays and about its children, but not about its
   grandchildren. Thus, after pruning, every ray should either
   intersect this processor, or should have an immediate parent who
   intersects this processor. To guarantee this, we use this
   algorithm: 
   (1) Descend to rays all of whose children are leaves.
   (2) Check if either that ray or any of its children intsect this
   processors domain. If not, remove the leaves, but not the parent.
   (3) Repeat this process until we don't remove any leaves.
*/
void prune_tree(Ray_Tree *tree, Grid *pGrid) {
  Ray *this_ray;
  int i, p, nremoved;
  int nremovetot=0;

  /* Do a loop of removals */
  nremoved = 1;
  while (nremoved > 0) {
    nremoved = 0;

    /* Walk the tree down to the leaves */
    this_ray = &(tree->rays[0]);
    while (this_ray) {

      /* Are we a base ray? Never prune these. */
      if (this_ray->lev == 0) {
	if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
	else this_ray = next_non_child_ray(this_ray);
	continue;
      }

      /* Are we a leaf? We never prune these directly, we only prune
	 their parents, so go to the next ray */
      if (this_ray->leaf) {
	if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
	else this_ray = next_non_child_ray(this_ray);
	continue;
      }	

      /* Are our children allocated yet? If not, we can't know whether
	 to prune here or not, so go on. */
      if (this_ray->lev >= tree->max_level-1) {
	this_ray = next_non_child_ray(this_ray);
	continue;
      }

      /* We're not a leaf, so check if all of our children are
	 leaves. If not, continue descending the tree. */
      for (i=0; i<4; i++) {
	if (this_ray->child[i]) {
	  if (this_ray->child[i]->leaf != 1) break;
	}
      }
      if (i != 4) {
	if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
	else this_ray = next_non_child_ray(this_ray);
	continue;
      }
	
      /* All our children are leaves, so check whether we intersect
	 this processor. If so, do nothing and continue to next
	 non-child ray. */ 
      for (p=0; p<this_ray->nproc; p++) {
	if (this_ray->proc_list[p] == pGrid->my_id) break;
      }
      if (p != this_ray->nproc) {
	if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
	else this_ray = next_non_child_ray(this_ray);
	continue;
      }

      /* We do not intersect this processor, so check all our
	 children. If any of them do not intersect this processor,
	 remove them. */
      for (i=0; i<4; i++) {
	if (this_ray->child[i]) {
	  for (p=0; p<this_ray->child[i]->nproc; p++) {
	    if (this_ray->child[i]->proc_list[p] == pGrid->my_id) break;
	  }
	  if (p == this_ray->child[i]->nproc) {
	    /* This ray does not intersect this processor, so remove
	       it */
	    free_ray(this_ray->child[i]);
	    free(this_ray->child[i]);
	    this_ray->child[i] = NULL;
	    nremoved++;
	  }
	}
      }

      /* Have all the children of this ray been removed? If so, mark
	 it as a leaf. */
      for (i=0; i<4; i++) {
	if (this_ray->child[i]) break;
      }
      if (i == 4) this_ray->leaf = 1;

      /* Continue to next ray */
      if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
      else this_ray = next_non_child_ray(this_ray);
    }

    nremovetot += nremoved;
  }
}
#endif


void fill_tree_cell_list(Ray_Tree *tree, Grid *pGrid, int min_lev)
{
  int i, j, k;
  int faceptr;
  Ray *this_ray;
  Real x1_cur, x2_cur, x3_cur;
  Real x1clo, x2clo, x3clo, x1chi, x2chi, x3chi;
  Real x1_int, x2_int, x3_int, dist_int;
  Real ray_length2, distance2;
  int nface, face[3];
  int nexclude, exclude[3];
#ifdef MPI_PARALLEL
  int p;
#endif

  /* Traverse the tree, starting from the first base ray */
  this_ray = &(tree->rays[0]);
  while (this_ray) {

    /* Is this ray below the minimum level we're starting at? If so,
       go to next ray. */
    if (this_ray->lev < min_lev) {
      if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
      else this_ray = next_non_child_ray(this_ray);
      continue;
    }

    /* Create storage space for cell lists */
#ifdef MPI_PARALLEL
    if (this_ray->nproc > 0) {
      if (!(this_ray->i_list = calloc(this_ray->nproc, sizeof(int *))))
	goto on_error;
      if (!(this_ray->j_list = calloc(this_ray->nproc, sizeof(int *))))
	goto on_error;
      if (!(this_ray->k_list = calloc(this_ray->nproc, sizeof(int *))))
	goto on_error;
      if (!(this_ray->ncell = calloc(this_ray->nproc, sizeof(int))))
	goto on_error;
      if (!(this_ray->cell_length_list = 
	    calloc(this_ray->nproc, sizeof(Real *))))
	goto on_error;
    }
#else
    if (!(this_ray->i_list = malloc(sizeof(int))))
      goto on_error;
    if (!(this_ray->j_list = malloc(sizeof(int))))
      goto on_error;
    if (!(this_ray->k_list = malloc(sizeof(int))))
      goto on_error;
    if (!(this_ray->cell_length_list = malloc(sizeof(Real))))
      goto on_error;
#endif

#ifdef MPI_PARALLEL
    /* Loop over processor list */
    for (p=0; p<this_ray->nproc; p++) {

      /* If this processor in the list does not match my ID, do
	 nothing */
      if (this_ray->proc_list[p] != pGrid->my_id) continue;

      /* Initialize pointers to cell lists */
      this_ray->i_list[p] = NULL;
      this_ray->j_list[p] = NULL;
      this_ray->k_list[p] = NULL;
      this_ray->cell_length_list[p] = NULL;

      /* Set current position to where this ray enters this processor
	 domain */
      x1_cur = this_ray->x1_dom[p];
      x2_cur = this_ray->x2_dom[p];
      x3_cur = this_ray->x3_dom[p];
#else
      x1_cur = this_ray->x1_0;
      x2_cur = this_ray->x2_0;
      x3_cur = this_ray->x3_0;
#endif

      /* Initialze to exclude no faces from consideration, since
	 generally rays start in cell interiors */
      nexclude = 0;

      /* Find the cell we start in. If we start exactly at a cell
       * edge, then decide which cell to go into based on the
       * direction of the ray.
       */
      i = floor((x1_cur - pGrid->x1_0) / pGrid->dx1) - pGrid->idisp;
      j = floor((x2_cur - pGrid->x2_0) / pGrid->dx2) - pGrid->jdisp;
      k = floor((x3_cur - pGrid->x3_0) / pGrid->dx3) - pGrid->kdisp;
      if ((x1_cur - pGrid->x1_0) / pGrid->dx1 -
	  floor((x1_cur - pGrid->x1_0) / pGrid->dx1) 
	  < EPS) {
	if (this_ray->n1 < 0) {
	  i--;
	  exclude[nexclude++] = 1;
	} else {
	  exclude[nexclude++] = 0;
	}
      }
      if ((x2_cur - pGrid->x2_0) / pGrid->dx2 -
	  floor((x2_cur - pGrid->x2_0) / pGrid->dx2) 
	  < EPS) {
	if (this_ray->n2 < 0) {
	  j--;
	  exclude[nexclude++] = 3;
	} else {
	  exclude[nexclude++] = 2;
	}
      }
      if ((x3_cur - pGrid->x3_0) / pGrid->dx3 -
	  floor((x3_cur - pGrid->x3_0) / pGrid->dx3) 
	  < EPS) {
	if (this_ray->n3 < 0) {
	  k--;
	  exclude[nexclude++] = 5;
	} else {
	  exclude[nexclude++] = 4;
	}
      }

      /* Initialize cell counter */
#ifdef MPI_PARALLEL
      this_ray->ncell[p] = 0;
#else
      this_ray->ncell = 0;
#endif

      /* Start walking the ray until we get to the end of the
       * processor domain. 
       */
#ifdef MPI_PARALLEL
      ray_length2 = this_ray->proc_length_list[p] * 
	this_ray->proc_length_list[p];
      distance2 = 
	(x1_cur-this_ray->x1_dom[p])*(x1_cur-this_ray->x1_dom[p]) +
	(x2_cur-this_ray->x2_dom[p])*(x2_cur-this_ray->x2_dom[p]) +
	(x3_cur-this_ray->x3_dom[p])*(x3_cur-this_ray->x3_dom[p]);	
#else
      ray_length2 = this_ray->length * this_ray->length;
      distance2 = 
	(x1_cur-this_ray->x1_0)*(x1_cur-this_ray->x1_0) +
	(x2_cur-this_ray->x2_0)*(x2_cur-this_ray->x2_0) +
	(x3_cur-this_ray->x3_0)*(x3_cur-this_ray->x3_0);	
#endif
      while (ray_length2-distance2 > sqrt(sqrt(EPS))*ray_length2) {

	/* Increment cell counter */
#ifdef MPI_PARALLEL
	this_ray->ncell[p]++;
#else
	this_ray->ncell++;
#endif

	/* Allocate memory to store info about this cell */
#ifdef MPI_PARALLEL
	if (this_ray->i_list[p]) {
	  if (!(this_ray->i_list[p] = 
		realloc(this_ray->i_list[p], 
			this_ray->ncell[p]*sizeof(int))))
	    goto on_error;
	  if (!(this_ray->j_list[p] = 
		realloc(this_ray->j_list[p], 
			this_ray->ncell[p]*sizeof(int))))
	    goto on_error;
	  if (!(this_ray->k_list[p] = 
		realloc(this_ray->k_list[p], 
			this_ray->ncell[p]*sizeof(int))))
	    goto on_error;
	  if (!(this_ray->cell_length_list[p] = 
		realloc(this_ray->cell_length_list[p], 
			this_ray->ncell[p]*sizeof(Real))))
	    goto on_error;
	} else {
	  if (!(this_ray->i_list[p] = malloc(sizeof(int))))
	    goto on_error;
	  if (!(this_ray->j_list[p] = malloc(sizeof(int))))
	    goto on_error;
	  if (!(this_ray->k_list[p] = malloc(sizeof(int))))
	    goto on_error;
	  if (!(this_ray->cell_length_list[p] = malloc(sizeof(Real))))
	    goto on_error;
	}
#else
	if (!(this_ray->i_list = 
	      realloc(this_ray->i_list, 
		      this_ray->ncell*sizeof(int))))
	  goto on_error;
	if (!(this_ray->j_list = 
	      realloc(this_ray->j_list, 
		      this_ray->ncell*sizeof(int))))
	  goto on_error;
	if (!(this_ray->k_list = 
	      realloc(this_ray->k_list, 
		      this_ray->ncell*sizeof(int))))
	  goto on_error;
	if (!(this_ray->cell_length_list = 
	      realloc(this_ray->cell_length_list, 
		      this_ray->ncell*sizeof(Real))))
	  goto on_error;
#endif

	/* Record cell index */
#ifdef MPI_PARALLEL
	this_ray->i_list[p][this_ray->ncell[p]-1] = i;
	this_ray->j_list[p][this_ray->ncell[p]-1] = j;
	this_ray->k_list[p][this_ray->ncell[p]-1] = k;
#else
	this_ray->i_list[this_ray->ncell-1] = i;
	this_ray->j_list[this_ray->ncell-1] = j;
	this_ray->k_list[this_ray->ncell-1] = k;
#endif

	/* Get location of the walls of the cell edge.
	 */
	cc_pos(pGrid, i, j, k, &x1clo, &x2clo, &x3clo);
	x1clo -= 0.5*pGrid->dx1;
	x2clo -= 0.5*pGrid->dx2;
	x3clo -= 0.5*pGrid->dx3;
	x1chi = x1clo+pGrid->dx1;
	x2chi = x2clo+pGrid->dx2;
	x3chi = x3clo+pGrid->dx3;

	/* Find intersection with face of this cell */
	find_ray_box_intersect(x1_cur, x2_cur, x3_cur,
			       this_ray->n1,
			       this_ray->n2,
			       this_ray->n3,
			       x1clo, x2clo, x3clo,
			       x1chi, x2chi, x3chi,
			       &x1_int, &x2_int, &x3_int, &dist_int, 
			       &nface, face, nexclude, exclude);
	if (nface==0) {
	  ath_error("[fill_tree_cell_list]: invalid box-ray intersection\n");
	}

	/* Are we done? If so, record info, allocate memory for later
	 * use (now that we know how many cells there are) and finish
	 * with ray. 
	 */
#ifdef MPI_PARALLEL
	distance2 = 
	  (x1_int-this_ray->x1_dom[p])*(x1_int-this_ray->x1_dom[p]) +
	  (x2_int-this_ray->x2_dom[p])*(x2_int-this_ray->x2_dom[p]) +
	  (x3_int-this_ray->x3_dom[p])*(x3_int-this_ray->x3_dom[p]);	
	if (ray_length2 - distance2 <= sqrt(sqrt(EPS))*ray_length2) {
	  this_ray->cell_length_list[p][this_ray->ncell[p]-1] = 
	    this_ray->proc_length_list[p] -
	    sqrt((x1_cur-this_ray->x1_dom[p])*(x1_cur-this_ray->x1_dom[p]) +
		 (x2_cur-this_ray->x2_dom[p])*(x2_cur-this_ray->x2_dom[p]) +
		 (x3_cur-this_ray->x3_dom[p])*(x3_cur-this_ray->x3_dom[p]));
	  break;
	}
#else
	distance2 = 
	  (x1_int-this_ray->x1_0)*(x1_int-this_ray->x1_0) +
	  (x2_int-this_ray->x2_0)*(x2_int-this_ray->x2_0) +
	  (x3_int-this_ray->x3_0)*(x3_int-this_ray->x3_0);	
	if (ray_length2 - distance2 <= sqrt(EPS)*ray_length2) {
	  this_ray->cell_length_list[this_ray->ncell-1] = 
	    this_ray->length -
	    sqrt((x1_cur-this_ray->x1_0)*(x1_cur-this_ray->x1_0) +
		 (x2_cur-this_ray->x2_0)*(x2_cur-this_ray->x2_0) +
		 (x3_cur-this_ray->x3_0)*(x3_cur-this_ray->x3_0));
	  break;
	}
#endif

	/* Get index of next cell, and exclude the face we enter that
	   cell on from consideration as the next intersection point. */
	nexclude = 0;
	for (faceptr=0; faceptr<nface; faceptr++) {
	  switch (face[faceptr]) {
	  case 0: { 
	    i--; 
	    exclude[nexclude++] = 1;
	    break;
	  }
	  case 1: {
	    i++;
	    exclude[nexclude++] = 0;
	    break;
	  }
	  case 2: {
	    j--;
	    exclude[nexclude++] = 3;
	    break;
	  }
	  case 3: {
	    j++;
	    exclude[nexclude++] = 2;
	    break;
	  }
	  case 4: {
	    k--;
	    exclude[nexclude++] = 5;
	    break;
	  }
	  case 5: {
	    k++;
	    exclude[nexclude++] = 4;
	    break;
	  }
	  }
	}

	/* Special check: have i, j, or k wandered outside of the
	   allowed range? If we had infinite precision this would
	   never happen, because of ray length checking. In practice
	   roundoff errors can mean that we keep going for one more
	   cell when we should stop, which causes seg faults. Here we
	   check to prevent this. */
	if ((i<pGrid->is) || (i > pGrid->ie) ||
	    (j<pGrid->js) || (j > pGrid->je) ||
	    (k<pGrid->ks) || (k > pGrid->ke)) {
	  this_ray->cell_length_list[p][this_ray->ncell[p]-1] = 
	    this_ray->proc_length_list[p] -
	    sqrt((x1_cur-this_ray->x1_dom[p])*(x1_cur-this_ray->x1_dom[p]) +
		 (x2_cur-this_ray->x2_dom[p])*(x2_cur-this_ray->x2_dom[p]) +
		 (x3_cur-this_ray->x3_dom[p])*(x3_cur-this_ray->x3_dom[p]));
	  break;
	}

	/* Move position */
	x1_cur = x1_int;
	x2_cur = x2_int;
	x3_cur = x3_int;

	/* Record distance through this cell */
#ifdef MPI_PARALLEL
	this_ray->cell_length_list[p][this_ray->ncell[p]-1] = dist_int;
#else
	this_ray->cell_length_list[this_ray->ncell-1] = dist_int;
#endif
      }

#ifdef MPI_PARALLEL
    } /* End loop over processors */
#endif

    /* Proceed to next ray */
    if (this_ray->lev < tree->max_level) this_ray = next_ray(this_ray);
    else this_ray = next_non_child_ray(this_ray);
  }

  return;
  
 on_error:
  ath_error("[fill_tree_cell_list]: malloc returned a NULL pointer\n");
}


void dump_tree(Ray_Tree *tree, Grid *pGrid) {

  char fname[80];
  FILE *fp;
  int lev, n, m;
  int dummy;
  Ray *this_ray;

  for (lev=0; lev<=tree->max_level; lev++) {
    sprintf(fname, "tree.p%d.l%d.bin", pGrid->my_id, lev);
    fp=fopen(fname, "w");
    this_ray = &(tree->rays[0]);
    while (this_ray) {
      if (this_ray->lev == lev) {
	fwrite(&this_ray->mynum, sizeof(int), 1, fp);
	fwrite(&this_ray->levnum, sizeof(int), 1, fp);
	fwrite(&this_ray->n1, sizeof(Real), 1, fp);
	fwrite(&this_ray->n2, sizeof(Real), 1, fp);
	fwrite(&this_ray->n3, sizeof(Real), 1, fp);
	fwrite(&this_ray->x1_0, sizeof(Real), 1, fp);
	fwrite(&this_ray->x2_0, sizeof(Real), 1, fp);
	fwrite(&this_ray->x3_0, sizeof(Real), 1, fp);
	fwrite(&this_ray->x1_1, sizeof(Real), 1, fp);
	fwrite(&this_ray->x2_1, sizeof(Real), 1, fp);
	fwrite(&this_ray->x3_1, sizeof(Real), 1, fp);
	fwrite(&this_ray->nproc, sizeof(int), 1, fp);
	for (n=0; n<this_ray->nproc; n++) {
	  fwrite(&(this_ray->proc_list[n]), sizeof(int), 1, fp);
	  fwrite(&(this_ray->x1_dom[n]), sizeof(Real), 1, fp);
	  fwrite(&(this_ray->x2_dom[n]), sizeof(Real), 1, fp);
	  fwrite(&(this_ray->x3_dom[n]), sizeof(Real), 1, fp);
	  fwrite(&(this_ray->ncell[n]), sizeof(int), 1, fp);	  
	  for (m=0; m<this_ray->ncell[n]; m++) {
	    dummy=this_ray->i_list[n][m]+pGrid->idisp;
	    fwrite(&dummy, sizeof(int), 1, fp);	  
	    dummy=this_ray->j_list[n][m]+pGrid->jdisp;
	    fwrite(&dummy, sizeof(int), 1, fp);	  
	    dummy=this_ray->k_list[n][m]+pGrid->kdisp;
	    fwrite(&dummy, sizeof(int), 1, fp);	  
	  }
	}
	if (n > 0) {
	  fwrite(&(this_ray->x1_dom[n]), sizeof(Real), 1, fp);
	  fwrite(&(this_ray->x2_dom[n]), sizeof(Real), 1, fp);
	  fwrite(&(this_ray->x3_dom[n]), sizeof(Real), 1, fp);
	}
	this_ray = next_non_child_ray(this_ray);
      } else this_ray = next_ray(this_ray);
    }
    fclose(fp);
  }
}

/* Build a tree */
void build_tree(Ray_Tree *tree, Real x1, Real x2, Real x3, Grid *pGrid) {

  /* Initialize the tree */
  init_tree(tree, x1, x2, x3);

  /* Assemble the rays out to the minimum level */
  fill_tree_ray_list(tree, pGrid, min_tree_level);

#ifdef MPI_PARALLEL

  /* Find the list of processors through which the rays pass */
  fill_tree_proc_list(tree, pGrid, 0);

  /* Prune the tree to eliminate parts of it we don't need on this processor */
  prune_tree(tree, pGrid);
#endif

  /* Find the list of cells through which the rays pass */
  fill_tree_cell_list(tree, pGrid, 0);
}

/* Build a tree with a specified rotation matrix */
void build_tree_rotation(Ray_Tree *tree, Real x1, Real x2, Real x3, 
			 Grid *pGrid, float rotation[3][3]) {

  /* Initialize the tree */
  init_tree_rotation(tree, x1, x2, x3, rotation);

  /* Assemble the rays out to the minimum level */
  fill_tree_ray_list(tree, pGrid, min_tree_level);

#ifdef MPI_PARALLEL

  /* Find the list of processors through which the rays pass */
  fill_tree_proc_list(tree, pGrid, 0);

  /* Prune the tree to eliminate parts of it we don't need on this processor */
  prune_tree(tree, pGrid);
#endif

  /* Find the list of cells through which the rays pass */
  fill_tree_cell_list(tree, pGrid, 0);
}

/* Expand an existing tree to a new maximum level */
void expand_tree(Ray_Tree *tree, Grid *pGrid, int new_max_lev) {
  int cur_max_lev = tree->max_level;

  /* Add rays out to new max level, then fill in the lists of
     processors and cells for them. */
  fill_tree_ray_list(tree, pGrid, new_max_lev);
#ifdef MPI_PARALLEL
  fill_tree_proc_list(tree, pGrid, cur_max_lev+1);
#endif
  fill_tree_cell_list(tree, pGrid, cur_max_lev+1);
}


/* Refresh trees at specified intervals, and build new ones when new
   radiation sources are added */
void refresh_trees(Grid *pGrid) {
  int n;
  for (n=0; n < pGrid->nradpoint; n++) {

    if (pGrid->radpointlist[n].tree.rebuild_ctr == -1) {
      /* This is a new tree, so build it */
      build_tree(&(pGrid->radpointlist[n].tree), 
		 pGrid->radpointlist[n].x1,
		 pGrid->radpointlist[n].x2, 
		 pGrid->radpointlist[n].x3, 
		 pGrid);
      continue;
    } else if (pGrid->radpointlist[n].tree.rebuild_ctr == 0) {
      /* It's time to rebuild an existing tree, so destroy it and then
	 build a new one */
      free_tree(&(pGrid->radpointlist[n].tree));
      build_tree(&(pGrid->radpointlist[n].tree), 
		 pGrid->radpointlist[n].x1,
		 pGrid->radpointlist[n].x2,
		 pGrid->radpointlist[n].x3,
		 pGrid);
      continue;
    } else {
      /* Don't rebuild, just increment the counter */
      if (rebuild_interval > 0)
	pGrid->radpointlist[n].tree.rebuild_ctr =
	  (pGrid->radpointlist[n].tree.rebuild_ctr+1) % rebuild_interval;
      else
	pGrid->radpointlist[n].tree.rebuild_ctr = 1;
    }

  } /* End of tree rebuild */
}

/* ------------------------------------------------------------
 * Photoionization routines
 * ------------------------------------------------------------
 *
 * These routines do the work of computing the photoionization rates
 * for point sources.
 *
 */

void trace_ray_segment(Ray *this_ray, int p, Real s, Real ***ph_rate,
		       Real *max_flux_frac, Grid *pGrid)
{
  Real cell_vol, dummy_r, cum_length, last_length;
  Real max_flux, flux_frac;
  Real n_H, tau, ray_area, kph;
  int l;
  long nray;

  /* Compute some useful quantities */
  nray = 12*(1<<(2*this_ray->lev));
  cell_vol = pGrid->dx1 * pGrid->dx2 * pGrid->dx3;
  dummy_r = sqrt(EPS) * pGrid->dx1;

  /* Initialize ray length */
  if (this_ray->lev != min_tree_level)
    cum_length = this_ray->parent->cum_length;
  else cum_length = dummy_r;

  /* Walk through cells on this processor */
  for (l=0; l<this_ray->ncell[p]; l++) {

    /* Find flux left on this ray. If it's less than MINFLUXFAC
       of the original, stop. */
    max_flux = s / (4*PI*cum_length*cum_length);
    flux_frac = this_ray->flux / max_flux;
    if (flux_frac < MINFLUXFRAC) {
      this_ray->flux = 0.0; /* Zero out remaining flux for safety */
      break;
    }

    /* Compute neutral density and opacity in this cell */
    n_H = pGrid->U
      [this_ray->k_list[p][l]]
      [this_ray->j_list[p][l]]
      [this_ray->i_list[p][l]]
      .s[0] / m_H;
    tau = sigma_ph * n_H * this_ray->cell_length_list[p][l];

    /* Compute area of this ray */
    ray_area = 4*PI*cum_length*cum_length / nray;

    /* Compute ionization rate coefficient */
    kph = this_ray->flux * (1.0 - exp(-tau)) * ray_area
      / (n_H * cell_vol);

    /* Add the ionization rate to the grid */
    ph_rate
      [this_ray->k_list[p][l]]
      [this_ray->j_list[p][l]]
      [this_ray->i_list[p][l]]	    
      += kph;

    if (((this_ray->i_list[p][l]+pGrid->idisp==192) ||
         (this_ray->i_list[p][l]+pGrid->idisp==63)) &&
        (this_ray->k_list[p][l]+pGrid->kdisp==128) &&
        (this_ray->j_list[p][l]+pGrid->jdisp>170)) {
      printf("ray %d %d, proc %d, cell %d %d %d\n",
	     this_ray->lev, this_ray->levnum, pGrid->my_id,
	     this_ray->i_list[p][l]+pGrid->idisp,
	     this_ray->j_list[p][l]+pGrid->jdisp,
	     this_ray->k_list[p][l]+pGrid->kdisp);
    }

    /* Step to next length */
    last_length = cum_length;
    cum_length += this_ray->cell_length_list[p][l];

    /* Reduce the flux */
    this_ray->flux *= (last_length/cum_length) * 
      (last_length/cum_length) * exp(-tau);

  } /* End walk along ray segment */

  /* Record maximum flux left along ray */
  if (flux_frac > *max_flux_frac) *max_flux_frac = flux_frac;  
}

#ifdef MPI_PARALLEL
#  define INIT_BUFSIZE 4096
void get_ph_rate_point(Real s, Ray_Tree *tree, Real ***ph_rate, 
		       Grid *pGrid)
{
  int i, p, proc, err, lev, nbuf, nray_todo, receive_done;
  Ray *this_ray;
  Real max_flux_frac, max_flux_frac_glob, dummy_r;
  MPI_Status stat;
  int mpi_buf_used, msg_size;
  static MPI_Request *send_req_buf = NULL, *recv_req_buf = NULL;
  static MPI_Status *status_buf;
  static Ray **ray_buf;
  static int *proc_buf;
  static long bufsize = 0;

  /* On first call, allocate the array for MPI requests. This is
     distinct from the MPI buffer, which is handled in the main
     ion_radtransfer routine so that it can be attached and detached at
     the end of the radiation integration and we don't interfere with
     any other module that may want to use a buffer. */
  if (bufsize==0) {
    bufsize=INIT_BUFSIZE;
    recv_req_buf = (MPI_Request *) calloc(bufsize, sizeof(MPI_Request));
    send_req_buf = (MPI_Request *) calloc(bufsize, sizeof(MPI_Request));
    ray_buf = (Ray **) calloc(bufsize, sizeof(Ray *));
    proc_buf = (int *) calloc(bufsize, sizeof(int));
    status_buf = (MPI_Status *) calloc(bufsize, sizeof(MPI_Status));
  }

  /* Initialize dummy radius */
  dummy_r = sqrt(EPS) * pGrid->dx1;

  /* Loop over levels. */
  for (lev=min_tree_level; lev<=tree->max_level; lev++) {

    /* Initialize maximum flux fraction remaining */
    max_flux_frac = 0.0;

    /* Load all the valid rays on this level for this processor into a   
       buffer */
    nbuf = 0;
    this_ray = &(tree->rays[0]);
    while (this_ray) {
      if ((this_ray->lev == lev) && (this_ray->nproc > 0)) {
	if (nbuf >= bufsize) {
	  ray_buf=(Ray **) realloc(ray_buf, 2*bufsize*sizeof(Ray *));
	  proc_buf=(int *) realloc(proc_buf, 2*bufsize*sizeof(int));
	  recv_req_buf=(MPI_Request *)
	    realloc(recv_req_buf, 2*bufsize*sizeof(MPI_Request));
	  send_req_buf=(MPI_Request *)
	    realloc(send_req_buf, 2*bufsize*sizeof(MPI_Request));
	  status_buf=(MPI_Status *)
	    realloc(status_buf, 2*bufsize*sizeof(MPI_Status));
	  bufsize *= 2;
	}
	ray_buf[nbuf] = this_ray;
	proc_buf[nbuf] = 0;
	recv_req_buf[nbuf] = MPI_REQUEST_NULL;
	send_req_buf[nbuf] = MPI_REQUEST_NULL;
	nbuf++;
	this_ray = next_non_child_ray(this_ray);
      } else {
	this_ray = next_ray(this_ray);
      }
    }

    /* Now loop over the buffered rays until they're all done */
    nray_todo = nbuf;
    mpi_buf_used = 0;
    while (nray_todo) {
      for (i=0; i<nbuf; i++) {

	/* Is this ray already done? If so, continue to next. */
	if (ray_buf[i] == NULL) continue;
	else this_ray = ray_buf[i];

	/* Loop over processors on this ray */
	for (p=proc_buf[i]; p<this_ray->nproc; p++) {

	  /* Now handle various cases */
	  /* Case 1: I own this ray segment */
	  if (this_ray->proc_list[p] == pGrid->my_id) {

	    /* Case 1A: This ray is at the source, so initialize its
	       flux, trace it, and go on to the next processor. */
	    if ((this_ray->lev == min_tree_level) && (p == 0)) {
	      this_ray->flux = s / (4.0 * PI * dummy_r*dummy_r);
	      trace_ray_segment(this_ray, p, s, ph_rate, 
				&max_flux_frac, pGrid);
	      continue;
	    }

	    /* Case 1B: This ray is not at the source, but this is the
	       first processor on the ray, and the parent ray is also
	       on this processor. In this case no message passing is
	       needed, and we can get the flux locally and trace the
	       ray. */
	    if ((p == 0) && 
		(this_ray->parent->proc_list[this_ray->parent->nproc-1] ==
		 pGrid->my_id)) {
	      this_ray->flux = this_ray->parent->flux;
	      trace_ray_segment(this_ray, p, s, ph_rate,
				&max_flux_frac, pGrid);
	      continue;
	    }

	    /* Case 1C: This ray is not at the source, and the flux on
	       it is not available on this processor. In this case,
	       make an MPI request to get the flux from the processor
	       that has it. If that request cannot be fulfilled
	       immediately, skip the ray and do it later. */
	    if (p > 0) proc = this_ray->proc_list[p-1];
	    else proc = this_ray->parent->proc_list[this_ray->parent->nproc-1];
	    if (recv_req_buf[i] == MPI_REQUEST_NULL) {
	      /* No request has been made to get this flux yet, so
		 make it. */
	      err = MPI_Irecv(&(this_ray->flux), 1, 
			      MP_RL, proc, this_ray->levnum, 
			      MPI_COMM_WORLD, recv_req_buf+i);
	      if (err) ath_error("[get_ph_rate_point]: MPI_Irecv error = %d\n",err);
	    }
	    /* Check if the requested communication has been
	       completed. */
	    err = MPI_Test(recv_req_buf+i, &receive_done, &stat);
	    if (err) ath_error("[get_ph_rate_point]: MPI_Test error = %d\n",err);
	    if (receive_done) {
	      /* Communication done, so trace this ray */
	      trace_ray_segment(this_ray, p, s, ph_rate,
				&max_flux_frac, pGrid);
	      continue;
	    } else {
	      /* Communication not done, so mark the processor we're
		 up to on this ray, then skip to next ray. */
	      proc_buf[i] = p;
	      break;
	    }
	    /* End of case 1: we're on this processor's segment */

	  } else {

	    /* Case 2: We're not on this processor's segment. In this
	       case, just check if we own the segment that comes
	       before this one. */

	    /* Case 2A: This segment the first on this ray, so check
	       if we own the last segment of this ray's parent. If we
	       do, send it. */
	    if (p == 0) {
	      if (lev == min_tree_level) continue;
	      proc = this_ray->parent->proc_list[this_ray->parent->nproc-1];
	      if (proc == pGrid->my_id) {
		/* Make sure there's room in the buffer before
		   sending, and it not allocate more space. */
		MPI_Pack_size(1, MP_RL, MPI_COMM_WORLD, &msg_size);
	        mpi_buf_used += msg_size;
		if (mpi_buf_used > mpi_bufsize) {
		  MPI_Buffer_detach(&mpi_buffer, &mpi_bufsize);
		  free(mpi_buffer);
		  mpi_bufsize = 2*mpi_buf_used;
		  mpi_buffer = malloc(mpi_bufsize);
		  MPI_Buffer_attach(mpi_buffer, mpi_bufsize);
		}
		/* Send message */
		err = MPI_Ibsend(&(this_ray->parent->flux), 1,
				 MP_RL, this_ray->proc_list[p], 
				 this_ray->levnum,
				 MPI_COMM_WORLD, send_req_buf+i);
		if (err) 
		  ath_error("[get_ph_rate_point]: MPI_Isend error = %d\n",err);
	      }
	      continue;
	      /* End of case 2A: this segment is the first on the ray */
	    } else {

	      /* Case 2B: This segment is not the first on the ray, so
		 check if we own the segment on this ray that came before
		 it. If so, send it. */
	      proc = this_ray->proc_list[p-1];
	      if (proc == pGrid->my_id) {
		/* Make sure there's room in the buffer before
		   sending, and it not allocate more space. */
		MPI_Pack_size(1, MP_RL, MPI_COMM_WORLD, &msg_size);
		mpi_buf_used += msg_size;
		if (mpi_buf_used > mpi_bufsize) {
		  MPI_Buffer_detach(&mpi_buffer, &mpi_bufsize);
		  free(mpi_buffer);
		  mpi_bufsize = 2*mpi_buf_used;
		  mpi_buffer = malloc(mpi_bufsize);
		  MPI_Buffer_attach(mpi_buffer, mpi_bufsize);
		}
		/* Send message */
		err = MPI_Ibsend(&(this_ray->flux), 1,
				 MP_RL, this_ray->proc_list[p],
				 this_ray->levnum,
				 MPI_COMM_WORLD, send_req_buf+i);
		if (err) 
		  ath_error("[get_ph_rate_point]: MPI_Isend error = %d\n",err);
	      }
	    } /* End of case 2B: this segment is not first on ray */
	  } /* End of case 2: I don't own this segment */

	} /* End of loop over processors */

	/* If we've finished all the processors on this ray, decrease
	   the number of rays left to do and mark this ray as done. */
	if (p >= this_ray->nproc) {
	  ray_buf[i]=NULL;
	  nray_todo--;
	}

      } /* End loop over rays in the ray buffer */

    } /* End loop over rays on this level */

    /* Wait until all remaining send buffers are done. This also
       flushes the buffer. */
    err = MPI_Waitall(nbuf, send_req_buf, status_buf);
    if (err) 
      ath_error("[get_ph_rate_point]: MPI_Waitall error = %d\n",err);

    /* Check what the largest flux fraction remaining on any ray
       is. If it is smaller than MINFLUXFRAC of the original, we can
       stop looping over levels and end here. */
    err = MPI_Allreduce(&max_flux_frac, &max_flux_frac_glob, 1, MP_RL,
			MPI_MAX, MPI_COMM_WORLD);
    if (err) ath_error("[get_ph_rate_point]: MPI_Allreduce error = %d\n", err);
    if (max_flux_frac_glob < MINFLUXFRAC) break;

    /* If we're here, we need to go to the next tree level. Create
       that level if it doesn't exist yet. This will automatically
       update tree->max_level. */
    if (lev == tree->max_level) 
      expand_tree(tree, pGrid, tree->max_level + 1);

  } /* End loop over levels */
}
#endif


/* --------------------------------------------------------------
 * Routines to add new radiators or restore ones from checkpoints
 * --------------------------------------------------------------
 */

void add_radpoint_3d(Grid *pGrid, Real x1, Real x2, Real x3, Real s) {
  int n;

  /* Add radiator to pgrid structure */
  pGrid->nradpoint++;
  if (pGrid->nradpoint==1) {
    if (!(pGrid->radpointlist = malloc(sizeof(Radpoint))))
      goto on_error;
  } else {
    if (!(pGrid->radpointlist = realloc(pGrid->radpointlist, 
					pGrid->nradpoint*sizeof(Radpoint))))
      goto on_error;
  }
  n = pGrid->nradpoint-1;
  pGrid->radpointlist[n].x1 = x1;
  pGrid->radpointlist[n].x2 = x2;
  pGrid->radpointlist[n].x3 = x3;
  pGrid->radpointlist[n].s = s;
  pGrid->radpointlist[n].tree.max_level = -1;
  pGrid->radpointlist[n].tree.rebuild_ctr = -1;

  return;

 on_error:
  ath_error("[add_radpoint_3d]: malloc returned a NULL pointer\n");
}

void restore_radpoint_3d(Grid *pGrid, Real x1, Real x2, Real x3, Real s, 
			 int rebuild_ctr, float rotation[3][3]) {
  int n;

  /* Add radiator to pgrid structure */
  pGrid->nradpoint++;
  if (pGrid->nradpoint==1) {
    if (!(pGrid->radpointlist = malloc(sizeof(Radpoint))))
      goto on_error;
  } else {
    if (!(pGrid->radpointlist = realloc(pGrid->radpointlist, 
					pGrid->nradpoint*sizeof(Radpoint))))
      goto on_error;
  }
  n = pGrid->nradpoint-1;
  pGrid->radpointlist[n].x1 = x1;
  pGrid->radpointlist[n].x2 = x2;
  pGrid->radpointlist[n].x3 = x3;
  pGrid->radpointlist[n].s = s;

  /* Build a tree now if the rebuild counter isn't set to -1. If it
     is -1, that means this tree hadn't been built yet at the point
     when the check file was written. */

  pGrid->radpointlist[n].tree.rebuild_ctr = rebuild_ctr;
  if (rebuild_ctr != -1) 
    build_tree_rotation(&(pGrid->radpointlist[n].tree),
			pGrid->radpointlist[n].x1,
			pGrid->radpointlist[n].x2, 
			pGrid->radpointlist[n].x3,
			pGrid,
			rotation);

  return;

 on_error:
  ath_error("[restore_radpoint_3d]: malloc returned a NULL pointer\n");
}

/* ------------------------------------------------------------
 * Initialize and store routines
 * ------------------------------------------------------------
 *
 * These are called at problem setup to allocate memory and read
 * setup parameters, or when checkpointing to dump internal data
 * needed on restart.
 *
 */

void ion_radpoint_init_domain_3d(Grid *pGrid, Domain *pDomain) {

  int bcflag;
  Real area1, area2, area3;

  /* Store parallel grid information for internal use */
#ifdef MPI_PARALLEL
  pD = pDomain;
  NGrid_x1 = par_geti("parallel","NGrid_x1");
  NGrid_x2 = par_geti("parallel","NGrid_x2");
  NGrid_x3 = par_geti("parallel","NGrid_x3");
#endif

  /* Are we in periodic geometry? */
  bcflag = par_geti("grid","ibc_x1");
  if (bcflag == 4) is_periodic1 = 1; else is_periodic1 = 0;
  bcflag = par_geti("grid","ibc_x2");
  if (bcflag == 4) is_periodic2 = 1; else is_periodic2 = 0;
  bcflag = par_geti("grid","ibc_x3");
  if (bcflag == 4) is_periodic3 = 1; else is_periodic3 = 0;

  /* Where are the physical edges of the complete grid? */
  x1lo = par_getd("grid", "x1min");
  x1hi = par_getd("grid", "x1max");
  x2lo = par_getd("grid", "x2min");
  x2hi = par_getd("grid", "x2max");
  x3lo = par_getd("grid", "x3min");
  x3hi = par_getd("grid", "x3max");

  /* Read input values needed to build trees */
  ray_number = par_getd("ionradiation", "ray_number");
  min_tree_level = par_getd("ionradiation", "min_tree_level");
  rebuild_interval = par_getd("ionradiation", "rebuild_interval");
}

Rad_Ran2_State ion_radpoint_get_ranstate() { return(ranstate); }

void ion_radpoint_set_ranstate(Rad_Ran2_State newstate) {
  ranstate = newstate;
}

void ion_radpoint_init_ranstate() {
    ranstate.seed = -2674;
    ranstate.idum2 = 123456789;
    ranstate.iy = 0;
    ion_radtransfer_ran2(&(ranstate.seed));
}

#endif /* ION_RADPOINT */
