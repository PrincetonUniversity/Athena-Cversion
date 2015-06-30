#include "../copyright.h"
/*==============================================================================
 * FILE: point_source.c
 *
 * PURPOSE: Contains functions for initializing the ray tracing transfer
 *          calculations of point source radiation.  Based on original
 *          implementation by Mark Krumholz of the algorithm of Abel &
 *          Wandelt (2002).  Reimplemented in the reivsed Athena framework
 *          by Shane Davis.
 *          Current implementation follows Krumholz and includes all
 *          functions in single file.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   init_pointsource()        - initialization of the ray tracing module
 *   void build_source_trees() - build each source tree on each grid
 *   point_source_transfer()   - transfer solver in main loop
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef POINT_SOURCE
#include "chealpix.h"

#ifdef DOUBLE_PREC
#  define EPS   DBL_EPSILON          /* Arithmetic precision */
#  define LARGE DBL_MAX              /* A big number */
#  define SMALL DBL_MIN              /* A small number */
#  define MP_RL MPI_DOUBLE
#else
#  define EPS   FLT_EPSILON
#  define LARGE FLT_MAX
#  define SMALL FLT_MIN
#  define MP_RL MPI_FLOAT
#endif
//#define FRAC_FLUX 1.0e-3
#define FRAC_FLUX 0.0

/* Static variables */
static int n_tree;               /* Number of ray trees */
static RayTreeS *trees=NULL;       /* Trees of rays */
static int is_periodic1=-1;        /* Periodic or not */
static int is_periodic2=-1;
static int is_periodic3=-1;
static Real min_area;              /* Smallest cell area possible */
static Real x1lo, x2lo, x3lo;      /* Physical lower corner of grid */
static Real x1hi, x2hi, x3hi;      /* Physical upper corner of grid */
static int ray_number;             /* Ray refinement factor */
static int min_tree_level;         /* Minimum level for ray tree */

/* Physical constants and parameters read from inputs file */
//static Real sigma_ph;              /* Photoionization cross section */


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * find_ray_box_intersect() - find intersection of ray with cartesian grid
 * next_non_child_ray() - find next non-child ray in tree
 * next_ray() - find next ray in tree
 * build_tree() - build the tree of rays
 * fill_tree_proc_list() - collect info for passing tree data between grids
 * fill_tree_cell_list() - collect info for passing tree data between cells
 * update_tree_radiation() - compute optical depth and attenuation of soruce
 *============================================================================*/

void find_ray_box_intersect(Real x1_0, Real x2_0, Real x3_0,
                            Real n1, Real n2, Real n3,
                            Real x1blo, Real x2blo, Real x3blo,
                            Real x1bhi, Real x2bhi, Real x3bhi,
                            Real *x1_int, Real *x2_int, Real *x3_int,
                            Real *dist_int, int *nface, int *face);
int next_non_child_ray(int rayptr);
int next_ray(int n, int rayptr);
void build_tree(int n, Real x1, Real x2, Real x3);
#ifdef MPI_PARALLEL
void fill_tree_proc_list(int n, MeshS *pM);
#endif
void fill_tree_cell_list(int n, GridS *pGrid);
void update_tree_radiation(int n, GridS *pGrid);

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/

/* Initialize global variables, allocate memory */
void init_point_source(MeshS *pM) {

  int nsource;
  int nl, nd;
  int bcflag, i, j, k;
  Real area1, area2, area3;
  GridS *pGrid;

/* Are we in periodic geometry? */
  if (pM->BCFlag_ix1 == 4) is_periodic1 = 1; else is_periodic1 = 0;
  if (pM->BCFlag_ix2 == 4) is_periodic2 = 1; else is_periodic2 = 0;
  if (pM->BCFlag_ix3 == 4) is_periodic3 = 1; else is_periodic3 = 0;

/* Where are the physical edges of the complete grid? */
  x1lo = pM->RootMinX[0];
  x1hi = pM->RootMaxX[0];
  x2lo = pM->RootMinX[1];
  x2hi = pM->RootMaxX[1];
  x3lo = pM->RootMinX[2];
  x3hi = pM->RootMaxX[2];

/* What's the smallest area a cell face can have? */
/* This implementation uses the cell size of the root Domain, which is
 * containted in Mesh.  An SMR based treatment would probably want to
 * use a more sophisticated treatment. */ 
  area1 = SQR(pM->dx[0]);
  area2 = SQR(pM->dx[1]);
  area3 = SQR(pM->dx[2]);
  if (area1 < area2) {
    if (area1 < area3) min_area = area1;
    else min_area = area3;
  } else {
    if (area2 < area3) min_area = area2;
    else min_area = area3;
  }

/* Read input values */
  ray_number = par_geti("pointsource", "ray_number");
  min_tree_level = par_geti("pointsource", "min_tree_level");
  //sigma_ph = par_getd("pointsource", "sigma_ph");
  nsource = par_geti("pointsource", "nsource");

/* Allocate memory for Source list on each grid */
/* These arguably should be intialized in init_grid instead.  -- SWD*/
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pGrid = pM->Domain[nl][nd].Grid;
        pGrid->nsource = nsource;
        pGrid->sourcelist = (RadiatorS*)calloc_1d_array(nsource,
                                                       sizeof(RadiatorS));
        if (pGrid->sourcelist == NULL) goto on_error;
      }
    }}

/* Allocate memory for trees */
  n_tree = nsource;
  trees = (RayTreeS*) calloc_1d_array(nsource, sizeof(RayTreeS));
  if (trees == NULL) goto on_error;

  return;

 on_error:
  ath_error("[init_point_source]: malloc returned a NULL pointer\n");

}

/* Loop over Grids and sources to create initial trees.  Each tree
 * is built on each Grid.  Added by SWD 6/1/14. */
void build_source_trees(MeshS *pM) {

  int n, nd, nl;
  int i, j, k;
  GridS *pG;

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;
/* Loop over radiatiors to create initial trees */
        for (n = 0; n < n_tree; n++) {

/* Construct the ray tree. Note that every processor does this. We're
 * creating a tree where every every ray has a specified start and
 * stop point, and the ray level is set based on the distance from the
 * radiator.
 */
          build_tree(n,pG->sourcelist[n].x1,pG->sourcelist[n].x2,
                     pG->sourcelist[n].x3);
        }
      }
    }}

 /* Initialize ray_weight to 0 in all cells*/
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        if (pM->Domain[nl][nd].Grid != NULL){
          pG = pM->Domain[nl][nd].Grid;
          for (k=pG->ks; k<=pG->ke; k++) {
            for (j=pG->js; j<=pG->je; j++) {
              for (i=pG->is; i<=pG->ie; i++) {
		for (n = 0; n < n_tree; n++) {
		  pG->ray_weight[n][k][j][i] = 0.0;          
		}}}}
        }
      }
    }

  for (n = 0; n < n_tree; n++) {
#ifdef MPI_PARALLEL
/* Find the list of processors through which the rays pass */
    fill_tree_proc_list(n, pM);
#endif

/* Find the list of cells through which the rays pass */
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        if (pM->Domain[nl][nd].Grid != NULL){
          pG = pM->Domain[nl][nd].Grid;
          fill_tree_cell_list(n, pG);
        }
      }}
  }
}

/* ------------------------------------------------------------
 * Main integration routine
 * ------------------------------------------------------------
 * This is the driver routine for the radiation integration step.
 */

void point_source_transfer(DomainS *pD)
{
  GridS *pGrid;
  int i,j,k,n,ipf;
  Real opac, dt;

  if (pD->Grid != NULL){
    pGrid = pD->Grid;
/* Initialize moments */
    for (ipf=0; ipf < pGrid->npf; ipf++) { 
      for (k=pGrid->ks; k<=pGrid->ke; k++) {
	for (j=pGrid->js; j<=pGrid->je; j++) {
	  for (i=pGrid->is; i<=pGrid->ie; i++) {
	    pGrid->Jps[ipf][k][j][i] = 0.0;
	    for (n=0; n<3; n++) {
	      pGrid->Hps[ipf][k][j][i][n] = 0.0;
	    }
	  /*          for (n=0; n<6; n++) {
            pGrid->Kps[k][j][i][n] = 0.0;
	    }*/
	  }}}}

/* Compute radiation transfer for all sources */
    for (n=0; n<n_tree; n++) {
      update_tree_radiation(n, pGrid );
    }
  

/* Compute radiation forces */
    dt = pGrid->dt;
    for (ipf=0; ipf < pGrid->npf; ipf++) { 
      for (k=pGrid->ks; k<=pGrid->ke; k++) {
	for (j=pGrid->js; j<=pGrid->je; j++) {
	  for (i=pGrid->is; i<=pGrid->ie; i++) {
	    opac = pGrid->U[k][j][i].d * PSOpacity(pGrid,ipf,i,j,k);						   
	    pGrid->U[k][j][i].M1 += dt * opac * pGrid->Hps[ipf][k][j][i][0];
	    pGrid->U[k][j][i].M2 += dt * opac * pGrid->Hps[ipf][k][j][i][1];
	    pGrid->U[k][j][i].M3 += dt * opac * pGrid->Hps[ipf][k][j][i][2];
	    //printf("%d %d %d %g %g %g\n",k,j,i,opac,Prat,pGrid->Hps[k][j][i][0]);
	  }}}}
  }
}

/* ------------------------------------------------------------
 * Utility routines
 * ------------------------------------------------------------
 *
 * This is a somewhat catchall category of routines.
 *
 */

/* Routine to find the intersection of a ray with a box. The
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
                            Real *dist_int, int *nface, int *face)
{
  Real dist[6];
  Real maxdist;
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

/* Find maximum distance that is not LARGE -- used for getting
 * numerical precision.
 */
  maxdist = 0;
  for (i=0; i<6; i++)
    if ((fabs(dist[i]) > maxdist) && (dist[i] < LARGE))
      maxdist = fabs(dist[i]);

/* Pick minimum distance */
  *dist_int = LARGE;
  *nface = 0;
  for (i=0; i<6; i++) {
    if ((*dist_int > dist[i]) && (dist[i] > sqrt(EPS)*maxdist)) {
      //if ((*dist_int > dist[i]) && (fabs(dist[i]) > sqrt(EPS)*maxdist)) {
      *nface = 1;
      *dist_int = dist[i];
      face[0] = i;
    } else if (fabs(*dist_int-dist[i]) < sqrt(EPS)*maxdist) {
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


/* ------------------------------------------------------------
 * Tree routines
 * ------------------------------------------------------------
 *
 * The trees used in these routines have a complex indexing scheme
 * whose primary virtue is that finding parents, children, siblings,
 * etc, can all be done arithmetically. As a result, we don't have
 * to do any de-references when traversing a tree, which makes it fast
 * to search through. At each level there are 12*4^l rays, so there
 * are a total of 4^(l+2) - 4 rays over the first l levels. This gives
 * rise to simple formulas for the level l and index within a level j
 * of a ray whose index within the overall ray list is i:
 *    l = floor(log_4 (i + 4)) - 1
 *    j = i - l
 * Most importantly, the index of a ray and the index of its leftmost
 * child are related by
 *    l_lc = l + 1
 *    j_lc = 4 j
 *    i_lc = 4 i + 12.
 * The indices of the leftmost child's siblings are i_lc + 1, i_lc +
 * 2, and i_lc + 3. Note that this implies that the position of a ray
 * from leftmost to rightmost within its sibling group is simply its
 * index i mod 4. To test if a ray is in the 0th to 3th position, just
 * evaluate i && 0, i && 1, i && 2, or i && 3. This only applies on
 * levels > 0, however. The index i of a ray and its parent are
 * related by
 *    l_p = l - 1
 *    j_p = j / 4
 *    i_p = i / 4 - 3.
 * Again, this is restricted to the case when l > 0.
 *
 * For simplicity we define some macros and routines to enact these
 * operations below.
 */
#define LEFTCHILD(i)   ( ((i)<<2) + 12 )
#define PARENT(i)      ( ((i)<12) ? -1 : (i)/4 - 3 )
#define SIBLING(i)     ( (i)+1 )
#define RIGHTMOST(i)   ( ((i) % 4) == 3 )
#define NRAYTOT(l)     ( (1<<(2*((l)+2))) - 4 )
#define NSIDE(l)       ( 1<<(l) )  /* HealPix nside parameter */
#define LEVELPTR(i,l)  ( ((l)==0) ? (i) : ((i) - NRAYTOT((l)-1)) )

/* Routine to return the level for a given index */
int ray_level(int rayptr) {
  int lev=-2;
  rayptr+=4;
  while (rayptr > 0) { lev++; rayptr >>= 2; }
  return(lev);
}

/* Routine to find the next ray that is not a child */
int next_non_child_ray(int rayptr) {

  while (rayptr > 11) {        /* Handle level 0 case separately */
    if (!RIGHTMOST(rayptr))    /* Is this ray a rightmost child? */
      return(SIBLING(rayptr)); /* No, so just point to its sibling */
    else
      rayptr = PARENT(rayptr); /* Ray is rightmost, so point to its
                                  parent, and continue searching. */
  }

/* If we're here, then we've gotten to level 0 (i.e. rayptr
 * < 12). Look for a sibling on that level. If it exists return
 * it. Otherwise return -1.
 */
  if (rayptr<11) return(SIBLING(rayptr));
  else return(-1);
}

/* Routine to find the next ray for a given ray. */
int next_ray(int n, int rayptr) {

/* If this ray has a child, the next is its child */
  if (trees[n].rays[rayptr].leaf == 0) return(LEFTCHILD(rayptr));
  else return(next_non_child_ray(rayptr));
}

/* Build a tree of rays */
void build_tree(int n, Real x1, Real x2, Real x3) {

  Real vec[3];
  Real ray_length;
  Real x1_int, x2_int, x3_int, dist_int;
  Real domsize, domsize1, domsize2, domsize3;
  Real x1lotree, x2lotree, x3lotree;
  Real x1hitree, x2hitree, x3hitree;
  int nface, face[3];
  int nray, rayptr, lev;
  int npf;
  Ray *this_ray, *parent;
  float rotation[3][3];
/* Rotation matrix,
   { cos(GAMMA),             sin(GAMMA),            0.0       }
   {-cos(ALPHA)*sin(GAMMA)   cos(ALPHA)*cos(GAMMA)  sin(ALPHA)}
   { sin(ALPHA)*sin(GAMMA)  -sin(ALPHA)*cos(GAMMA)  cos(ALPHA)}
   used to rotate the tree to avoid problems coming
   from the fact that some healpix vectors always parallel the
   coordinate axes */
#define ALPHA 0.4
#define GAMMA 0.6
  //#define GAMMA 0.2
  rotation[0][0] = cos(GAMMA);
  rotation[1][0] = sin(GAMMA);
  rotation[2][0] = 0.0;
  rotation[0][1] = -cos(ALPHA)*sin(GAMMA);
  rotation[1][1] = cos(ALPHA)*cos(GAMMA);
  rotation[2][1] = sin(ALPHA);
  rotation[0][2] = sin(ALPHA)*sin(GAMMA);
  rotation[1][2] = -sin(ALPHA)*cos(GAMMA);
  rotation[2][2] = cos(ALPHA);
#undef ALPHA
#undef GAMMA

/* Read nfp from input */
    npf = par_geti_def("pointsource", "nf",1);
    
/* Set tree center location. */
  trees[n].x1 = x1;
  trees[n].x2 = x2;
  trees[n].x3 = x3;

/* Initialize levels */
  trees[n].max_level = -1;
  trees[n].rays = NULL;

/* Find extent of tree -- this is just the extent of the
 * computational box, except in periodic geometry, when it is the
 * position of the tree center +- half the domain size. */
  if (is_periodic1) {
    domsize1 = x1hi - x1lo;
    x1lotree = trees[n].x1 - domsize1/2;
    x1hitree = trees[n].x1 + domsize1/2;
  } else {
    x1lotree = x1lo;
    x1hitree = x1hi;
  }
  if (is_periodic2) {
    domsize2 = x2hi - x2lo;
    x2lotree = trees[n].x2 - domsize2/2;
    x2hitree = trees[n].x2 + domsize2/2;
  } else {
    x2lotree = x2lo;
    x2hitree = x2hi;
  }
  if (is_periodic3) {
    domsize3 = x3hi - x3lo;
    x3lotree = trees[n].x3 - domsize3/2;
    x3hitree = trees[n].x3 + domsize3/2;
  } else {
    x3lotree = x3lo;
    x3hitree = x3hi;
  }
  domsize = domsize1 < domsize2 ? domsize1 : domsize2;
  domsize = domsize < domsize3 ? domsize : domsize3;

/* Start walking the ray tree. We're not finding cell lists yet,
   just figuring out where rays start, branch, and end */
  rayptr = 0;
  while (rayptr>=0) {

/* Figure out number of rays for this level
 * lev   nray
 * ---   ----
 *   0     12
 *   1     48
 *   2    192
 *   3    768  etc. */
    lev = ray_level(rayptr);
    nray = 12*(1<<(2*lev));

/* Does this ray level exist yet? If not, allocate memory for it. */
 
    if (lev > trees[n].max_level) {
      if (trees[n].rays) {
        if (!(trees[n].rays =
              realloc(trees[n].rays, NRAYTOT(lev)*sizeof(Ray))))
          goto on_error;
      } else {
        if (!(trees[n].rays =
              calloc(NRAYTOT(lev), sizeof(Ray))))
          goto on_error;
      }
      trees[n].max_level = lev;
    }

/* Set up convenient pointers to this ray and its parent */
    this_ray = &(trees[n].rays[rayptr]);
    if (lev != 0) parent = &(trees[n].rays[PARENT(rayptr)]);
    else parent = NULL;

/* Allocate memory for flux array */
    this_ray->flux = calloc(npf,sizeof(Real));
    if (this_ray->flux == NULL) goto on_error;
    //if (myID_Comm_world == 0 )printf("%d %g  ",itmp++,this_ray->flux[0]);

/* Find the start position and direction for this ray, and its
 * maximum length. For lengths, set length of rays on levels below
 * MINLEVEL to small values to guarantee that we always trace out
 * enough rays to avoid introducing artificial asymmetries. Also
 * rotate all coordinates by 45 degrees about the z axis and then
 * about the x axis, to avoid introducing asymmetries when a
 * source lies exactly along cell boundaries, so some rays would
 * move parallel to cell boundaries.
 */
    pix2vec_nest(NSIDE(lev), LEVELPTR(rayptr,lev), vec);
    this_ray->n1 = rotation[0][0]*vec[0] + rotation[0][1]*vec[1] +
      rotation[0][2]*vec[2];
    this_ray->n2 = rotation[1][0]*vec[0] + rotation[1][1]*vec[1] +
      rotation[1][2]*vec[2];
    this_ray->n3 = rotation[2][0]*vec[0] + rotation[2][1]*vec[1] +
      rotation[2][2]*vec[2];
    if (lev >= min_tree_level)
      ray_length = sqrt(nray*min_area/(4*PI*ray_number));
    else
      ray_length = 0.0;
    if (lev != 0) {
      this_ray->x1_0 = trees[n].x1 + parent->cum_length * this_ray->n1;
      this_ray->x2_0 = trees[n].x2 + parent->cum_length * this_ray->n2;
      this_ray->x3_0 = trees[n].x3 + parent->cum_length * this_ray->n3;
      ray_length -= parent->cum_length;
    } else {
      this_ray->x1_0 = trees[n].x1;
      this_ray->x2_0 = trees[n].x2;
      this_ray->x3_0 = trees[n].x3;
    }
/* Is start of new ray outside of box? If so, set this ray to have
 * length -1 and be a leaf, and continue to next ray.
 */
    if ((this_ray->x1_0 < x1lotree) ||
        (this_ray->x2_0 < x2lotree) ||
        (this_ray->x3_0 < x3lotree) ||
        (this_ray->x1_0 > x1hitree) ||
        (this_ray->x2_0 > x2hitree) ||
        (this_ray->x3_0 > x3hitree)) {
      this_ray->x1_1 = this_ray->x1_0;
      this_ray->x2_1 = this_ray->x2_0;
      this_ray->x3_1 = this_ray->x3_0;
      this_ray->length = -1.0;
      this_ray->leaf = 1;
      rayptr = next_ray(n, rayptr);
      continue;
    }

/* Is the start of the ray exactly on a domain wall? If so, and
 *  the ray is headed out of the domain, skip.
 */
    if (
        ((fabs(this_ray->x1_0-x1lotree)<EPS*(x1hitree-x1lotree)) &&
         (this_ray->n1 < 0))
        ||
        ((fabs(this_ray->x2_0-x2lotree)<EPS*(x2hitree-x2lotree)) &&
         (this_ray->n2 < 0))
        ||
        ((fabs(this_ray->x3_0-x3lotree)<EPS*(x3hitree-x3lotree)) &&
         (this_ray->n3 < 0))
        ||
        ((fabs(this_ray->x1_0-x1hitree)<EPS*(x1hitree-x1lotree)) &&
         (this_ray->n1 > 0))
        ||
        ((fabs(this_ray->x2_0-x2hitree)<EPS*(x2hitree-x2lotree)) &&
         (this_ray->n2 > 0))
        ||
        ((fabs(this_ray->x3_0-x3hitree)<EPS*(x3hitree-x3lotree)) &&
         (this_ray->n3 > 0))
        ) {
      this_ray->x1_1 = this_ray->x1_0;
      this_ray->x2_1 = this_ray->x2_0;
      this_ray->x3_1 = this_ray->x3_0;
      this_ray->length = 0.0;
      this_ray->leaf = 1;

      rayptr = next_ray(n, rayptr);
      continue;
    }

/* Find the intersection of this ray with the edge of the
 * computational domain */
    find_ray_box_intersect(this_ray->x1_0,
                           this_ray->x2_0,
                           this_ray->x3_0,
                           this_ray->n1,
                           this_ray->n2,
                           this_ray->n3,
                           x1lotree, x2lotree, x3lotree,
                           x1hitree, x2hitree, x3hitree,
                           &x1_int, &x2_int, &x3_int, &dist_int,
                           &nface, face);
    if (nface==0) {
/*
  fprintf(stderr,"[build_tree]: invalid box-ray intersection\n");
  fprintf(stderr,"   ray = %d\n", rayptr);
  fprintf(stderr,"   x_0 = (%e, %e, %e)\n",
  this_ray->x1_0,
  this_ray->x2_0,
  this_ray->x3_0);
  fprintf(stderr,"   n   = (%f, %f, %f)\n",
  this_ray->n1,
  this_ray->n2,
  this_ray->n3);
  fprintf(stderr,"   blo = (%e, %e, %e)\n",
  x1lotree, x2lotree, x3lotree);
  fprintf(stderr,"   bhi = (%e, %e, %e)\n",
  x1hitree, x2hitree, x3hitree);
  fflush(stderr);
*/
      ath_error("[build_tree]: invalid box-ray intersection\n");
    }

/* Is this ray long enough to get to edge of computational
 * domain?
 */

    if (ray_length > dist_int) {

/* Yes, so set end position and call this ray a leaf */
      this_ray->x1_1 = x1_int;
      this_ray->x2_1 = x2_int;
      this_ray->x3_1 = x3_int;
      this_ray->length = dist_int;
      if (parent) this_ray->cum_length = parent->cum_length + dist_int;
      else this_ray->cum_length = dist_int;
      this_ray->leaf = 1;

    } else {

/* No, this ray isn't long enough */
      this_ray->x1_1 = this_ray->x1_0 + this_ray->n1 * ray_length;
      this_ray->x2_1 = this_ray->x2_0 + this_ray->n2 * ray_length;
      this_ray->x3_1 = this_ray->x3_0 + this_ray->n3 * ray_length;
      this_ray->length = ray_length;
      if (parent) this_ray->cum_length = parent->cum_length + ray_length;
      else this_ray->cum_length = ray_length;
      this_ray->leaf = 0;

    }

/* Continue to next ray. */
    rayptr = next_ray(n, rayptr);

  }

/* If we have periodic bc's, we now need to go through the tree and
 * wrap around any points that are outside the computational domain.
 */
  if (is_periodic1 || is_periodic2 || is_periodic3) {

/* Walk tree */
    rayptr = 0;
    while (rayptr>=0) {
      this_ray = &(trees[n].rays[rayptr]);
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
      rayptr=next_ray(n, rayptr);
    }

  }

  return;

 on_error:
  ath_error("[build_tree]: malloc returned a NULL pointer\n");
}


#ifdef MPI_PARALLEL
/* Updated by SWD.  Assumes uniform Cartesian grid.*/
void fill_tree_proc_list(int n, MeshS *pM) {

  DomainS *pD;
  int i, j, k, i_0, j_0, k_0;
  int rayptr, faceptr;
  Ray *this_ray;
  Real x1_cur, x2_cur, x3_cur;
  Real x1plo, x2plo, x3plo, x1phi, x2phi, x3phi;
  Real x1_int, x2_int, x3_int, dist_int;
  Real ray_length;
  int nface, face[3];


/* Traverse the tree, starting from the first base ray */
  rayptr = 0;
  while (rayptr>=0) {

/* Set pointer to this ray */
    this_ray = &(trees[n].rays[rayptr]);

/* Set current position to ray base and length to ray length */
    x1_cur = this_ray->x1_0;
    x2_cur = this_ray->x2_0;
    x3_cur = this_ray->x3_0;
    ray_length = this_ray->length;

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

/* Does this ray start outside the problem domain? If so, stop
 * dealing with it and continue to next ray */
    if (ray_length==-1.0) {
      rayptr = next_ray(n, rayptr);
      continue;
    }

/* Find index of cell containing this position. Handle the special
 * case where the position exactly coincides with a cell edge by
 * using the direction of the ray to pick whether to go left or
 * right.
 */

/* SWD -- currently assumes there is a single root Domain ie no SMR */
    pD = &(pM->Domain[0][0]);
    i_0 = floor((x1_cur - pD->MinX[0]) / pD->dx[0]);
    j_0 = floor((x2_cur - pD->MinX[1]) / pD->dx[1]);
    k_0 = floor((x3_cur - pD->MinX[2]) / pD->dx[2]);
    if ((x1_cur - pD->MinX[0]) / pD->dx[0] -
        floor((x1_cur - pD->MinX[0]) / pD->dx[0])
        < EPS) {
      if (this_ray->n1 < 0) i_0--;
    }
    if ((x2_cur - pD->MinX[1]) / pD->dx[1] -
        floor((x2_cur - pD->MinX[1]) / pD->dx[1])
        < EPS) {
      if (this_ray->n2 < 0) j_0--;
    }
    if ((x3_cur - pD->MinX[2]) / pD->dx[2] -
        floor((x3_cur - pD->MinX[2]) / pD->dx[2])
        < EPS) {
      if (this_ray->n3 < 0) k_0--;
    }

/* Find processor on which this index resides */
    for (k=0; k < pD->NGrid[2]; k++)
      if (pD->GData[k][0][0].Disp[2]+pD->GData[k][0][0].Nx[2] > k_0) break;
    for (j=0; j < pD->NGrid[1]; j++)
      if (pD->GData[0][j][0].Disp[1]+pD->GData[0][j][0].Nx[1] > j_0) break;
    for (i=0; i < pD->NGrid[0]; i++)
      if (pD->GData[0][0][i].Disp[0]+pD->GData[0][0][i].Nx[0] > i_0) break;

/* Walk along ray */
    while (ray_length > 0) {

/* Increment number of processors encountered */
      this_ray->nproc++;

/* Allocate memory to store information associated with this
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
      this_ray->proc_list[this_ray->nproc-1] = pD->GData[k][j][i].ID_Comm_world;

/* Get location of the walls of this processor's domain in
 * physical space, adjusting for ghost cells
 */
      x1plo = x1lo + pD->GData[k][j][i].Disp[0]*pD->dx[0];
      x2plo = x2lo + pD->GData[k][j][i].Disp[1]*pD->dx[1];
      x3plo = x3lo + pD->GData[k][j][i].Disp[2]*pD->dx[2];
      x1phi = x1lo + (pD->GData[k][j][i].Disp[0]+
                      pD->GData[k][j][i].Nx[0])*pD->dx[0];
      x2phi = x2lo + (pD->GData[k][j][i].Disp[1]+
                      pD->GData[k][j][i].Nx[1])*pD->dx[1];
      x3phi = x3lo + (pD->GData[k][j][i].Disp[2]+
                      pD->GData[k][j][i].Nx[2])*pD->dx[2];

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
                             &nface, face);

      if (nface==0) {
        fprintf(stderr,"[build_tree]: invalid box-ray intersection\n");
        fprintf(stderr,"mx,my,mz = (%g, %g, %g)\n",pD->MinX[0],pD->MinX[1],pD->MinX[2]);
        fprintf(stderr,"dx,dy,dz = (%g, %g, %g)\n",pD->dx[0],pD->dx[1],pD->dx[2]);
        fprintf(stderr,"i,j,k = (%d, %d, %d)\n",i,j,k);
        fprintf(stderr,"i0,j0,k0 = (%d, %d, %d)\n",i_0,j_0,k_0);
        fprintf(stderr,"   ray = %d\n", rayptr);
        fprintf(stderr,"   x_0 = (%e, %e, %e)\n",
                x1_cur,
                x2_cur,
                x3_cur);
        fprintf(stderr,"   n   = (%f, %f, %f)\n",
                this_ray->n1,
                this_ray->n2,
                this_ray->n3);
        fprintf(stderr,"   blo = (%e, %e, %e)\n",
                x1plo, x2plo, x3plo);
        fprintf(stderr,"   bhi = (%e, %e, %e)\n",
                x1phi, x2phi, x3phi);
                fflush(stderr);
        ath_error("[fill_tree_proc_list]: invalid box-ray intersection\n");
      }

/* Are we done? If so, record info, allocate memory for future
 * use, and finish with ray.
 */
      if (ray_length - dist_int < sqrt(EPS)*this_ray->length) {
        this_ray->x1_dom[this_ray->nproc] = this_ray->x1_1;
        this_ray->x2_dom[this_ray->nproc] = this_ray->x2_1;
        this_ray->x3_dom[this_ray->nproc] = this_ray->x3_1;
        this_ray->proc_length_list[this_ray->nproc-1] = ray_length;
        if (!(this_ray->etau_list=calloc(this_ray->nproc, sizeof(Real **))))
          goto on_error;
        break;
      }

/* Subtract distance to domain wall from ray length */
      ray_length -= dist_int;

/* Get index of next processor */
      for (faceptr=0; faceptr<nface; faceptr++) {
        switch (face[faceptr]) {
        case 0: { i--; break; }
        case 1: { i++; break; }
        case 2: { j--; break; }
        case 3: { j++; break; }
        case 4: { k--; break; }
        case 5: { k++; break; }
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
          i = pD->NGrid[0] - 1;
          x1_cur = x1hi;
        } else if (i >= pD->NGrid[0]) {
          i = 0;
          x1_cur = x1lo;
        }
      }
      if (is_periodic2) {
        if (j < 0) {
          j = pD->NGrid[1] - 1;
          x2_cur = x2hi;
        } else if (j >= pD->NGrid[1]) {
          j = 0;
          x2_cur = x2lo;
        }
      }
      if (is_periodic3) {
        if (k < 0) {
          k = pD->NGrid[2] - 1;
          x3_cur = x3hi;
        } else if (k >= pD->NGrid[2]) {
          k = 0;
          x3_cur = x3lo;
        }
      }
    } /* while (ray_length > 0) */

/* Proceed to next ray */
    rayptr=next_ray(n, rayptr);
  }

  return;

 on_error:
  ath_error("[fill_tree_proc_list]: malloc returned a NULL pointer\n");
}
#endif


void fill_tree_cell_list(int n, GridS *pGrid)
{
  DomainS *pD;
  int i, j, k;
  int rayptr, faceptr;
  Ray *this_ray;
  Real x1_cur, x2_cur, x3_cur;
  Real x1clo, x2clo, x3clo, x1chi, x2chi, x3chi;
  Real x1_int, x2_int, x3_int, dist_int;
  Real ray_length;
  int nface, face[3];
  int lev, nray;
  Real ray_sa;

#ifdef MPI_PARALLEL
  int p;
#endif

/* Traverse the tree, starting from the first base ray */
  rayptr = 0;
  while (rayptr>=0) {

/* Level, number of rays, and ray solid angle, Used for computing
   ray_weight, which is used for weighting different ray contributions
   within cells */
    lev = ray_level(rayptr);
    nray = 12*(1<<(2*lev));
    ray_sa = 4 * PI / nray;
    //ray_sa = 1.0 / nray;

/* Set pointer to this ray */
    this_ray = &(trees[n].rays[rayptr]);

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
      if (this_ray->proc_list[p] != myID_Comm_world) continue;

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

/* Find the cell we start in. If we start exactly at a cell
 * edge, then decide which cell to go into based on the
 * direction of the ray.
 */
      i = floor((x1_cur - pGrid->MinX[0]) / pGrid->dx1) + pGrid->is;
      j = floor((x2_cur - pGrid->MinX[1]) / pGrid->dx2) + pGrid->js;
      k = floor((x3_cur - pGrid->MinX[2]) / pGrid->dx3) + pGrid->ks;
      if ((x1_cur - pGrid->MinX[0]) / pGrid->dx1 -
          floor((x1_cur - pGrid->MinX[0]) / pGrid->dx1)
          < EPS) {
        if (this_ray->n1 < 0) i--;
      }
      if ((x2_cur - pGrid->MinX[1]) / pGrid->dx2 -
          floor((x2_cur - pGrid->MinX[1]) / pGrid->dx2)
          < EPS) {
        if (this_ray->n2 < 0) j--;
      }
      if ((x3_cur - pGrid->MinX[2]) / pGrid->dx3 -
          floor((x3_cur - pGrid->MinX[2]) / pGrid->dx3)
          < EPS) {
        if (this_ray->n3 < 0) k--;
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
      ray_length = this_ray->proc_length_list[p];
#else
      ray_length = this_ray->length;
#endif

#ifdef MPI_PARALLEL
      while (ray_length > sqrt(EPS)*this_ray->proc_length_list[p]) {
#else
      while (ray_length > sqrt(EPS)*this_ray->length) {
#endif

/* Increment cell counter */
#ifdef MPI_PARALLEL
        (this_ray->ncell[p])++;
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
                               &nface, face);
        if (nface==0) {
          ath_error("[fill_tree_cell_list]: invalid box-ray intersection\n");
        }

/* Are we done? If so, record info, allocate memory for later
 * use (now that we know how many cells there are, end finish
 * with ray.
 */
#ifdef MPI_PARALLEL
        //        if (ray_length - dist_int < EPS*this_ray->proc_length_list[p]) {
        if (ray_length - dist_int < EPS*this_ray->length) {
          this_ray->cell_length_list[p][this_ray->ncell[p]-1] = ray_length;
          pGrid->ray_weight[n][k][j][i] += ray_length * ray_sa;
          break;
        }
#else
        if (ray_length - dist_int < EPS*this_ray->length) {
          this_ray->cell_length_list[this_ray->ncell-1] = ray_length;
          pGrid->ray_weight[n][k][j][i] += ray_length * ray_sa;
          break;
        }
#endif
/* Used for normalization radiation moments, which are weighted by solid angle
 * and distance through the cell */
        pGrid->ray_weight[n][k][j][i] += dist_int * ray_sa;

/* Subtract distance to cell wall from ray length */
        ray_length -= dist_int;

/* Get index of next cell */
        for (faceptr=0; faceptr<nface; faceptr++) {
          switch (face[faceptr]) {
          case 0: { i--; break; }
          case 1: { i++; break; }
          case 2: { j--; break; }
          case 3: { j++; break; }
          case 4: { k--; break; }
          case 5: { k++; break; }
          }
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

/* Allocate memory to hold list of  e^(-tau)
 * factors along the ray on this processor. */
#ifdef MPI_PARALLEL
      if (!(this_ray->etau_list[p]=(Real**)calloc_2d_array(pGrid->npf,this_ray->ncell[p],
							   sizeof(Real))))
        goto on_error;
      }
#else
      if (!(this_ray->etau_list=(Real**)calloc_2d_array(pGrid->npf,this_ray->ncell,
							sizeof(Real))))
        goto on_error;
#endif

/* Proceed to next ray */
    rayptr = next_ray(n, rayptr);
  }

  return;

 on_error:
  ath_error("[fill_tree_cell_list]: malloc returned a NULL pointer\n");
}


/* Routine to compute radiation field. The code required to do this
 * in serial versus in parallel, and the opmtimum strategy for doing
 * it, is sufficiently different that it's easiest just to write the
 * routine twice, once for each case.
 */
#ifdef MPI_PARALLEL
void update_tree_radiation(int n, GridS *pGrid)
{
  int rayptr, lev, l;
  int proc, err;
  Ray *this_ray, *parent;
  Real kph, cum_length, last_length;
  Real dummy_r, flux_frac, max_frac, max_frac_glob;
  Real density;
  MPI_Status stat;
  int p;
  int i,j,k,ipf;
  Real Jray, ray_sa;
  int nray;
  Real x1c,x2c,x3c,rcell2;
  int MPI_counter=0;

/* Dummy radius, used in initializing fluxes. Its value doesn't
 * matter, it is just used for coding convenience.
 */
  dummy_r = 1.0e-6 * pGrid->dx1;
  //printf("i: %d %p %p\n",myID_Comm_world,&max_frac,&max_frac_glob);
/* Traverse the tree. Note that we traverse it by level, so that we
 * do all the level 0 rays before any of the level 1 rays, etc. This
 * facilitates parallelism.
 */
  for (lev = min_tree_level; lev <= trees[n].max_level; lev++) {

    nray = 12*(1<<(2*lev));
    ray_sa = 4 * PI / nray;
    //ray_sa = 1.0 / nray;

/* Initialize fraction of flux left */
    
    for (rayptr = NRAYTOT(lev-1); rayptr < NRAYTOT(lev); rayptr++) {
      max_frac = 0.0;
/* Pointer to this ray */
      this_ray = &(trees[n].rays[rayptr]);

/* Loop over processor list */
      for (p=0; p<this_ray->nproc; p++) {

/* If this processor in the list does not match my ID, do
   nothing */
        if (this_ray->proc_list[p] != myID_Comm_world) continue;

/* Walk along the ray to get opacities. We do the opacities
 * separately because this part is completely parallel.
 */
        for (l=0; l<this_ray->ncell[p]; l++) {

          i = this_ray->i_list[p][l];
          j = this_ray->j_list[p][l];
          k = this_ray->k_list[p][l];

/* Compute exp of optical depth in cell */
          density = pGrid->U[k][j][i].d;
	  for(ipf=0; ipf < pGrid->npf; ipf++) {
	    this_ray->etau_list[p][ipf][l] = exp(-PSOpacity(pGrid,ipf,i,j,k)*
	      density * this_ray->cell_length_list[p][l]);
	  }
        }
      }
    } /* End loop over rays on this level */

/* Now traverse rays on this level again, updating moments */
    for (rayptr = NRAYTOT(lev-1); rayptr < NRAYTOT(lev); rayptr++) {

/* Pointer to this ray and its parent*/
      this_ray = &(trees[n].rays[rayptr]);
      if (lev != 0) parent = &(trees[n].rays[PARENT(rayptr)]);

/* Initialize ray length */
      if (lev != min_tree_level)
        cum_length = parent->cum_length;
      else
        cum_length = dummy_r;

/* Loop over processor list on this ray */
      for (p=0; p<this_ray->nproc; p++) {

/* ----------------------------- Initialize flux on ray ----------------------------- */
        if (this_ray->proc_list[p] == myID_Comm_world) {
/* This ray is on my processor */

          if ((lev == min_tree_level) && (p == 0)) {
/* We  are at source */
	    for(ipf=0; ipf < pGrid->npf; ipf++)
	      this_ray->flux[ipf] = 1.0;
	
          } else {
/* If we are the first processor on this ray, get the flux
 * from the parent. This may be on another processor. */
            if (p == 0) {
              proc = parent->proc_list[parent->nproc-1];
              if (proc == myID_Comm_world) {
		for(ipf=0; ipf < pGrid->npf; ipf++) 
		  this_ray->flux[ipf] = parent->flux[ipf];		
              } else {
		err = MPI_Recv(this_ray->flux, pGrid->npf,
                               MP_RL, proc, flux_tag,
                               MPI_COMM_WORLD, &stat);
                if (err) ath_error("[update_tree_radiation]: MPI_Recv error = %d\n",err);
              }
            } else { /* We are not the first processor, so get the
                        flux from the previous processor */
              proc = this_ray->proc_list[p-1];
              err = MPI_Recv(this_ray->flux, pGrid->npf,
                             MP_RL, proc, flux_tag,
                             MPI_COMM_WORLD, &stat);
              if (err) ath_error("[update_tree_radiation]: MPI_Recv error = %d\n",err);
            }
/* End we are not at source */
          }
        } else {
/* I may own the parent of this ray segment, in which case I
 * need to pass it to the processor that is handling this
 * segment. */	  
          if (lev > min_tree_level) {
            if (p==0){
              if (parent->proc_list[parent->nproc-1]==myID_Comm_world) {
                err = MPI_Send(parent->flux, pGrid->npf,
                               MP_RL, this_ray->proc_list[0], flux_tag,
                               MPI_COMM_WORLD);
                if (err)
                  ath_error("[update_tree_radiation]: MPI_Send error = %d\n",err);
              }
            } else if (this_ray->proc_list[p-1]==myID_Comm_world) {
              err = MPI_Send(this_ray->flux, pGrid->npf,
                             MP_RL, this_ray->proc_list[p], flux_tag,
                             MPI_COMM_WORLD);
              if (err) ath_error("[update_tree_radiation]: MPI_Send error = %d\n",err);
              }
          }
        }
/* ------------------------- End initialize flux on ray ------------------------------ *


/* Walk along ray */
        for (l=0; l<this_ray->ncell[p]; l++) {
          /* Step to next */

          last_length = cum_length;
          cum_length += this_ray->cell_length_list[p][l];

/* Store moments */
          i = this_ray->i_list[p][l];
          j = this_ray->j_list[p][l];
          k = this_ray->k_list[p][l];
	  cc_pos(pGrid,i,j,k,&x1c,&x2c,&x3c);

/* Used distance from source to cell center to evaluate solid angle fraction */
          rcell2 = (SQR(x1c-trees[n].x1)+SQR(x2c-trees[n].x2)+SQR(x3c-trees[n].x3)) *
	           4 * PI * pGrid->ray_weight[n][k][j][i];

/* For multiple rays in cell, weight by ray length and solid angle.  Must
   normailze with ray_weight */
	  for(ipf=0; ipf < pGrid->npf; ipf++) {	    
	    Jray = pGrid->sourcelist[n].s * this_ray->flux[ipf] * this_ray->cell_length_list[p][l] * ray_sa / rcell2;
	    pGrid->Jps[ipf][k][j][i] += Jray;
	    pGrid->Hps[ipf][k][j][i][0] += Jray * this_ray->n1;
	    pGrid->Hps[ipf][k][j][i][1] += Jray * this_ray->n2;
	    pGrid->Hps[ipf][k][j][i][2] += Jray * this_ray->n3;
	  /*          pGrid->Kps[k][j][i][0] += Jray * SQR(this_ray->n1);
		      pGrid->Kps[k][j][i][1] += Jray * this_ray->n1 * this_ray->n2;
		      pGrid->Kps[k][j][i][2] += Jray * SQR(this_ray->n2);
		      pGrid->Kps[k][j][i][3] += Jray * this_ray->n1 * this_ray->n3;
		      pGrid->Kps[k][j][i][4] += Jray * this_ray->n2 * this_ray->n3;
		      pGrid->Kps[k][j][i][5] += Jray * SQR(this_ray->n3);*/
/* Reduce the flux */
	    if (last_length > dummy_r) {	   
	      this_ray->flux[ipf] *= this_ray->etau_list[p][ipf][l];	
	    }
	  } /*end loop over frequencies */
        } /* End walk along ray */
      } /* End loop over processors /*

/* Record fraction of flux left */
      for(ipf=0; ipf < pGrid->npf; ipf++) {
	if(this_ray->flux !=NULL)
	  flux_frac = this_ray->flux[ipf];
	//	flux_frac = this_ray->flux[ipf] / pGrid->sourcelist[n].s;
	if (flux_frac > max_frac) max_frac = flux_frac;
      }
  } /* End second loop over rays */

/* See if we need to proceed to the next level. We only bother
 * continuing if at least FRAC_FLUX of the original photons are
 * left.
 */
    err = MPI_Allreduce(&max_frac, &max_frac_glob, 1, MP_RL, MPI_MAX,
                        MPI_COMM_WORLD);
    if (err) ath_error("[update_tree_radiation]: MPI_Allreduce error = %d\n",err);    
    //if (myID_Comm_world == 0) printf("%g %g\n",max_frac_glob,FRAC_FLUX);
    if (max_frac_glob < FRAC_FLUX) break;
  } /* End loop over levels */
}
#else /* Not MPI_PARALLEL */
void update_tree_radiation(int n, GridS *pGrid)
{
  int rayptr, lev, l;
  Ray *this_ray, *parent;
  Real tau, cum_length, last_length;
  Real dummy_r, max_flux, flux_frac;
  Real density;
  int nray;
  int i,j,k,m,ipf;
  Real ray_area, ray_sa, proj;
  Real x1c,x2c,x3c, rcell2;
  Real Jray;

/* Dummy radius, used in initializing fluxes. Its value doesn't
 * matter, it is just used for coding convenience.
 */
  dummy_r = 1.0e-6 * pGrid->dx1;
  
/* Traverse tree */
  rayptr=0;
  while (rayptr >= 0) {

/* Level, number of rays, and ray solid angle */
    lev = ray_level(rayptr);
    nray = 12*(1<<(2*lev));
    ray_sa = 4 * PI / nray;
    //ray_sa = 1.0 / nray;
    //theta0 = sqrt(ray_sa);

/* Pointer to this ray and its parent*/
    this_ray = &(trees[n].rays[rayptr]);
    if (lev != 0) parent = &(trees[n].rays[PARENT(rayptr)]);

/* Initialize ray length and flux */
    if (lev == 0) {
      cum_length = dummy_r;
      //      this_ray->flux = pGrid->sourcelist[n].s / (4.0*PI*SQR(dummy_r));
      for(ipf=0; ipf < pGrid->npf; ipf++) {
	this_ray->flux[ipf] = 1.0;
      }
    } else {
      cum_length = parent->cum_length;
      for(ipf=0; ipf < pGrid->npf; ipf++) {
	this_ray->flux[ipf] = parent->flux[ipf];
      }
    }

    flux_frac = 1.0;
/* Walk ray */
    for (l=0; l<this_ray->ncell; l++) {

      i = this_ray->i_list[l];
      j = this_ray->j_list[l];
      k = this_ray->k_list[l];
      cc_pos(pGrid,i,j,k,&x1c,&x2c,&x3c);

/* Compute exp of optical depth in cell */
      density = pGrid->U[k][j][i].d;
      for(ipf=0; ipf < pGrid->npf; ipf++) {
	this_ray->etau_list[ipf][l] =
	  exp(-PSOpacity(pGrid,ipf,i,j,k) * density * this_ray->cell_length_list[l]);
      }
/* Wise et al. covering factor. Maybe useful for photoionization but
   not currently used. */
/*      last_length = cum_length;
      cum_length += this_ray->cell_length_list[l];
      lmax = this_ray->n1;
      midp = last_length + 0.5 * this_ray->cell_length_list[l];
      dxedge[0] = 0.5 * pGrid->dx1 - fabs(trees[n].x1 + this_ray->n1 * midp - x1c);
      dxedge[1] = 0.5 * pGrid->dx2 - fabs(trees[n].x2 + this_ray->n2 * midp - x2c);
      dxedge[2] = 0.5 * pGrid->dx3 - fabs(trees[n].x3 + this_ray->n3 * midp - x3c);
      dedge = dxedge[0];
      for (m=1;m<3;m++){
        if (dxedge[m] < dedge) dedge = dxedge[m];
      }
      
      cov_fac = MIN(0.5+dedge/(theta0 * last_length),1.0);*/


      last_length = cum_length;
      cum_length += this_ray->cell_length_list[l];

/* compute moments */  
      cc_pos(pGrid,i,j,k,&x1c,&x2c,&x3c);
      rcell2 = (SQR(x1c-trees[n].x1)+SQR(x2c-trees[n].x2)+SQR(x3c-trees[n].x3)) *
	       4 * PI * pGrid->ray_weight[n][k][j][i];
      //printf("%d %g\n",n,pGrid->sourcelist[n].s);
      for(ipf=0; ipf < pGrid->npf; ipf++) {
	Jray = pGrid->sourcelist[n].s * this_ray->flux[ipf] * this_ray->cell_length_list[l] * ray_sa / rcell2;
	pGrid->Jps[ipf][k][j][i] += Jray;
	pGrid->Hps[ipf][k][j][i][0] += Jray * this_ray->n1;
	pGrid->Hps[ipf][k][j][i][1] += Jray * this_ray->n2;
	pGrid->Hps[ipf][k][j][i][2] += Jray * this_ray->n3;
	/*pGrid->Kps[k][j][i][0] += Jray * SQR(this_ray->n1);
	  pGrid->Kps[k][j][i][1] += Jray * this_ray->n1 * this_ray->n2;
	  pGrid->Kps[k][j][i][2] += Jray * SQR(this_ray->n2);
	  pGrid->Kps[k][j][i][3] += Jray * this_ray->n1 * this_ray->n3;
	  pGrid->Kps[k][j][i][4] += Jray * this_ray->n2 * this_ray->n3;
	  pGrid->Kps[k][j][i][5] += Jray * SQR(this_ray->n3);*/

/* Reduce the flux */
	if (last_length > dummy_r)
	  this_ray->flux[ipf] *=  this_ray->etau_list[ipf][l];
      }
/* Find flux left on this ray. If it's less than FRAC_FLUX of the
 * original, stop following this ray. */

      flux_frac = this_ray->flux[0] / pGrid->sourcelist[n].s;
      if (flux_frac < FRAC_FLUX) break;
    }

/* Go to next ray */
    if (flux_frac >= FRAC_FLUX) rayptr = next_ray(n, rayptr);
    else rayptr = next_non_child_ray(rayptr);

  } /* End loop over rays */
}
#endif /* End not MPI_PARALLEL */




#endif /* POINT_SOURCE */
