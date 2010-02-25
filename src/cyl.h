/* NOTE:  LATER, WE CAN ADD THESE PROTOTYPES TO THE FILE prototypes.h AND DELETE THIS FILE... */

#ifndef CYL_H
#define CYL_H 

#ifdef CYLINDRICAL


/*============================================================================
 * CONVERSION FUNCTIONS
 *============================================================================*/
// void Cons_to_Prim(const ConsS *pU, PrimS *pW);
void Prim_to_Cons(ConsS *pU, const PrimS *pW);

/*============================================================================
 * BOUNDARY CONDITION FUNCTIONS
 *============================================================================*/
void do_nothing_bc(GridS *pG);
void diode_outflow_ix1(GridS *pGrid);
void diode_outflow_ox1(GridS *pGrid);
void diode_outflow_ix3(GridS *pGrid);
void diode_outflow_ox3(GridS *pGrid);

/*============================================================================
 * ERROR-ANALYSIS FUNCTIONS
 *============================================================================*/
Real compute_div_b(GridS *pG);
void compute_l1_error(char *problem, MeshS *pM, ConsS ***RootSoln);

/*============================================================================
 * ROOT-FINDING FUNCTIONS
 *============================================================================*/
int sign_change(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *a, Real *b);
int bisection(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *root);

#endif /* CYLINDRICAL */
#endif /* CYL_H */
