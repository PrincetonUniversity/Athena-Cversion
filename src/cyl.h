/* NOTE:  LATER, WE CAN ADD THESE PROTOTYPES TO THE FILE prototypes.h AND DELETE THIS FILE... */

#ifndef CYL_H
#define CYL_H 

#ifdef CYLINDRICAL


/*============================================================================
 * CONVERSION FUNCTIONS
 *============================================================================*/
void Gas_to_Prim(const Gas *pU, Prim *pW);
void Prim_to_Gas(Gas *pU, const Prim *pW);

/*============================================================================
 * BOUNDARY CONDITION FUNCTIONS
 *============================================================================*/
void do_nothing_bc(Grid *pG);
void diode_outflow_ix1(Grid *pGrid);
void diode_outflow_ox1(Grid *pGrid);
void diode_outflow_ix3(Grid *pGrid);
void diode_outflow_ox3(Grid *pGrid);

/*============================================================================
 * ERROR-ANALYSIS FUNCTIONS
 *============================================================================*/
Real compute_div_b(Grid *pG);
void compute_l1_error(char *problem, Grid *pG, Domain *pDomain, Gas ***Soln, const int errortype);

/*============================================================================
 * ROOT-FINDING FUNCTIONS
 *============================================================================*/
int sign_change(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *a, Real *b);
int bisection(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *root);

#endif /* CYLINDRICAL */
#endif /* CYL_H */
