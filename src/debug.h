// NOTE:  THIS FILE DOES NOT NEED TO BE ADDED TO THE MAIN DISTRIBUTION, SO 
//        WE WILL NOT PLACE THESE PROTOTYPES IN THE FILE prototypes.h

#ifndef DEBUG_H
#define DEBUG_H 

/*============================================================================
 * DEFINITIONS
 *============================================================================*/
#define IVIEW 68
#define JVIEW 4
#define KVIEW 4
#define VIEW1D (i==IVIEW) && (pG->my_id==0)
#define VIEW2D (i==IVIEW) && (j==JVIEW) && (pG->my_id==0)
#define VIEW3D (i==IVIEW) && (j==JVIEW) && (k==KVIEW) && (pG->my_id==0)

#define BAR_LENGTH 70
#define BAR_SYMBOL '='

#define FMT "%18.15f"

// #ifdef CYLINDRICAL

/*============================================================================
 * DISPLAY FUNCTIONS
 *============================================================================*/
void debug_header(const int level, char *section_name);
int print_cons1d(const int level, char *name, const Cons1D *U,
                 const int k, const int j, const int i, const int dir);
int print_prim1d(const int level, char *name, const Prim1D *W,
                 const int k, const int j, const int i, const int dir);
int print_gas(const int level, char *name, const Gas *U,
              const int k, const int j, const int i);

// #endif /* CYLINDRICAL */

#endif /* DEBUG_H */
