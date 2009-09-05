#ifndef PROTOTYPES_H
#define PROTOTYPES_H 
#include "copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions from the /src directory,
 *   and all subdirectories in /src
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "athena.h"
#include "defs.h"
#include "config.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

/* Include prototypes in /src sub-directories */

#ifdef ION_RADIATION
#include "ionradiation/prototypes.h"
#endif

#ifdef FFT_ENABLED
#include "fftsrc/prototypes.h"
#endif

#include "microphysics/prototypes.h"
#include "integrators/prototypes.h"
#include "reconstruction/prototypes.h"
#include "rsolvers/prototypes.h"
#include "particles/prototypes.h"

/*----------------------------------------------------------------------------*/
/* main.c */
int athena_main(int argc, char *argv[]);

/*----------------------------------------------------------------------------*/
/* ath_array.c */
void*   calloc_1d_array(                      size_t nc, size_t size);
void**  calloc_2d_array(           size_t nr, size_t nc, size_t size);
void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);
void free_1d_array(void *array);
void free_2d_array(void *array);
void free_3d_array(void *array);

/*----------------------------------------------------------------------------*/
/* ath_log.c */
void ath_log_set_level(const int out, const int err);
void ath_log_open(const char *basename, const int lazy, const char *mode);
void ath_log_close(void);
FILE *athout_fp(void);
FILE *atherr_fp(void);
void ath_flush_out(void);
void ath_flush_err(void);
int ath_perr(const int level, const char *fmt, ...);
int ath_pout(const int level, const char *fmt, ...);

/*----------------------------------------------------------------------------*/
/* ath_files.c */
char *ath_fname(const char *path, const char *basename,
                const int dlen, const int idump, 
                const char *id, const char *ext);
FILE *ath_fopen(const char *path, const char *basename,
                const int dlen, const int idump, 
                const char *id, const char *ext, const char *mode);
size_t ath_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);

/*----------------------------------------------------------------------------*/
/* ath_signal.c */
void ath_sig_init(void);
int  ath_sig_act(int *piquit);

/*----------------------------------------------------------------------------*/
/* baton.c */
void baton_start(const int Nb, const int tag);
void baton_stop(const int Nb, const int tag);

/*----------------------------------------------------------------------------*/
/* cc_pos.c */
void cc_pos(const Grid *pG, const int i, const int j,const int k,
            Real *px1, Real *px2, Real *px3);
void vc_pos(const Grid *pG, const int i, const int j,const int k,
            Real *px1, Real *px2, Real *px3);
#ifdef PARTICLES
int celli(const Grid* pGrid, const Real x, const Real dx1_1, int *i, Real *a);
Real x1cc(const Grid* pGrid, const int i);
int cellj(const Grid* pGrid, const Real y, const Real dx2_1, int *j, Real *b);
Real x2cc(const Grid* pGrid, const int j);
int cellk(const Grid* pGrid, const Real z, const Real dx3_1, int *k, Real *c);
Real x3cc(const Grid* pGrid, const int k);
#endif

/*----------------------------------------------------------------------------*/
/* convert_var.c */
void Cons1D_to_Prim1D(const Cons1D *pU, Prim1D *pW MHDARG( , const Real *pBx));
void Prim1D_to_Cons1D(Cons1D *pU, const Prim1D *pW MHDARG( , const Real *pBx));
Real cfast(const Cons1D *U MHDARG( , const Real *Bx));

/*----------------------------------------------------------------------------*/
/* init_domain.c */
void init_domain(Grid *pG, Domain *pD);
void get_myGridIndex(Domain *pD, const int my_id, int *pi, int *pj, int *pk);

/*----------------------------------------------------------------------------*/
/* init_grid.c */
void init_grid(Grid *pGrid, Domain *pD);

/*----------------------------------------------------------------------------*/
/* new_dt.c */
void new_dt(Grid *pGrid);
#ifdef MPI_PARALLEL
void sync_dt(Grid *pG);
#endif

/*----------------------------------------------------------------------------*/
/* output.c - and related files */
void init_output(Grid *pGrid);
void data_output(Grid *pGrid, Domain *pD, const int flag);
int  add_output(Output *new_out);
void add_rst_out(Output *new_out);
void data_output_destruct(void);
void dump_history_enroll(const Gasfun_t pfun, const char *label);
void data_output_enroll(Real time, Real dt, int num, const VGFunout_t fun,
			const char *fmt, const Gasfun_t expr, int n,
			const Real dmin, const Real dmax, int sdmin, int sdmax
#ifdef PARTICLES
			, const int out_pargrid, PropFun_t par_prop
#endif
);
float ***subset3(Grid *pGrid, Output *pout);
float  **subset2(Grid *pGrid, Output *pout);
float   *subset1(Grid *pGrid, Output *pout);

void output_fits (Grid *pGrid, Domain *pD, Output *pOut);
void output_pdf  (Grid *pGrid, Domain *pD, Output *pOut);
void output_pgm  (Grid *pGrid, Domain *pD, Output *pOut);
void output_ppm  (Grid *pGrid, Domain *pD, Output *pOut);
void output_vtk  (Grid *pGrid, Domain *pD, Output *pOut);
void output_tab  (Grid *pGrid, Domain *pD, Output *pOut);

void dump_binary (Grid *pGrid, Domain *pD, Output *pOut);
void dump_dx     (Grid *pGrid, Domain *pD, Output *pOut);
void dump_history(Grid *pGrid, Domain *pD, Output *pOut);
void dump_tab_cons(Grid *pGrid, Domain *pD, Output *pOut);
void dump_tab_prim(Grid *pGrid, Domain *pD, Output *pOut);
void dump_vtk    (Grid *pGrid, Domain *pD, Output *pOut);

/*----------------------------------------------------------------------------*/
/* par.c */
void   par_open(char *filename);
void   par_cmdline(int argc, char *argv[]);
int    par_exist(char *block, char *name);

char  *par_gets(char *block, char *name);
int    par_geti(char *block, char *name);
double par_getd(char *block, char *name);

char  *par_gets_def(char *block, char *name, char   *def);
int    par_geti_def(char *block, char *name, int    def);
double par_getd_def(char *block, char *name, double def);

void   par_sets(char *block, char *name, char *sval, char *comment);
void   par_seti(char *block, char *name, char *fmt, int ival, char *comment);
void   par_setd(char *block, char *name, char *fmt, double dval, char *comment);

void   par_dump(int mode, FILE *fp);
void   par_close(void);

#ifdef MPI_PARALLEL
void par_dist_mpi(const int mytid, MPI_Comm comm);
#endif

/*----------------------------------------------------------------------------*/
/* prob/PROBLEM.c ; linked to problem.c */
void problem(Grid *pgrid, Domain *pDomain);
void Userwork_in_loop(Grid *pgrid, Domain *pDomain);
void Userwork_after_loop(Grid *pgrid, Domain *pDomain);
void problem_read_restart(Grid *pG, Domain *pD, FILE *fp);
void problem_write_restart(Grid *pG, Domain *pD, FILE *fp);
Gasfun_t get_usr_expr(const char *expr);
VGFunout_t get_usr_out_fun(const char *name);
#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name);
void gasvshift(const Real x1, const Real x2, const Real x3, Real *u1, Real *u2, Real *u3);
void Userforce_particle(Vector *ft, const Real x1, const Real x2, const Real x3, Real *w1, Real *w2, Real *w3);
#endif

/*----------------------------------------------------------------------------*/
/* restart.c  */
void dump_restart(Grid *pG, Domain *pD, Output *pout);
void restart_grid_block(char *res_file, Grid *pGrid, Domain *pDomain);

/*----------------------------------------------------------------------------*/
/* self_gravity.c  */
#ifdef SELF_GRAVITY
VGDFun_t selfg_init(Grid *pG, Domain *pD);
void selfg_flux_correction(Grid *pG);
#endif
void selfg_by_multig_1d(Grid *pG, Domain *pD);
void selfg_by_multig_2d(Grid *pG, Domain *pD);
void selfg_by_multig_3d(Grid *pG, Domain *pD);
void selfg_by_multig_3d_init(Grid *pG, Domain *pD);
#if defined(FFT_ENABLED) && defined(SELF_GRAVITY_USING_FFT)
void selfg_by_fft_1d(Grid *pG, Domain *pD);
void selfg_by_fft_2d(Grid *pG, Domain *pD);
void selfg_by_fft_3d(Grid *pG, Domain *pD);
void selfg_by_fft_2d_init(Grid *pG, Domain *pD);
void selfg_by_fft_3d_init(Grid *pG, Domain *pD);
#endif /* FFT_ENABLED */

/*----------------------------------------------------------------------------*/
/* set_bvals_mhd.c  */
void set_bvals_mhd_init(Grid *pG, Domain *pD);
void set_bvals_mhd_fun(enum Direction dir, VBCFun_t prob_bc);
void set_bvals_mhd(Grid *pGrid, Domain *pDomain);

/*----------------------------------------------------------------------------*/
/* set_bvals_grav.c  */
#ifdef SELF_GRAVITY
void set_bvals_grav_init(Grid *pG, Domain *pD);
void set_bvals_grav_fun(enum Direction dir, VBCFun_t prob_bc);
void set_bvals_grav(Grid *pGrid, Domain *pDomain);
#endif

/*----------------------------------------------------------------------------*/
/* set_bvals_shear.c  */
#ifdef SHEARING_BOX
void ShearingSheet_ix1(Grid *pG, Domain *pD);
void ShearingSheet_ox1(Grid *pG, Domain *pD);
void RemapEy_ix1(Grid *pG, Domain *pD, Real ***emfy, Real **remapEyiib);
void RemapEy_ox1(Grid *pG, Domain *pD, Real ***emfy, Real **remapEyoib);
void set_bvals_shear_init(Grid *pG, Domain *pD);
void set_bvals_shear_destruct(void);
#ifdef FARGO
void Fargo(Grid *pG, Domain *pD);
#endif
#endif /* SHEARING_BOX */

/*----------------------------------------------------------------------------*/
/* show_config.c */
void show_config(void);
void show_config_par(void);

/*----------------------------------------------------------------------------*/
/* utils.c */
char *ath_strdup(const char *in);
int ath_gcd(int a, int b);
int ath_big_endian(void);
void ath_bswap(void *vdat, int sizeof_len, int cnt);
void ath_error(char *fmt, ...);
void minmax1(float *data, int nx1, float *dmin, float *dmax);
void minmax2(float **data, int nx2, int nx1, float *dmin, float *dmax);
void minmax3(float ***data, int nx3, int nx2, int nx1,float *dmin,float *dmax);

#endif /* PROTOTYPES_H */
