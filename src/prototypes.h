#ifndef PROTOTYPES_H
#define PROTOTYPES_H 
#include "copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions from the following files:
 *   main.c
 *   ath_array.c
 *   ath_files.c
 *   ath_log.c
 *   ath_signal.c
 *   baton.c
 *   cc_pos.c
 *   convert_var.c
 *   esystem_prim.c, esystem_roe.c
 *   flux_force.c, flux_hllc.c, flux_hlld.c, flux_hlle.c, flux_roe.c, 
 *     flux_2shock.c, flux_exact
 *   init_domain.c
 *   init_grid.c
 *   integrate.c
 *   integrate_1d.c, integrate_2d.c, integrate_3d-vl., integrate_3d-ctu.c
 *   lr_states_prim1.c, lr_states_prim2.c, lr_states_prim3.c
 *   new_dt.c
 *   output.c
 *   output_fits.c, output_pdf.c output_pgm.c, output_ppm.c, output_tab.c
 *   dump_binary.c, dump_dx.c, dump_history.c, dump_table.c, dump_vtk.c
 *   par.c
 *   restart.c
 *   set_bvals_mhd.c
 *   set_bvals_grav.c
 *   set_bvals_shear.c
 *   show_config.c
 *   utils.c
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "athena.h"
#include "defs.h"

#include "config.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

/* Include prototypes for ionization files if ionization is turned on */
#ifdef ION_RADIATION
#include "ionrad_prototypes.h"
#endif

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

/*----------------------------------------------------------------------------*/
/* convert_var.c */
void Cons1D_to_Prim1D(const Cons1D *pU, Prim1D *pW MHDARG( , const Real *pBx));
void Prim1D_to_Cons1D(Cons1D *pU, const Prim1D *pW MHDARG( , const Real *pBx));
Real cfast(const Cons1D *U MHDARG( , const Real *Bx));

/*----------------------------------------------------------------------------*/
/* esystem_*.c */
void esys_prim_iso_hyd(const Real d, const Real v1,
  Real eigenvalues[],
  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4]);

void esys_prim_adb_hyd(const Real d, const Real v1, const Real p,
  Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5]);

void esys_prim_iso_mhd(const Real d, const Real v1, const Real b1,
  const Real b2, const Real b3, Real eigenvalues[],
  Real right_eigenmatrix[][6], Real left_eigenmatrix[][6]);

void esys_prim_adb_mhd(const Real d, const Real v1, const Real p,
  const Real b1, const Real b2, const Real b3, Real eigenvalues[],
  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7]);

void esys_roe_iso_hyd(const Real v1, const Real v2, const Real v3,
  Real eigenvalues[],
  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4]);

void esys_roe_adb_hyd(const Real v1, const Real v2, const Real v3,
  const Real h, Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5]);

void esys_roe_iso_mhd(const Real d, const Real v1, const Real v2,
  const Real v3, const Real b1, const Real b2, const Real b3,
  const Real x, const Real y, Real eigenvalues[],
  Real right_eigenmatrix[][6], Real left_eigenmatrix[][6]);

void esys_roe_adb_mhd(const Real d, const Real v1, const Real v2,
  const Real v3, const Real h, const Real b1, const Real b2, const Real b3,
  const Real x, const Real y, Real eigenvalues[],
  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7]);

/*----------------------------------------------------------------------------*/
/* flux_*.c */
void flux_force (const Cons1D Ul, const Cons1D Ur,
                 const Prim1D Wl, const Prim1D Wr,
                 MHDARG(const Real Bxi,) Cons1D *pF);
void flux_hllc  (const Cons1D Ul, const Cons1D Ur,
                 const Prim1D Wl, const Prim1D Wr, Cons1D *pF);
void flux_hlld  (const Cons1D Ul, const Cons1D Ur,
                 const Prim1D Wl, const Prim1D Wr, const Real Bxi, Cons1D *pF);
void flux_hlle  (const Cons1D Ul, const Cons1D Ur,
                 const Prim1D Wl, const Prim1D Wr,
                 MHDARG(const Real Bxi,) Cons1D *pF);
void flux_roe   (const Cons1D Ul, const Cons1D Ur,
                 const Prim1D Wl, const Prim1D Wr,
                 MHDARG(const Real Bxi,) Cons1D *pF);
void flux_2shock(const Cons1D Ul, const Cons1D Ur,
                 const Prim1D Wl, const Prim1D Wr,
                 MHDARG(const Real Bxi,) Cons1D *pF);
void flux_exact (const Cons1D Ul, const Cons1D Ur,
                 const Prim1D Wl, const Prim1D Wr,
                 MHDARG(const Real Bxi,) Cons1D *pF);

/*----------------------------------------------------------------------------*/
/* init_domain.c */
void init_domain(Grid *pG, Domain *pD);
void get_myGridIndex(Domain *pD, const int my_id, int *pi, int *pj, int *pk);

/*----------------------------------------------------------------------------*/
/* init_grid.c */
void init_grid(Grid *pGrid, Domain *pD);

/*----------------------------------------------------------------------------*/
/* integrate.c */
VGDFun_t integrate_init(int Nx1, int Nx2, int Nx3);
void integrate_destruct(void);

/*----------------------------------------------------------------------------*/
/* integrate_1d.c */
void integrate_destruct_1d(void);
void integrate_init_1d(int Nx1);
void integrate_1d(Grid *pG, Domain *pD);

/*----------------------------------------------------------------------------*/
/* integrate_2d.c */
void integrate_destruct_2d(void);
void integrate_init_2d(int Nx1, int Nx2);
void integrate_2d(Grid *pG, Domain *pD);

/*----------------------------------------------------------------------------*/
/* integrate_3d.c */
void integrate_destruct_3d(void);
void integrate_init_3d(int Nx1, int Nx2, int Nx3);
void integrate_3d_vl(Grid *pG, Domain *pD);
void integrate_3d_ctu(Grid *pG, Domain *pD);

/*----------------------------------------------------------------------------*/
/* lr_states.c */
void lr_states_destruct(void);
void lr_states_init(int nx1, int nx2, int nx3);
void lr_states(const Prim1D W[], MHDARG( const Real Bxc[] , )
               const Real dt, const Real dtodx, const int is, const int ie,
               Prim1D Wl[], Prim1D Wr[]);

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
			const Real dmin, const Real dmax, int sdmin, int sdmax);
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
void dump_tab    (Grid *pGrid, Domain *pD, Output *pOut);
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

/*----------------------------------------------------------------------------*/
/* restart.c  */
void dump_restart(Grid *pG, Domain *pD, Output *pout);
void restart_grid_block(char *res_file, Grid *pGrid, Domain *pDomain);

/*----------------------------------------------------------------------------*/
/* self_gravity.c  */
VGDFun_t selfg_init(Grid *pG, Domain *pD);
void selfg_flux_correction(Grid *pG);
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
