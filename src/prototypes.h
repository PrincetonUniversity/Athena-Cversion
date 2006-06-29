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
 *   ath_parallel.c
 *   ath_signal.c
 *   dump_binary.c, dump_dx.c, dump_history.c, dump_table.c, dump_vtk.c
 *   init_grid_block.c
 *   output.c
 *   output_fits.c, output_pdf.c output_pgm.c, output_ppm.c, output_tab.c
 *   par.c
 *   restart.c
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

/*----------------------------------------------------------------------------*/
/* main.c */

int athena_main(int argc, char *argv[]);

/*----------------------------------------------------------------------------*/
/* ath_array.c */

void*   calloc_1d_array(                      size_t nc, size_t size);
void**  calloc_2d_array(           size_t nr, size_t nc, size_t size);
void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);

void free_1d_array(void   *array);
void free_2d_array(void  **array);
void free_3d_array(void ***array);

/*----------------------------------------------------------------------------*/
/* ath_files.c */

char *fname_construct(const char *basename, const int dlen, const int idump, 
		      const char *type, const char *ext);

FILE *ath_fopen(const char *basename, const int dlen, const int idump, 
		const char *type, const char *ext, const char *mode);

/*----------------------------------------------------------------------------*/
/* ath_parallel.c */

#ifdef MPI_PARALLEL
void domain_partition(Grid *pG);
void domain_destruct(void);
void domain_ijk(const int my_id, int *pi, int *pj, int *pk);
#endif

/*----------------------------------------------------------------------------*/
/* ath_signal.c */

void ath_sig_init(void);
int  ath_sig_act(Grid *pG);

/*----------------------------------------------------------------------------*/
/* dump_*.c (multiple files) */

void dump_binary(Grid *pGrid, Output *pOut);
void dump_dx(Grid *pGrid, Output *pOut);
void dump_history(Grid *pGrid, Output *pOut);
void dump_history_enroll(const Gasfun_t pfun, const char *label);
void dump_table(Grid *pGrid, Output *pOut);
void dump_vtk(Grid *pGrid, Output *pOut);

/*----------------------------------------------------------------------------*/
/* init_grid_block.c */

void init_grid_block(Grid *pGrid);

/*----------------------------------------------------------------------------*/
/* output.c - and related files */

void init_output(Grid *pGrid);
void data_output(Grid *pgrid, const int flag);
int  add_output(Output *new_out);
void add_rst_out(Output *new_out);
void data_output_destruct(void);
void data_output_enroll(Real time, Real dt, int num, const VGFunout_t fun,
	                const char *fmt,  const Gasfun_t expr, int n,
                        const Real dmin, const Real dmax, int sdmin, int sdmax);
float **subset2(Grid *pgrid, Output *pout);
float  *subset1(Grid *pgrid, Output *pout);

void output_fits(Grid *pGrid, Output *pOut);
void output_pdf(Grid *pGrid, Output *pOut);
void output_pgm(Grid *pGrid, Output *pOut);
void output_ppm(Grid *pGrid, Output *pOut);
void output_tab(Grid *pGrid, Output *pOut);

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

void problem(Grid *pgrid);
void Userwork_before_loop(Grid *pgrid);
void Userwork_in_loop(Grid *pgrid);
void Userwork_after_loop(Grid *pgrid);
void problem_read_restart(Grid *pG, FILE *fp);
void problem_write_restart(Grid *pG, FILE *fp);
Gasfun_t get_usr_expr(const char *expr);

/*----------------------------------------------------------------------------*/
/* restart.c  */

void dump_restart(Grid *pG, Output *pout);
void restart_grid_block(char *res_file, Grid *pGrid);

/*----------------------------------------------------------------------------*/
/* show_config.c */

void show_config(void);

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
