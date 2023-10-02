#include "ex5.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "struct.h"

void calc_first_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***p_next_3D, double ***u_inter_3D, double ***v_inter_3D, double ***w_inter_3D, double ***u_prev_3D, double ***v_prev_3D, double ***w_prev_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, double *x_grid, double *y, double *z, double dt, double time, int myid, int num_procs, HYPRE_IJVector *x,
					 HYPRE_ParVector *par_x, HYPRE_IJVector *b, HYPRE_ParVector *par_b, HYPRE_IJMatrix *AA, HYPRE_ParCSRMatrix *parcsr_A);
HYPRE_ParCSRMatrix calc_first_step_matrix(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_inter_3D, double ***v_inter_3D, double ***w_inter_3D , double ***u_prev_3D, double ***v_prev_3D, double ***w_prev_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, double *x_grid, double *y, double *z, double dt, double time, int myid, int num_procs);

double calc_x_nonlinear_term(double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz );

double calc_y_nonlinear_term(double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz );

double calc_z_nonlinear_term(double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz );

double calc_x_visc_term(struct ParametersCommon paramsc, double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz );

double calc_y_visc_term(struct ParametersCommon paramsc, double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz );

double calc_z_visc_term(struct ParametersCommon paramsc, double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz );

void fin_first_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***p_next_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, int myid, int num_procs, HYPRE_IJVector *x);
