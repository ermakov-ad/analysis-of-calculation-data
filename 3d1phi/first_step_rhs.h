#include "ex5.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "struct.h"
#include "first_step.h"
void calc_first_step_rhs(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***p_next_3D, double ***u_inter_3D, double ***v_inter_3D, double ***w_inter_3D , double ***u_prev_3D, double ***v_prev_3D, double ***w_prev_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, double *x_grid, double *y, double *z, double dt, double time, int myid, int num_procs, HYPRE_IJVector *x,
					 HYPRE_ParVector *par_x, HYPRE_IJVector *b, HYPRE_ParVector *par_b, double ***x_nonlin, double ***y_nonlin, double ***z_nonlin, double ***x_force,
					 double ***y_force, double ***z_force);
void init_nonlin( double ***x_nonlin, double ***y_nonlin, double ***z_nonlin, int nx, int ny, int nz);
void init_force(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***x_force, double ***y_force, double ***z_force, double N, double *x_grid, double *y, double *z, double dt);
