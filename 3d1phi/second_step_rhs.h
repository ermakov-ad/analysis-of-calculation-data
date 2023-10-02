#include "ex5.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "struct.h"
#include "second_step.h"
void calc_second_step_rhs(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_first_step_3D, double ***v_first_step_3D,
                      double ***w_first_step_3D, double *p_prev_3D, double *x_grid, double *y_grid, double *z_grid, double dt, int myid, int num_procs, HYPRE_IJVector *x,
						 HYPRE_ParVector *par_x, HYPRE_IJVector *b, HYPRE_ParVector *par_b);
