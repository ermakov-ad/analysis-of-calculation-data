//
// Created by petchu on 13.05.2022.
//
#include "ex5.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "ex.h"
#include "struct.h"
#ifndef CODE_SECOND_STEP_H
#define CODE_SECOND_STEP_H
void calc_second_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_first_step_3D, double ***v_first_step_3D,
                      double ***w_first_step_3D, double *p_prev_3D, double *x_grid, double *y_grid, double *z_grid, double dt, int myid, int num_procs,
					  HYPRE_IJMatrix *AA, HYPRE_ParCSRMatrix *parcsr_A);
HYPRE_ParCSRMatrix calc_second_step_matrix(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_first_step_3D, double ***v_first_step_3D,
                      double ***w_first_step_3D, double *p_prev_3D, double *x_grid, double *y_grid, double *z_grid, double dt, int myid, int num_procs);
void fin_second_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double *p_prev_3D, int myid, int num_procs, HYPRE_IJVector *x);
#endif //CODE_SECOND_STEP_H
