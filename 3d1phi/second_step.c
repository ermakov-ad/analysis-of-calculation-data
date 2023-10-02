//
// Created by petchu on 13.05.2022.
//

#include "second_step.h"
#define my_min(a,b)  (((a)<(b)) ? (a) : (b))
void calc_second_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_first_step_3D, double ***v_first_step_3D,
                      double ***w_first_step_3D, double *p_prev_3D, double *x_grid, double *y_grid, double *z_grid, double dt, int myid, int num_procs,
					  HYPRE_IJMatrix *AA, HYPRE_ParCSRMatrix *parcsr_A){
    double hx = x_grid[1] - x_grid[0];
    double hy = y_grid[1] - y_grid[0];
    double hz = z_grid[1] - z_grid[0];
    int i;
    int N, n;

    int ilower, iupper;
    int local_size, extra;

    int solver_id;
    int vis, print_system;

    if (BOUN_TYPE )
        N = ( params2d.cells_number_x - 2 ) * ( params2d.cells_number_y - 2 ) * ( params2d.cells_number_z - 2 );
    else
        N = ( params2d.cells_number_x ) * ( params2d.cells_number_y ) * ( params2d.cells_number_z );
    solver_id = 0;
    vis = 0;
    print_system = 0;

    int arg_index = 0;
    int print_usage = 0;

    local_size = N / num_procs;
    extra = N - local_size * num_procs;

    ilower = local_size * myid;
    ilower += my_min(myid, extra);

    iupper = local_size * (myid + 1);
    iupper += my_min(myid + 1, extra);
    iupper = iupper - 1;

    /* How many rows do I have? */
    local_size = iupper - ilower + 1;
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, AA);

    /* Choose a parallel csr format storage (see the User's Manual) */
    HYPRE_IJMatrixSetObjectType(*AA, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(*AA);

    {
        int nnz;


        double *values = (double *) malloc(7 * sizeof(double));
        int *cols = (int *) malloc(7 * sizeof(int));
        int *tmp = (int *) malloc(2 * sizeof(int));
        //printf("\ncheck\n");
        if (BOUN_TYPE == 1) {
            for (i = ilower; i <= iupper; i++) {

                nnz = 0;
                int k_center = (int) i / ((params2d.cells_number_x - 2) * (params2d.cells_number_y - 2));
                int remainder_k = (int) i % ((params2d.cells_number_x - 2) * (params2d.cells_number_y - 2));
                int j_center = remainder_k / (params2d.cells_number_x - 2);
                int i_center = remainder_k % (params2d.cells_number_x - 2);
                //printf("\n %d %d %d\n", i_center, j_center, k_center);
                if (k_center == ( params2d.cells_number_z - 2 ) / 2 && i_center == ( params2d.cells_number_x - 2 ) / 2 && j_center == ( params2d.cells_number_y - 2 ) / 2) {
                    cols[nnz] = i;
                    values[nnz] = 1;
                    nnz++;
                } else {

                    if (k_center > 0) {
                        cols[nnz] = (k_center - 1) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center) * (params2d.cells_number_x - 2) + i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hz / hz;
                        nnz++;
                    } else {
                        cols[nnz] = (params2d.cells_number_z - 3) * (params2d.cells_number_x - 2) *
                                    (params2d.cells_number_y - 2) + (j_center) * (params2d.cells_number_x - 2) +
                                    i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hz / hz;
                        nnz++;
                    }

                    if (j_center > 0) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center - 1) * (params2d.cells_number_x - 2) + i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hy / hy;
                        nnz++;
                    } else {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (params2d.cells_number_y - 3) * (params2d.cells_number_x - 2) + i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hy / hy;
                        nnz++;
                    }

                    if (i_center > 0) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center) * (params2d.cells_number_x - 2) + i_center - 1;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hx / hx;
                        nnz++;
                    } else {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center) * (params2d.cells_number_x - 2) + (params2d.cells_number_x - 3);
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hx / hx;
                        nnz++;
                    }

                    cols[nnz] = i;
                    values[nnz] = -2 * dt / hx / hx - 2 * dt / hy / hy - 2 * dt / hz / hz;
                    nnz++;


                    if (k_center < params2d.cells_number_z - 3) {
                        cols[nnz] = (k_center + 1) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center) * (params2d.cells_number_x - 2) + i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hz / hz;
                        nnz++;
                    } else {
                        cols[nnz] = (0) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center) * (params2d.cells_number_x - 2) + i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hz / hz;
                        nnz++;
                    }


                    if (j_center < params2d.cells_number_y - 3) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center + 1) * (params2d.cells_number_x - 2) + i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hy / hy;
                        nnz++;
                    } else {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (0) * (params2d.cells_number_x - 2) + i_center;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hy / hy;
                        nnz++;
                    }

                    if (i_center < params2d.cells_number_x - 3) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center) * (params2d.cells_number_x - 2) + i_center + 1;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hx / hx;
                        nnz++;
                    } else {
                        cols[nnz] = (k_center) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j_center) * (params2d.cells_number_x - 2) + 0;
                        //printf("\n%d\n", cols[nnz]);
                        values[nnz] = dt / hx / hx;
                        nnz++;
                    }
                }
                /* Set the values for row i */
                tmp[0] = nnz;
                tmp[1] = i;
                HYPRE_IJMatrixSetValues(*AA, 1, &tmp[0], &tmp[1], cols, values);
            }
        }
        if (BOUN_TYPE == 0){

            for (i = ilower; i <= iupper; i++) {
                //b_values[i] = 0.0;
                nnz = 0;
                int k_center = (int) i / ((params2d.cells_number_x ) * (params2d.cells_number_y ));
                int remainder_k = (int) i % ((params2d.cells_number_x) * (params2d.cells_number_y));
                int j_center = remainder_k / (params2d.cells_number_x);
                int i_center = remainder_k % (params2d.cells_number_x);
                double x_rhs, y_rhs, z_rhs;
                //printf("\n %d %d %d\n", i_center, j_center, k_center);

                if (k_center == 1 && j_center == 1 && i_center == 1) {
                    cols[nnz] = i;
                    values[nnz] = 1.0;
                    nnz++;
                    //b_values[i] += values[nnz];
                } else {

                    if (k_center > 0) {
                        cols[nnz] = (k_center - 1) * (params2d.cells_number_x) * (params2d.cells_number_y) +
                                    (j_center) * (params2d.cells_number_x) + i_center;
                        values[nnz] = dt / hz / hz;
                        //b_values[i] += values[nnz];
                        nnz++;
                    }


                    if (j_center > 0) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x) * (params2d.cells_number_y) +
                                    (j_center - 1) * (params2d.cells_number_x) + i_center;
                        values[nnz] = dt / hy / hy;
                       // b_values[i] += values[nnz];
                        nnz++;
                    }

                    if (i_center > 0) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x) * (params2d.cells_number_y) +
                                    (j_center) * (params2d.cells_number_x) + i_center - 1;
                        values[nnz] = dt / hx / hx;
                        //b_values[i] += values[nnz];
                        nnz++;
                    }

                    cols[nnz] = i;
                    values[nnz] = 0.0;
                    if (i_center == 0)
                        values[nnz] += -dt / hx / hx;

                    else
                    if (i_center == params2d.cells_number_x - 1)
                        values[nnz] += -dt / hx / hx;
                    else
                        values[nnz] += -2 * dt / hx / hx;

                    if (j_center == 0)
                        values[nnz] += -dt / hy / hy;
                    else
                    if (j_center == params2d.cells_number_y - 1)
                        values[nnz] += -dt / hy / hy;
                    else
                        values[nnz] += -2 * dt / hy / hy;

                    if (k_center == 0)
                        values[nnz] += -dt / hz / hz;
                    else
                    if (k_center == params2d.cells_number_z - 1)
                        values[nnz] += -dt / hz / hz;
                    else
                        values[nnz] += -2 * dt / hz / hz;
                    //b_values[i] += values[nnz];
                    nnz++;


                    if (k_center < params2d.cells_number_z - 1) {
                        cols[nnz] = (k_center + 1) * (params2d.cells_number_x) * (params2d.cells_number_y) +
                                    (j_center) * (params2d.cells_number_x) + i_center;
                        values[nnz] = dt / hz / hz;
                        //b_values[i] += values[nnz];
                        nnz++;
                    }


                    if (j_center < params2d.cells_number_y - 1) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x) * (params2d.cells_number_y) +
                                    (j_center + 1) * (params2d.cells_number_x) + i_center;
                        values[nnz] = dt / hy / hy;
                        //b_values[i] += values[nnz];
                        nnz++;
                    }

                    if (i_center < params2d.cells_number_x - 1) {
                        cols[nnz] = (k_center) * (params2d.cells_number_x) * (params2d.cells_number_y) +
                                    (j_center) * (params2d.cells_number_x) + i_center + 1;
                        values[nnz] = dt / hx / hx;
                        //b_values[i] += values[nnz];
                        nnz++;
                    }
                }

                tmp[0] = nnz;
                tmp[1] = i;
                HYPRE_IJMatrixSetValues(*AA, 1, &tmp[0], &tmp[1], cols, values);
            }
        }

        free(values);
        free(cols);
        free(tmp);
    }

    HYPRE_IJMatrixAssemble(*AA);
    HYPRE_IJMatrixGetObject(*AA, (void**) parcsr_A);

}

void fin_second_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double *p_prev_3D, int myid, int num_procs, HYPRE_IJVector *x){

    int i;
    int N, n;

    int ilower, iupper;
    int local_size, extra;

    int solver_id;
    int vis, print_system;

    if (BOUN_TYPE )
        N = ( params2d.cells_number_x - 2 ) * ( params2d.cells_number_y - 2 ) * ( params2d.cells_number_z - 2 );
    else
        N = ( params2d.cells_number_x ) * ( params2d.cells_number_y ) * ( params2d.cells_number_z );
    solver_id = 0;
    vis = 0;
    print_system = 0;

    int arg_index = 0;
    int print_usage = 0;

    local_size = N / num_procs;
    extra = N - local_size * num_procs;

    ilower = local_size * myid;
    ilower += my_min(myid, extra);

    iupper = local_size * (myid + 1);
    iupper += my_min(myid + 1, extra);
    iupper = iupper - 1;

    /* How many rows do I have? */
    local_size = iupper - ilower + 1;

    int nvalues = local_size;
    int *rows = (int*) calloc(nvalues, sizeof(int));
    double *values =  (double*) calloc(nvalues, sizeof(double));

    for (i = 0; i < nvalues; i++)
    {
        rows[i] = ilower + i;
    }

    /* get the local solution */
    HYPRE_IJVectorGetValues(*x, nvalues, rows, values);
    for ( i = ilower; i <= iupper; i ++) {
        p_prev_3D[i] = values[i];
    }
    free (values);
    free (rows);
}
