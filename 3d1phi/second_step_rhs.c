#include"second_step_rhs.h"
#define my_min(a,b)  (((a)<(b)) ? (a) : (b))
void calc_second_step_rhs(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_first_step_3D, double ***v_first_step_3D,
                      double ***w_first_step_3D, double *p_prev_3D, double *x_grid, double *y_grid, double *z_grid, double dt, int myid, int num_procs,HYPRE_IJVector *x,
						 HYPRE_ParVector *par_x, HYPRE_IJVector *b, HYPRE_ParVector *par_b){


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

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, b);
    HYPRE_IJVectorSetObjectType(*b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(*b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, x);
    HYPRE_IJVectorSetObjectType(*x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(*x);
    /* Set the rhs values to h^2 and the solution to zero */
    {
        double *rhs_values, *x_values;
        int    *rows;

        rhs_values =  (double*) calloc(local_size, sizeof(double));
        x_values =  (double*) calloc(local_size, sizeof(double));
        rows = (int*) calloc(local_size, sizeof(int));
        if (BOUN_TYPE == 0){
            for (i = ilower; i <= iupper; i++) {
                int k_center = (int) i / ((params2d.cells_number_x) * (params2d.cells_number_y));
                int remainder_k = (int) i % ((params2d.cells_number_x) * (params2d.cells_number_y));
                int j_center = remainder_k / (params2d.cells_number_x);
                int i_center = remainder_k % (params2d.cells_number_x);

                double b_x, b_y, b_z;
                if (i_center == 0)
                    b_x = (u_first_step_3D[i_center + 1][j_center][k_center]) / hx;
                else
                if (i_center == params2d.cells_number_x - 1)
                    b_x = ( -u_first_step_3D[i_center][j_center][k_center]) / hx;
                else
                    b_x = (u_first_step_3D[i_center + 1][j_center][k_center] -
                           u_first_step_3D[i_center][j_center][k_center]) / hx;

                if (j_center == 0)
                    b_y = (v_first_step_3D[i_center][j_center + 1][k_center]) / hy;
                else
                if (j_center == params2d.cells_number_y - 1)
                    b_y = ( -v_first_step_3D[i_center][j_center][k_center]) / hy;
                else
                    b_y = (v_first_step_3D[i_center][j_center + 1][k_center] -
                           v_first_step_3D[i_center][j_center][k_center]) / hy;

                if (k_center == 0)
                    b_z = (w_first_step_3D[i_center][j_center][k_center + 1]) / hz;
                else
                if (k_center == params2d.cells_number_z - 1)
                    b_z = (-w_first_step_3D[i_center][j_center][k_center]) / hz;
                else
                    b_z = (w_first_step_3D[i_center][j_center][k_center + 1] -
                           w_first_step_3D[i_center][j_center][k_center]) / hz;
                if (k_center == 1 && j_center == 1 && i_center == 1)
                    rhs_values[i] = 0;
                else
                    rhs_values[i] = b_x + b_y + b_z;
                //rhs_values[i] = b_values[i];
                //printf("%lf %lf %lf %d %d %d \n",b_x, b_y, b_z, i_center, j_center, k_center);

                x_values[i] = p_prev_3D[i];//0.0;

                rows[i] = ilower + i;
            }
        }

        if (BOUN_TYPE == 1) {
            for (i = ilower; i <= iupper; i++) {
                int k_center = (int) i / ((params2d.cells_number_x - 2) * (params2d.cells_number_y - 2)) + 1;
                int remainder_k = (int) i % ((params2d.cells_number_x - 2) * (params2d.cells_number_y - 2));
                int j_center = remainder_k / (params2d.cells_number_x - 2) + 1;
                int i_center = remainder_k % (params2d.cells_number_x - 2) + 1;

                double b_x, b_y, b_z;
                if (j_center == 1)
                    b_y = (v_first_step_3D[i_center][j_center + 1][k_center] -
                        v_first_step_3D[i_center][params2d.cells_number_y - 1][k_center]) / hy;
                    else
                        b_y = (v_first_step_3D[i_center][j_center + 1][k_center] -
                               v_first_step_3D[i_center][j_center][k_center]) / hy;

                if (i_center == 1)
                    b_x = (u_first_step_3D[i_center + 1][j_center][k_center] -
                           u_first_step_3D[params2d.cells_number_x - 1][j_center][k_center]) / hx;
                else
                    b_x = (u_first_step_3D[i_center + 1][j_center][k_center] -
                           u_first_step_3D[i_center][j_center][k_center]) / hx;

                if (k_center == 1)
                    b_z = (w_first_step_3D[i_center][j_center][k_center + 1] -
                           w_first_step_3D[i_center][j_center][params2d.cells_number_z - 1]) / hz;
                else
                    b_z = (w_first_step_3D[i_center][j_center][k_center + 1] -
                           w_first_step_3D[i_center][j_center][k_center]) / hz;


                //if (i == 0){
                  //  printf("%e %e %e %e %e %e\n",v_first_step_3D[i_center][j_center + 1][k_center], v_first_step_3D[i_center][j_center][k_center],
                    //       u_first_step_3D[i_center + 1][j_center][k_center], u_first_step_3D[i_center][j_center][k_center], w_first_step_3D[i_center][j_center][k_center + 1],
                      //     w_first_step_3D[i_center][j_center][k_center]);
                //}
                rhs_values[i] =  b_x + b_y + b_z;
                if (k_center - 1 == ( params2d.cells_number_z - 2 ) / 2 && i_center - 1 == ( params2d.cells_number_x - 2 ) / 2 && j_center - 1 == ( params2d.cells_number_y - 2 ) / 2)
                    rhs_values[i] = 0.0;
                x_values[i] = 0.0;
                rows[i] = ilower + i;
            }
        }
        HYPRE_IJVectorSetValues(*b, local_size, rows, rhs_values);
        HYPRE_IJVectorSetValues(*x, local_size, rows, x_values);
        free(x_values);
        free(rhs_values);
        free(rows);
    }

    //printf(" ");
    HYPRE_IJVectorAssemble(*b);
    /*  As with the matrix, for testing purposes, one may wish to read in a rhs:
        HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD,
                                  HYPRE_PARCSR, &b );
        as an alternative to the
        following sequence of HYPRE_IJVectors calls:
        Create, SetObjectType, Initialize, SetValues, and Assemble
    */
    HYPRE_IJVectorGetObject(*b, (void **) par_b);

    HYPRE_IJVectorAssemble(*x);
    HYPRE_IJVectorGetObject(*x, (void **) par_x);
}
