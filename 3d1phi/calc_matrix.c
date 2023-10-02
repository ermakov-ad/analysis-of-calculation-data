#include "first_step.h"
#include "second_step.h"
#define my_min(a,b)  (((a)<(b)) ? (a) : (b))
HYPRE_ParCSRMatrix calc_first_step_matrix(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_inter_3D, double ***v_inter_3D, double ***w_inter_3D , double ***u_prev_3D, double ***v_prev_3D, double ***w_prev_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, double *x_grid, double *y, double *z, double dt, double time, int myid, int num_procs){
	double hx = x_grid[1] - x_grid[0];
	    double hy = y[1] - y[0];
	    double hz = z[1] - z[0];

	                int i;
	                int N, n;

	                int ilower, iupper;
	                int local_size, extra;

	                int solver_id;
	                int vis, print_system;

	                HYPRE_IJMatrix AA;
	                HYPRE_ParCSRMatrix parcsr_A;
	                HYPRE_Solver solver, precond;


	                /* Initialize HYPRE */

	            #if defined(HYPRE_USING_GPU)
	                /* use cuSPARSE for SpGEMM */
	               HYPRE_SetSpGemmUseCusparse(0);
	            #endif

	               	N = ( params2d.cells_number_x + 1 ) * ( params2d.cells_number_y) * ( params2d.cells_number_z ) +
	               		( params2d.cells_number_x ) * ( params2d.cells_number_y + 1) * ( params2d.cells_number_z ) +
				( params2d.cells_number_x ) * ( params2d.cells_number_y ) * ( params2d.cells_number_z + 1 );
	               	int cx = params2d.cells_number_x;
	               	int cy = params2d.cells_number_y;
	               	int cz = params2d.cells_number_z;
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

	                /* Create the matrix.
	                   Note that this is a square matrix, so we indicate the row partition
	                   size twice (since number of rows = number of cols) */
	                HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &AA);

	                /* Choose a parallel csr format storage (see the User's Manual) */
	                HYPRE_IJMatrixSetObjectType(AA, HYPRE_PARCSR);

	                /* Initialize before setting coefficients */
	                HYPRE_IJMatrixInitialize(AA);
	                {
	                    int nnz;
	                    /* OK to use constant-length arrays for CPUs
	                    double values[5];
	                    int cols[5];
	                    */

	                    double *values = (double *) malloc(13 * sizeof(double));
	                    int *cols = (int *) malloc(13 * sizeof(int));
	                    int *tmp = (int *) malloc(2 * sizeof(int));
	                    for (i = ilower; i <= iupper; i++) {
	                    	nnz = 0;
	                    	int k_center;
	                    	int remainder_k;
	                    	int j_center;
	                    	int i_center;
	                    	int eq;
	                    	if ( i < (cx + 1)*cy*cz){
	                    		k_center = (int) i / ((params2d.cells_number_x + 1 ) * (params2d.cells_number_y ));
	                    		remainder_k = (int) i % ((params2d.cells_number_x + 1) * (params2d.cells_number_y));
	                    		j_center = remainder_k / (params2d.cells_number_x + 1);
	                    		i_center = remainder_k % (params2d.cells_number_x + 1);
	                    		eq = 0;
	                    	}
	                    	else if ( i < (cx + 1)*cy*cz + cx*(cy+1)*cz){
	                    		k_center = (int) ( i - (cx + 1)*cy*cz ) / ((params2d.cells_number_x ) * (params2d.cells_number_y + 1 ));
	                    		remainder_k = (int) ( i - (cx + 1)*cy*cz ) % ((params2d.cells_number_x ) * (params2d.cells_number_y + 1 ));
	                    		j_center = remainder_k / (params2d.cells_number_x );
	                    		i_center = remainder_k % (params2d.cells_number_x );
	                    		eq = 1;
	                    	}
	                    	else{
	                    		k_center = (int) ( i - (cx + 1)*cy*cz - cx*(cy+1)*cz ) / ((params2d.cells_number_x ) * (params2d.cells_number_y ));
	                    		remainder_k = (int) ( i - (cx + 1)*cy*cz - cx*(cy+1)*cz ) % ((params2d.cells_number_x ) * (params2d.cells_number_y ));
	                    		j_center = remainder_k / (params2d.cells_number_x );
	                    		i_center = remainder_k % (params2d.cells_number_x );
	                    		eq = 2;
	                    	}
	                    	//if (eq == 1)
	                    	//printf("\n %d %d %d %d\n", i_center, j_center, k_center, eq);
	                    	if (eq == 0){
	                    		if (i_center == 0){
	                    			cols[nnz] = i;
	                    			values[nnz] = 1.0;
	                    			nnz++;
	                    		}else
	                    		if (i_center == cx){
	                    			cols[nnz] = i;
	                    			values[nnz] = 1.0;
	                    			nnz++;
	                    		}else
	                    	    if (j_center == 0){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (j_center == cy-1){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (k_center == 0){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (k_center == cz-1){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else{
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0 + 2.0 * dt * paramsc.mu / hy / hy + 2.0 * dt * paramsc.mu / hz / hz;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * ( j_center + 1 ) + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = - dt * paramsc.mu / hy / hy;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * ( j_center - 1 ) + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = - dt * paramsc.mu / hy / hy;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * ( k_center + 1 );
	                    	        values[nnz] = - dt * paramsc.mu / hz / hz;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * ( k_center - 1 );
	                    	        values[nnz] = - dt * paramsc.mu / hz / hz;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * ( j_center + 1 ) + ( cy + 1 ) * cx * k_center;
	                    	        values[nnz] = dt * paramsc.mu / hx / hy - dt * paramsc.omega / 2.0;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center - 1 + cx * ( j_center + 1 ) + ( cy + 1 ) * cx * k_center;
	                    	        values[nnz] = - dt * paramsc.mu / hx / hy - dt * paramsc.omega / 2.0;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * j_center + ( cy + 1 ) * cx * k_center;
	                    	        values[nnz] = -dt * paramsc.mu / hx / hy - dt * paramsc.omega / 2.0;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center - 1 + cx * j_center + ( cy + 1 ) * cx * k_center;
	                    	        values[nnz] = dt * paramsc.mu / hx / hy - dt * paramsc.omega / 2.0;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * j_center + cy * cx * ( k_center + 1 );
	                    	        values[nnz] = dt * paramsc.mu / hx / hy;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center - 1 + cx * j_center + cy * cx * ( k_center + 1 );
	                    	        values[nnz] = - dt * paramsc.mu / hx / hy;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * j_center + cy * cx * k_center;
	                    	        values[nnz] = - dt * paramsc.mu / hx / hy;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center - 1 + cx * j_center + cy * cx * k_center ;
	                    	        values[nnz] = dt * paramsc.mu / hx / hy;
	                    	        nnz++;
	                    	    }

	                    	}

	                    	if (eq == 1){
	                    		if (i_center == 0){
	                    			cols[nnz] = i;
	                    			values[nnz] = 1.0;
	                    			nnz++;
	                    		}else
	                    		if (i_center == cx - 1){
	                    			cols[nnz] = i;
	                    			values[nnz] = 1.0;
	                    			nnz++;
	                    		}else
	                    	    if (j_center == 0){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (j_center == cy){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (k_center == 0){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (k_center == cz-1){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else{
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0 + 2.0 * dt * paramsc.mu / hx / hx + 2.0 * dt * paramsc.mu / hz / hz;
	                    	        nnz++;

	                    	        //V
	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * j_center + ( cy + 1 ) * cx * ( k_center - 1 );
	                    	        values[nnz] = - dt * paramsc.mu / hz / hz;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * j_center + ( cy + 1 ) * cx * ( k_center + 1 );
	                    	        values[nnz] = - dt * paramsc.mu / hz / hz;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + 1 + cx * j_center + ( cy + 1 ) * cx * k_center;
	                    	        values[nnz] = - dt * paramsc.mu / hx / hx;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center - 1 + cx * j_center + ( cy + 1 ) * cx * k_center;
	                    	        values[nnz] = - dt * paramsc.mu / hx / hx;
	                    	        nnz++;

	                    	        //W
	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * j_center + cy * cx * ( k_center + 1 ) ;
	                    	        values[nnz] = dt * paramsc.mu / hy / hz;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center - 1 ) + cy * cx * ( k_center + 1 ) ;
	                    	        values[nnz] = - dt * paramsc.mu / hy / hz;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * j_center + cy * cx * k_center ;
	                    	        values[nnz] = - dt * paramsc.mu / hy / hz;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center - 1 ) + cy * cx * k_center ;
	                    	        values[nnz] = dt * paramsc.mu / hy / hz;
	                    	        nnz++;

	                    	        //U
	                    	        cols[nnz] = i_center + 1 + ( cx + 1 ) * ( j_center + 1 ) + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = dt * paramsc.mu / hx / hy + dt * paramsc.omega / 2.0;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + 1 + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = -dt * paramsc.mu / hx / hy + dt * paramsc.omega / 2.0;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * ( j_center + 1 ) + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = -dt * paramsc.mu / hx / hy + dt * paramsc.omega / 2.0;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = dt * paramsc.mu / hx / hy + dt * paramsc.omega / 2.0;
	                    	        nnz++;


	                    	    }

	                    	}

	                    	if (eq == 2){
	                    		if (i_center == 0){
	                    			cols[nnz] = i;
	                    			values[nnz] = 1.0;
	                    			nnz++;
	                    		}else
	                    		if (i_center == cx - 1){
	                    			cols[nnz] = i;
	                    			values[nnz] = 1.0;
	                    			nnz++;
	                    		}else
	                    	    if (j_center == 0){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (j_center == cy - 1){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (k_center == 0){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else
	                    	    if (k_center == cz){
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0;
	                    	        nnz++;
	                    	    }else{
	                    	        cols[nnz] = i;
	                    	        values[nnz] = 1.0 + 2.0 * dt * paramsc.mu / hy / hy + 2.0 * dt * paramsc.mu / hx / hx;
	                    	        nnz++;

	                    	        //W
	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + 1 + cx * j_center + cy * cx * k_center ;
	                    	        values[nnz] = - dt * paramsc.mu / hx / hx;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center - 1 + cx * j_center + cy * cx * k_center ;
	                    	        values[nnz] = - dt * paramsc.mu / hx / hx;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center + 1 ) + cy * cx * k_center ;
	                    	        values[nnz] = - dt * paramsc.mu / hy / hy;
	                    	        nnz++;

	                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center - 1 ) + cy * cx * k_center ;
	                    	        values[nnz] = - dt * paramsc.mu / hy / hy;
	                    	        nnz++;
	                    	        //U
	                    	        cols[nnz] = i_center + 1 + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = dt * paramsc.mu / hx / hz;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + 1 + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * ( k_center - 1 );
	                    	        values[nnz] = - dt * paramsc.mu / hx / hz;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * k_center;
	                    	        values[nnz] = - dt * paramsc.mu / hx / hz;
	                    	        nnz++;

	                    	        cols[nnz] = i_center + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * ( k_center - 1 );
	                    	        values[nnz] = dt * paramsc.mu / hx / hz;
	                    	        nnz++;

	                    	    }

	                    	}
	                        tmp[0] = nnz;
	                        tmp[1] = i;
	                        HYPRE_IJMatrixSetValues(AA, 1, &tmp[0], &tmp[1], cols, values);
	                    }
	                    free(values);
	                    free(cols);
	                    free(tmp);
	                }
	                printf("\n matrix is written, begin assembly\n");
	                //printf("\ncheck\n");
	                /* Assemble after setting the coefficients */
	                HYPRE_IJMatrixAssemble(AA);
	                /* Get the parcsr matrix object to use */
	                HYPRE_IJMatrixGetObject(AA, (void**) &parcsr_A);
	                //HYPRE_IJMatrixDestroy(AA);
	                return parcsr_A;

}
HYPRE_ParCSRMatrix calc_second_step_matrix(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_first_step_3D, double ***v_first_step_3D,
                      double ***w_first_step_3D, double *p_prev_3D, double *x_grid, double *y_grid, double *z_grid, double dt, int myid, int num_procs){
    double hx = x_grid[1] - x_grid[0];
    double hy = y_grid[1] - y_grid[0];
    double hz = z_grid[1] - z_grid[0];
    int i;
    //int myid, num_procs;
    int N, n;

    int ilower, iupper;
    int local_size, extra;

    int solver_id;
    int vis, print_system;


    HYPRE_IJVector x;
    HYPRE_ParVector par_x;
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;

    HYPRE_IJMatrix AA;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_Solver solver, precond;


    /* Initialize HYPRE */

#if defined(HYPRE_USING_GPU)
    /* use cuSPARSE for SpGEMM */
   HYPRE_SetSpGemmUseCusparse(0);
#endif

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

    /* Create the matrix.
       Note that this is a square matrix, so we indicate the row partition
       size twice (since number of rows = number of cols) */
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &AA);

    /* Choose a parallel csr format storage (see the User's Manual) */
    HYPRE_IJMatrixSetObjectType(AA, HYPRE_PARCSR);

    /* Initialize before setting coefficients */
    HYPRE_IJMatrixInitialize(AA);
    // printf("%lf\n", u_first_step_3D[1][1][1]);
    /* Now go through my local rows and set the matrix entries.
       Each row has at most 5 entries. For example, if n=3:

       A = [M -I 0; -I M -I; 0 -I M]
       M = [4 -1 0; -1 4 -1; 0 -1 4]

       Note that here we are setting one row at a time, though
       one could set all the rows together (see the User's Manual).
    */
    //double *b_values =  (double*) calloc(local_size, sizeof(double));
    {
        int nnz;
        /* OK to use constant-length arrays for CPUs
        double values[5];
        int cols[5];
        */

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
                HYPRE_IJMatrixSetValues(AA, 1, &tmp[0], &tmp[1], cols, values);
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
                HYPRE_IJMatrixSetValues(AA, 1, &tmp[0], &tmp[1], cols, values);
            }
        }

        free(values);
        free(cols);
        free(tmp);
    }
    //printf("\ncheck\n");
    /* Assemble after setting the coefficients */
    HYPRE_IJMatrixAssemble(AA);
    /* Get the parcsr matrix object to use */
    HYPRE_IJMatrixGetObject(AA, (void**) &parcsr_A);
    //HYPRE_IJMatrixDestroy(AA);
    return parcsr_A;
}
