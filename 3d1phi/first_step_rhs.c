#include "first_step_rhs.h"
#define my_min(a,b)  (((a)<(b)) ? (a) : (b))
void calc_first_step_rhs(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***p_next_3D, double ***u_inter_3D, double ***v_inter_3D, double ***w_inter_3D , double ***u_prev_3D, double ***v_prev_3D, double ***w_prev_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, double *x_grid, double *y, double *z, double dt, double time, int myid, int num_procs, HYPRE_IJVector *x,
					 HYPRE_ParVector *par_x, HYPRE_IJVector *b, HYPRE_ParVector *par_b, double ***x_nonlin, double ***y_nonlin, double ***z_nonlin, double ***x_force,
					 double ***y_force, double ***z_force){

    double hx = x_grid[1] - x_grid[0];
    double hy = y[1] - y[0];
    double hz = z[1] - z[0];

	int i;
	int N, n;

	int ilower, iupper;
	int local_size, extra;

	int solver_id;
	int vis, print_system;

    double mu = paramsc.mu / 2.0;
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
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, b);
    HYPRE_IJVectorSetObjectType(*b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(*b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, x);
    HYPRE_IJVectorSetObjectType(*x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(*x);


    {
        double *rhs_values, *x_values;
        int    *rows;

        rhs_values =  (double*) calloc(local_size, sizeof(double));
        x_values =  (double*) calloc(local_size, sizeof(double));
        rows = (int*) calloc(local_size, sizeof(int));
            for (i = ilower; i <= iupper; i++) {
            	x_values[i] = 0.0;
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

                if (eq == 0){
                	if (i_center == 0 || i_center == cx)
                		rhs_values[i] = 0.0;
                	else if ( j_center == 0 )
                		rhs_values[i] = u_prev_3D[i_center][0][k_center] / 2.0 + u_prev_3D[i_center][1][k_center] / 6.0;
                	else if ( j_center == cy - 1 )
                		rhs_values[i] = u_prev_3D[i_center][cy - 1][k_center] / 2.0 + u_prev_3D[i_center][cy - 2][k_center] / 6.0;
                	else if ( k_center == 0 )
                		rhs_values[i] = u_prev_3D[i_center][j_center][0] / 2.0 + u_prev_3D[i_center][j_center][1] / 6.0;
                	else if ( k_center == cz - 1 )
                		rhs_values[i] = u_prev_3D[i_center][j_center][cz - 1] / 2.0 + u_prev_3D[i_center][j_center][cz - 2] / 6.0;
                	else{

		double nonlin_n_1 = x_nonlin[i_center][j_center][k_center];;//calc_x_nonlinear_term( u_inter_3D, v_inter_3D, w_inter_3D, i_center, j_center, k_center, hx, hy, hz );
		double nonlin_n = calc_x_nonlinear_term( u_prev_3D, v_prev_3D, w_prev_3D, i_center, j_center, k_center, hx, hy, hz );
		x_nonlin[i_center][j_center][k_center] = nonlin_n;
		double visc = calc_x_visc_term( paramsc, u_prev_3D, v_prev_3D, w_prev_3D, i_center, j_center, k_center, hx, hy, hz );
                		rhs_values[i] =  u_prev_3D[i_center][j_center][k_center] +
		    dt * ( 3.0 * nonlin_n - nonlin_n_1 ) / 2.0  + dt * visc / 2.0 + 2 * paramsc.omega * dt * (v_prev_3D[i_center][j_center][k_center] +
		    v_prev_3D[i_center - 1][j_center][k_center] + v_prev_3D[i_center][j_center + 1][k_center] + v_prev_3D[i_center - 1][j_center + 1][k_center] ) / 4.0 +
			x_force[i_center][j_center][k_center];
                		x_values[i] = u_prev_3D[i_center][j_center][k_center];
                	}
                }
                else
                if (eq == 1){
                	if (j_center == 0 || j_center == cy)
                		rhs_values[i] = 0.0;
                	else if ( i_center == 0 )
                		rhs_values[i] = v_prev_3D[0][j_center][k_center] / 2.0 + u_prev_3D[1][j_center][k_center] / 6.0;
                	else if ( i_center == cx - 1 )
                		rhs_values[i] = v_prev_3D[cx - 1][j_center][k_center] / 2.0 + u_prev_3D[cx - 2][j_center][k_center] / 6.0;
                	else if ( k_center == 0 )
                		rhs_values[i] = v_prev_3D[i_center][j_center][0] / 2.0 + v_prev_3D[i_center][j_center][1] / 6.0;
                	else if ( k_center == cz - 1 )
                		rhs_values[i] = v_prev_3D[i_center][j_center][cz - 1] / 2.0 + v_prev_3D[i_center][j_center][cz - 2] / 6.0;
                	else{

		double nonlin_n_1 = y_nonlin[i_center][j_center][k_center];//calc_y_nonlinear_term( u_inter_3D, v_inter_3D, w_inter_3D, i_center, j_center, k_center, hx, hy, hz );
		double nonlin_n = calc_y_nonlinear_term( u_prev_3D, v_prev_3D, w_prev_3D, i_center, j_center, k_center, hx, hy, hz );
		y_nonlin[i_center][j_center][k_center] = nonlin_n;
		double visc = calc_y_visc_term( paramsc, u_prev_3D, v_prev_3D, w_prev_3D, i_center, j_center, k_center, hx, hy, hz );
                		rhs_values[i] = v_prev_3D[i_center][j_center][k_center] +
		    dt * ( 3.0 * nonlin_n - nonlin_n_1 ) / 2.0  + dt * visc / 2.0 - 2 * paramsc.omega * dt * (u_prev_3D[i_center][j_center][k_center] +
		    u_prev_3D[i_center + 1][j_center][k_center] + u_prev_3D[i_center][j_center - 1][k_center] + u_prev_3D[i_center + 1][j_center - 1][k_center] ) / 4.0 +
			y_force[i_center][j_center][k_center];
                		x_values[i] = v_prev_3D[i_center][j_center][k_center];
                	}
                }
                else
                if (eq == 2){
                	if (k_center == 0 || k_center == cz)
            	rhs_values[i] = 0.0;
                    else if ( i_center == 0 )
                    	rhs_values[i] = w_prev_3D[0][j_center][k_center] / 2.0 + w_prev_3D[1][j_center][k_center] / 6.0;
                    else if ( i_center == cx - 1 )
                    		rhs_values[i] = w_prev_3D[cx - 1][j_center][k_center] / 2.0 + w_prev_3D[cx - 2][j_center][k_center] / 6.0;
            	else if ( j_center == 0 )
                    	rhs_values[i] = w_prev_3D[i_center][0][k_center] / 2.0 + w_prev_3D[i_center][1][k_center] / 6.0;
                    else if ( j_center == cy - 1 )
                    	rhs_values[i] = w_prev_3D[i_center][cy - 1][k_center] / 2.0 + w_prev_3D[i_center][cy - 2][k_center] / 6.0;
                    else{

		double nonlin_n_1 = z_nonlin[i_center][j_center][k_center];;//calc_z_nonlinear_term( u_inter_3D, v_inter_3D, w_inter_3D, i_center, j_center, k_center, hx, hy, hz );
		double nonlin_n = calc_z_nonlinear_term( u_prev_3D, v_prev_3D, w_prev_3D, i_center, j_center, k_center, hx, hy, hz );
		z_nonlin[i_center][j_center][k_center] = nonlin_n;
		double visc = calc_z_visc_term( paramsc, u_prev_3D, v_prev_3D, w_prev_3D, i_center, j_center, k_center, hx, hy, hz );
                		rhs_values[i] = w_prev_3D[i_center][j_center][k_center] +  dt * ( 3.0 * nonlin_n - nonlin_n_1 ) / 2.0 + dt * visc / 2.0 +
                		z_force[i_center][j_center][k_center];
                    	x_values[i] = w_prev_3D[i_center][j_center][k_center];
                    }
                    }


                rows[i] = ilower + i;
            }
            HYPRE_IJVectorSetValues(*b, local_size, rows, rhs_values);
            HYPRE_IJVectorSetValues(*x, local_size, rows, x_values);
            free(x_values);
            free(rhs_values);
            free(rows);
    }
    printf("\n rhs is written, begin assembly\n");
    HYPRE_IJVectorAssemble(*b);

    HYPRE_IJVectorGetObject(*b, (void **) par_b);

    HYPRE_IJVectorAssemble(*x);
    HYPRE_IJVectorGetObject(*x, (void **) par_x);

}

void init_nonlin( double ***x_nonlin, double ***y_nonlin, double ***z_nonlin, int nx, int ny, int nz){
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int k = 0; k < nz; k++){
				x_nonlin[i][j][k] = 0.0;
				y_nonlin[i][j][k] = 0.0;
				z_nonlin[i][j][k] = 0.0;
			}

}



void init_force(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***x_force, double ***y_force, double ***z_force, double N, double *x_grid,
		double *y, double *z, double dt){
    double f0 = paramsc.f0;
    double Ax = paramsc.mu * paramsc.mu * 10000 * 4.3*4.3*4.3;
    double Ay = Ax;
    double Az = - 2.0 * Ax;
    double kx = 1.5;
    double ky = 2.5;
    double kz = 4.5;
    double kx2 = 2.5;
    double ky2 = 1.5;
    double kz2 = 3.5;
    double kx3 = 1.5;
    double ky3 = 1.5;
    double kz3 = 3.5;
	int cx = params2d.cells_number_x;
	int cy = params2d.cells_number_y;
	int cz = params2d.cells_number_z;
	double tmp = 0.0;
    for (int i = 0; i < N; i++) {
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
             if (eq == 0){
             	if (i_center == 0 || i_center == cx)
             		tmp = 0.0;
             	else if ( j_center == 0 )
             		tmp = 0.0;
             	else if ( j_center == cy - 1 )
             		tmp = 0.0;
             	else if ( k_center == 0 )
             		tmp = 0.0;
             	else if ( k_center == cz - 1 )
             		tmp = 0.0;
             	else{

             		x_force[i_center][j_center][k_center] = f0 * dt * ( Ax / kx * cos(kx * x_grid[i_center] ) * sin( ky * ( y[j_center] + y[j_center+1] ) / 2.0 ) * sin( kz * ( z[k_center] + z[k_center+1] ) / 2.0 ) +
				Ax / kx2 * cos(kx2 * x_grid[i_center] ) * sin( ky2 * ( y[j_center] + y[j_center+1] ) / 2.0 ) * sin( kz2 * ( z[k_center] + z[k_center+1] ) / 2.0 ) +
				Ax / kx3 * cos(kx3 * x_grid[i_center] ) * sin( ky3 * ( y[j_center] + y[j_center+1] ) / 2.0 ) * sin( kz3 * ( z[k_center] + z[k_center+1] ) / 2.0 ) );

             	}
             }
             else
             if (eq == 1){
             	if (j_center == 0 || j_center == cy)
             		tmp = 0.0;
             	else if ( i_center == 0 )
             		tmp = 0.0;
             	else if ( i_center == cx - 1 )
             		tmp = 0.0;
             	else if ( k_center == 0 )
             		tmp = 0.0;
             	else if ( k_center == cz - 1 )
             		tmp = 0.0;
             	else{

             		y_force[i_center][j_center][k_center] = f0 * dt * ( Ay / ky * sin(kx * ( x_grid[i_center + 1] + x_grid[i_center] ) / 2.0 ) * cos( ky * (y[j_center] ) ) *
				sin( kz * ( z[k_center] + z[k_center+1] ) / 2.0 ) + Ay / ky2 * sin(kx2 * ( x_grid[i_center + 1] + x_grid[i_center] ) / 2.0 ) * cos( ky2 * (y[j_center] ) ) *
				sin( kz2 * ( z[k_center] + z[k_center+1] ) / 2.0 ) + Ay / ky3 * sin(kx3 * ( x_grid[i_center + 1] + x_grid[i_center] ) / 2.0 ) * cos( ky3 * (y[j_center] ) ) *
				sin( kz3 * ( z[k_center] + z[k_center+1] ) / 2.0 ) );
             	}
             }
             else
             if (eq == 2){
             	if (k_center == 0 || k_center == cz)
             		tmp = 0.0;
                 else if ( i_center == 0 )
                	 tmp = 0.0;
                 else if ( i_center == cx - 1 )
                	 tmp = 0.0;
                 else if ( j_center == 0 )
                	 tmp = 0.0;
                 else if ( j_center == cy - 1 )
                	 tmp = 0.0;
                 else{

                	 z_force[i_center][j_center][k_center] = f0 * dt * ( Az / kz * sin(kx * ( x_grid[i_center + 1] + x_grid[i_center] ) / 2.0 ) *
				sin( ky * ( y[j_center] + y[j_center+1] ) / 2.0 ) * cos( kz * ( z[k_center] ) ) +  Az / kz2 * sin(kx2 * ( x_grid[i_center + 1] + x_grid[i_center] ) / 2.0 ) *
				sin( ky2 * ( y[j_center] + y[j_center+1] ) / 2.0 ) * cos( kz2 * ( z[k_center] ) ) +  Az / kz3 * sin(kx3 * ( x_grid[i_center + 1] + x_grid[i_center] ) / 2.0 ) *
				sin( ky3 * ( y[j_center] + y[j_center+1] ) / 2.0 ) * cos( kz3 * ( z[k_center] ) ) );
                 }
                 }

         }

}


