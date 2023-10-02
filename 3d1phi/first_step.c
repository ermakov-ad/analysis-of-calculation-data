#include "first_step.h"
#define my_min(a,b)  (((a)<(b)) ? (a) : (b))
void calc_first_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***p_next_3D, double ***u_inter_3D, double ***v_inter_3D, double ***w_inter_3D , double ***u_prev_3D, double ***v_prev_3D, double ***w_prev_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, double *x_grid, double *y, double *z, double dt, double time, int myid, int num_procs, HYPRE_IJVector *x,
					 HYPRE_ParVector *par_x, HYPRE_IJVector *b, HYPRE_ParVector *par_b, HYPRE_IJMatrix *AA, HYPRE_ParCSRMatrix *parcsr_A){

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
              HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, AA);

                /* Choose a parallel csr format storage (see the User's Manual) */
                HYPRE_IJMatrixSetObjectType(*AA, HYPRE_PARCSR);
                /* Initialize before setting coefficients */
                HYPRE_IJMatrixInitialize(*AA);
                {
                    int nnz;

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
                    	        values[nnz] = 1.0 + 2.0 * dt * mu * ( 1.0 / hx / hx + 1.0 / hy / hy + 1.0 / hz / hz );
                    	        nnz++;

                    	        cols[nnz] = i_center + 1 + ( cx + 1 ) * ( j_center ) + ( cx + 1 ) * cy * k_center;
                    	        values[nnz] = - dt * mu / hx / hx;
                    	        nnz++;

                    	        cols[nnz] = i_center - 1 + ( cx + 1 ) * ( j_center ) + ( cx + 1 ) * cy * k_center;
                    	        values[nnz] = - dt * mu / hx / hx;
                    	        nnz++;

                    	        cols[nnz] = i_center + ( cx + 1 ) * ( j_center + 1 ) + ( cx + 1 ) * cy * k_center;
                    	        values[nnz] = - dt * mu / hy / hy;
                    	        nnz++;

                    	        cols[nnz] = i_center + ( cx + 1 ) * ( j_center - 1 ) + ( cx + 1 ) * cy * k_center;
                    	        values[nnz] = - dt * mu / hy / hy;
                    	        nnz++;

                    	        cols[nnz] = i_center + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * ( k_center + 1 );
                    	        values[nnz] = - dt * mu / hz / hz;
                    	        nnz++;

                    	        cols[nnz] = i_center + ( cx + 1 ) * j_center + ( cx + 1 ) * cy * ( k_center - 1 );
                    	        values[nnz] = - dt * mu / hz / hz;
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
                    	        values[nnz] = 1.0 + 2.0 * dt * mu * ( 1.0 / hx / hx + 1.0 / hy / hy + 1.0 / hz / hz );
                    	        nnz++;

                    	        //V
                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * ( j_center + 1 ) + ( cy + 1 ) * cx * ( k_center );
                    	        values[nnz] = - dt * mu / hy / hy;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * ( j_center - 1 ) + ( cy + 1 ) * cx * ( k_center );
                    	        values[nnz] = - dt * mu / hy / hy;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * j_center + ( cy + 1 ) * cx * ( k_center - 1 );
                    	        values[nnz] = - dt * mu / hz / hz;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + cx * j_center + ( cy + 1 ) * cx * ( k_center + 1 );
                    	        values[nnz] = - dt * mu / hz / hz;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center + 1 + cx * j_center + ( cy + 1 ) * cx * k_center;
                    	        values[nnz] = - dt * mu / hx / hx;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz ) + i_center - 1 + cx * j_center + ( cy + 1 ) * cx * k_center;
                    	        values[nnz] = - dt * mu / hx / hx;
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
                    	        values[nnz] = 1.0 + 2.0 * dt * mu * ( 1.0 / hx / hx + 1.0 / hy / hy + 1.0 / hz / hz );
                    	        nnz++;

                    	        //W
                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + 1 + cx * j_center + cy * cx * k_center ;
                    	        values[nnz] = - dt * mu / hx / hx;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center - 1 + cx * j_center + cy * cx * k_center ;
                    	        values[nnz] = - dt * mu / hx / hx;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center + 1 ) + cy * cx * k_center ;
                    	        values[nnz] = - dt * mu / hy / hy;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center - 1 ) + cy * cx * k_center ;
                    	        values[nnz] = - dt * mu / hy / hy;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center ) + cy * cx * ( k_center + 1 ) ;
                    	        values[nnz] = - dt * mu / hz / hz;
                    	        nnz++;

                    	        cols[nnz] = ( ( cx + 1 ) * cy * cz  + cx * ( cy + 1 ) * cz ) + i_center + cx * ( j_center ) + cy * cx * ( k_center - 1 ) ;
                    	        values[nnz] = - dt * mu / hz / hz;
                    	        nnz++;

                    	    }

                    	}
                        tmp[0] = nnz;
                        tmp[1] = i;
                        HYPRE_IJMatrixSetValues(*AA, 1, &tmp[0], &tmp[1], cols, values);
                    }
                    free(values);
                    free(cols);
                    free(tmp);
                }
                printf("\n matrix is written, begin assembly\n");
                //printf("\ncheck\n");
                /* Assemble after setting the coefficients */
                HYPRE_IJMatrixAssemble(*AA);

                /* Get the parcsr matrix object to use */
                HYPRE_IJMatrixGetObject(*AA, (void**) parcsr_A);

}

double calc_x_nonlinear_term(double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz ){


   return  ( ( ( u[i][j][k] + u[i - 1][j][k] ) / 2.0) * ( ( u[i][j][k] + u[i - 1][j][k] ) / 2.0) -
    ( ( u[i + 1][j][k] + u[i][j][k] ) / 2.0 ) * ( ( u[i + 1][j][k] + u[i][j][k] ) / 2.0 ) ) / hx +
	( ( w[i][j][k] + w[i - 1][j][k] ) / 2.0 * ( u[i][j][k] + u[i][j][k - 1] ) / 2.0 -
	( w[i][j][k + 1] + w[i - 1][j][k + 1] ) / 2.0 * ( u[i][j][k + 1] + u[i][j][k] ) / 2.0 ) / hz +
	( ( v[i][j][k] + v[i - 1][j][k] ) / 2.0 * ( u[i][j][k] + u[i][j - 1][k] ) / 2.0 - 
	( v[i][j + 1][k] + v[i - 1][j + 1][k] ) / 2.0 * ( u[i][j + 1][k] + u[i][j][k] ) / 2.0 ) / hy;
}



double calc_y_nonlinear_term(double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz ){


    double dv2 = ( pow( ( v[i][j][k] + v[i][j - 1][k] ) / 2.0, 2.0 ) - pow( ( v[i][j + 1][k] + v[i][j][k] ) / 2.0, 2.0 ) ) / hy;
    return dv2 + ( ( v[i][j][k] + v[i - 1][j][k] ) / 2.0 * ( u[i][j - 1][k] + u[i][j][k] ) / 2.0 - 
	( v[i + 1][j][k] + v[i][j][k] ) / 2.0 * ( u[i + 1][j - 1][k] + u[i + 1][j][k] ) / 2.0 ) / hx +
	( ( w[i][j - 1][k] + w[i][j][k] ) / 2.0 * ( v[i][j][k] + v[i][j][k - 1] ) / 2.0 - 
	( w[i][j - 1][k + 1] + w[i][j][k + 1] ) / 2.0 * ( v[i][j][k + 1] + v[i][j][k] ) / 2.0 ) / hz;
                                        
}



double calc_z_nonlinear_term(double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz ){


    double dw2 = ( pow( ( w[i][j][k] + w[i][j][k-1] ) / 2.0, 2.0 ) - pow( ( w[i][j][k+1] + w[i][j][k] ) / 2.0, 2.0 ) ) / hz;
    return dw2 + ( ( u[i][j][k] + u[i][j][k - 1] ) / 2.0 * ( w[i][j][k] + w[i - 1][j][k] ) / 2.0 - 
	( u[i + 1][j][k] + u[i + 1][j][k - 1] ) / 2.0 * ( w[i + 1][j][k] + w[i][j][k] ) / 2.0 ) / hx +
	( ( v[i][j][k] + v[i][j][k - 1] ) / 2.0 * ( w[i][j][k] + w[i][j - 1][k] ) / 2.0 -
	( v[i][j + 1][k] + v[i][j + 1][k - 1] ) / 2.0 * ( w[i][j + 1][k] + w[i][j][k] ) / 2.0 ) / hy;
                                            
}

double calc_x_visc_term(struct ParametersCommon paramsc, double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz ){
    return paramsc.mu * ( ( u[i+1][j][k] - 2.0 * u[i][j][k] +  u[i-1][j][k] ) / hx / hx + 
	( u[i][j + 1][k] - 2.0 * u[i][j][k] + u[i][j - 1][k] ) / hy / hy +
        ( u[i][j][k + 1] - 2.0 * u[i][j][k] + u[i][j][k - 1] ) / hz / hz );
}

double calc_y_visc_term( struct ParametersCommon paramsc, double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz ){

    return paramsc.mu * ( ( v[i+1][j][k] - 2.0 * v[i][j][k] +  v[i-1][j][k] ) / hx / hx + 
	( v[i][j + 1][k] - 2.0 * v[i][j][k] + v[i][j - 1][k] ) / hy / hy +
        ( v[i][j][k + 1] - 2.0 * v[i][j][k] + v[i][j][k - 1] ) / hz / hz );
}
                                
double calc_z_visc_term(struct ParametersCommon paramsc, double ***u, double ***v, double ***w, int i, int j, int k, double hx, double hy, double hz ){

    return paramsc.mu * ( ( w[i+1][j][k] - 2.0 * w[i][j][k] +  w[i-1][j][k] ) / hx / hx +
	( w[i][j + 1][k] - 2.0 * w[i][j][k] + w[i][j - 1][k] ) / hy / hy +
        ( w[i][j][k + 1] - 2.0 * w[i][j][k] + w[i][j][k - 1] ) / hz / hz );
}

void fin_first_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***p_next_3D, double ***u_first_step_3D, double ***v_first_step_3D,
                     double ***w_first_step_3D, int myid, int num_procs, HYPRE_IJVector *x){
    int i;
    int N, n;

    int ilower, iupper;
    int local_size, extra;

    int solver_id;
    int vis, print_system;

#if defined(HYPRE_USING_GPU)
    /* use cuSPARSE for SpGEMM */
   HYPRE_SetSpGemmUseCusparse(0);
#endif
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
    int nvalues = local_size;
    int *rows = (int*) calloc(nvalues, sizeof(int));
    double *values =  (double*) calloc(nvalues, sizeof(double));

    for (int i = 0; i < nvalues; i++)
    {
        rows[i] = ilower + i;
    }

    /* get the local solution */
    HYPRE_IJVectorGetValues(*x, nvalues, rows, values);
    for ( int i = ilower; i <= iupper; i ++) {
		int k_center;
		int remainder_k;
		int j_center;
		int i_center;
		int eq;
        if ( i < (cx + 1)*cy*cz){
        	eq = 1;
			k_center = (int) i / ((params2d.cells_number_x + 1 ) * (params2d.cells_number_y ));
			remainder_k = (int) i % ((params2d.cells_number_x + 1) * (params2d.cells_number_y));
			j_center = remainder_k / (params2d.cells_number_x + 1);
			i_center = remainder_k % (params2d.cells_number_x + 1);
			u_first_step_3D[i_center][j_center][k_center] = values[i];
        }
        else if ( i < (cx + 1)*cy*cz + cx*(cy+1)*cz){
        	eq = 2;
			k_center = (int) ( i - (cx + 1)*cy*cz ) / ((params2d.cells_number_x ) * (params2d.cells_number_y + 1 ));
			remainder_k = (int) ( i - (cx + 1)*cy*cz ) % ((params2d.cells_number_x ) * (params2d.cells_number_y + 1 ));
			j_center = remainder_k / (params2d.cells_number_x );
			i_center = remainder_k % (params2d.cells_number_x );
			v_first_step_3D[i_center][j_center][k_center] = values[i];
        }
        else{
        	eq = 3;
			k_center = (int) ( i - (cx + 1)*cy*cz - cx*(cy+1)*cz ) / ((params2d.cells_number_x ) * (params2d.cells_number_y ));
			remainder_k = (int) ( i - (cx + 1)*cy*cz - cx*(cy+1)*cz ) % ((params2d.cells_number_x ) * (params2d.cells_number_y ));
			j_center = remainder_k / (params2d.cells_number_x );
			i_center = remainder_k % (params2d.cells_number_x );
			w_first_step_3D[i_center][j_center][k_center] = values[i];
        }
    }
    free (values);
    free (rows);
}
