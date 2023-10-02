#include "poisson.h"

void fill_matrix(struct ParametersCommon paramsc, struct Parameters2d params2d, double **matrix, double *b, double ***u_first_step_3D, double ***v_first_step_3D, 
                     double ***w_first_step_3D, double *x, double *y, double *z, double dt ){

		double hx = x[1] - x[0];
		double hy = y[1] - y[0];
		double hz = z[1] - z[0];
        printf("\n");
                printf("\npoisson begin   \n" );
		//������� ��������� ��������
                //zyxcxyz//
                int left, right, up, down, close, far;
                double x_minus, x_plus, y_minus, y_right, z_left, z_right, center, x_rhs, y_rhs, z_rhs, b_x, b_y, b_z;
		int row_counter = 0; // ������� ����� ������� �������
                int c_x = params2d.cells_number_x;
                int c_y = params2d.cells_number_y;
                int c_z = params2d.cells_number_z;
                if (BOUN_TYPE == 0 ){
                for ( int k = 0; k < params2d.cells_number_z; k++ ){
		    for ( int j = 0; j <  params2d.cells_number_y; j++ ){
		        for ( int i = 0; i < params2d.cells_number_x; i++ ){
                        
                            if ( j == 0){
                                y_rhs =  - dt / hy / hy;
                                matrix[row_counter][5] = dt / hy / hy;
                                b_y = v_first_step_3D[i][j + 1][k] / hy;
                            }
                            else{
                                if ( j == params2d.cells_number_y - 1 ){
                                    y_rhs =  - dt / hy / hy;
                                    matrix[row_counter][1] = dt / hy / hy;
                                    b_y = - v_first_step_3D[i][j][k] / hy;
                                }
                                else{
                                    y_rhs =  - 2 * dt / hy / hy;
                                    matrix[row_counter][5] = dt / hy / hy;
                                    matrix[row_counter][1] = dt / hy / hy;
                                    b_y = ( v_first_step_3D[i][j + 1][k] - v_first_step_3D[i][j][k] ) / hy;
                                }                                
                            }



                            if ( i == 0 ){
                                x_rhs =  - dt / hx / hx;
                                matrix[row_counter][4] = dt / hx / hx;
                                b_x = u_first_step_3D[i + 1][j][k] / hx;
                            }
                            else{
                                if ( i == params2d.cells_number_x - 1 ){
                                    x_rhs =  - dt / hx / hx;
                                    matrix[row_counter][2] = dt / hx / hx;
                                    b_x = - u_first_step_3D[i][j][k] / hx;
                                }
                                else{
                                    x_rhs =  - 2 * dt / hx / hx;
                                    matrix[row_counter][4] = dt / hx / hx;
                                    matrix[row_counter][2] = dt / hx / hx;
                                    b_x = ( u_first_step_3D[i + 1][j][k] - u_first_step_3D[i][j][k] ) / hx;
                                }                                
                            }



                            if ( k == 0 ){
                                z_rhs =  - dt / hz / hz;
                                matrix[row_counter][6] = dt / hz / hz;
                                b_z = w_first_step_3D[i][j][k + 1] / hz;
                            }
                            else{
                                if( k == params2d.cells_number_z - 1 ){
                                    z_rhs =  - dt / hz / hz;
                                    matrix[row_counter][0] = dt / hz / hz;
                                    b_z = - w_first_step_3D[i][j][k] / hz;
                                }
                                else{
                                    z_rhs =  - 2 * dt / hz / hz;
                                    matrix[row_counter][0] = dt / hz / hz;
                                    matrix[row_counter][6] = dt / hz / hz;
                                    b_z = ( w_first_step_3D[i][j][k + 1] - w_first_step_3D[i][j][k] ) / hz;
                                }                                
                            }
                            matrix[row_counter][3] = x_rhs + y_rhs + z_rhs;
                            b[row_counter] = b_x + b_y + b_z;
                            
                            
										    
		        row_counter++;
                        }
		    }
	        }
                }
                else
                if ( BOUN_TYPE == 1 ){
                for ( int k = 0; k < params2d.cells_number_z; k++ ){
		    for ( int j = 0; j <  params2d.cells_number_y; j++ ){
		        for ( int i = 0; i < params2d.cells_number_x; i++ ){
                        int cntr = k * params2d.cells_number_x * params2d.cells_number_y + j * params2d.cells_number_x + i;
                            if ( j == 0){
                                y_rhs =  0;
                                matrix[cntr][5] = 0;
                                matrix[cntr][1] = 0;
                                b_y = ( v_first_step_3D[i][j + 1][k] - v_first_step_3D[i][j][k] ) / hy;
                            }
                            else{
                                if ( j == params2d.cells_number_y - 1 ){
                                    y_rhs =  - dt / hy / hy;
                                    matrix[cntr][1] = dt / hy / hy;
                                    b_y = 0;//- v_first_step_3D[i][j][k] / hy;
                                }
                                else{
                                    y_rhs =  - 2 * dt / hy / hy;
                                    matrix[cntr][5] = dt / hy / hy;
                                    matrix[cntr][1] = dt / hy / hy;
                                    b_y = ( v_first_step_3D[i][j + 1][k] - v_first_step_3D[i][j][k] ) / hy;
                                }                                
                            }



                            if ( i == 0 ){
                                x_rhs =  - dt / hx / hx;
                                matrix[cntr][4] = dt / hx / hx;
                                b_x = 0;// u_first_step_3D[i + 1][j][k] / hx;
                            }
                            else{
                                if ( i == params2d.cells_number_x - 1 ){
                                    x_rhs =  - dt / hx / hx;
                                    matrix[cntr][2] = dt / hx / hx;
                                    b_x = 0;//- u_first_step_3D[i][j][k] / hx;
                                }
                                else{
                                    x_rhs =  - 2 * dt / hx / hx;
                                    matrix[cntr][4] = dt / hx / hx;
                                    matrix[cntr][2] = dt / hx / hx;
                                    b_x = ( u_first_step_3D[i + 1][j][k] - u_first_step_3D[i][j][k] ) / hx;
                                }                                
                            }



                            if ( k == 0 ){
                                z_rhs =  - dt / hz / hz;
                                matrix[cntr][6] = dt / hz / hz;
                                b_z = 0;//w_first_step_3D[i][j][k + 1] / hz;
                            }
                            else{
                                if( k == params2d.cells_number_z - 1 ){
                                    z_rhs =  - dt / hz / hz;
                                    matrix[cntr][0] = dt / hz / hz;
                                    b_z = 0;//- w_first_step_3D[i][j][k] / hz;
                                }
                                else{
                                    z_rhs =  - 2 * dt / hz / hz;
                                    matrix[cntr][0] = dt / hz / hz;
                                    matrix[cntr][6] = dt / hz / hz;
                                    b_z = ( w_first_step_3D[i][j][k + 1] - w_first_step_3D[i][j][k] ) / hz;
                                }                                
                            }
                            matrix[cntr][3] = x_rhs + y_rhs + z_rhs;
                            b[cntr] = b_x + b_y + b_z;
                            
                            
										    
		        row_counter++;
                        }
		    }
	        }
                }
                printf("\nmatrix is filled\n");
}