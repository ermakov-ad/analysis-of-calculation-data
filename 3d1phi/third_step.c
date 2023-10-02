#include "third_step.h"

void calc_third_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_next_3D, double ***v_next_3D, double ***w_next_3D, double ***u_first_step_3D, double ***v_first_step_3D, double ***w_first_step_3D,
double ***p_next_3D, double hx, double hy, double hz, double dt ){
    printf(" ");
    if ( BOUN_TYPE == 0 ){
printf("\n");
	for ( int j = 0; j <  params2d.cells_number_y; j++ ){
	    for ( int i = 0; i < params2d.cells_number_x; i++ ){
                for ( int k = 0; k < params2d.cells_number_z; k++ ){
                    if ( j == 0){
                        v_next_3D[i][j][k] = 0.0;
                    }
                    else{
                        v_next_3D[i][j][k] = v_first_step_3D[i][j][k] - dt / hy * ( p_next_3D[i][j][k] - p_next_3D[i][j - 1][k] );
                                                                
                    }
 //                   if (i==50&&j==1&&k==50)
  //                      printf("%e %e \n", v_first_step_3D[i][j][k],dt / hy * ( p_next_3D[i][j][k] - p_next_3D[i][j - 1][k]));


                    if ( i == 0 ){
                        u_next_3D[i][j][k] = 0.0;
                    }
                    else{
                        u_next_3D[i][j][k] = u_first_step_3D[i][j][k] - dt / hx * ( p_next_3D[i][j][k] - p_next_3D[i - 1][j][k] );
                    }



                    if ( k == 0 ){
                        w_next_3D[i][j][k] = 0.0;
                    }
                    else{
                        w_next_3D[i][j][k] = w_first_step_3D[i][j][k] - dt / hz * ( p_next_3D[i][j][k] - p_next_3D[i][j][k - 1] );
                    }
                    u_next_3D[params2d.cells_number_x][j][k] = 0.0;
                    v_next_3D[i][params2d.cells_number_y][k] = 0.0;
                    w_next_3D[i][j][params2d.cells_number_z] = 0.0;

                }
	    }
	}
    }
    else{                        
        for ( int j = 0; j <  params2d.cells_number_y; j++ ){
	        for ( int i = 0; i < params2d.cells_number_x; i++ ){
                for ( int k = 0; k < params2d.cells_number_z; k++ ){
                    if ( j == 0){
                        v_next_3D[i][j][k] = v_first_step_3D[i][j][k];
                        u_next_3D[i][j][k] = u_first_step_3D[i][j][k];
                        w_next_3D[i][j][k] = w_first_step_3D[i][j][k];
                    }
                    else{
                        if ( j == 1)
                            v_next_3D[i][j][k] = v_first_step_3D[i][params2d.cells_number_y - 1][k] - dt / hy * ( p_next_3D[i][1][k] - p_next_3D[i][params2d.cells_number_y - 2][k] );
                        else
                        if ( j == params2d.cells_number_y - 1 )
                            v_next_3D[i][j][k] = v_first_step_3D[i][j][k] - dt / hy * ( p_next_3D[i][1][k] - p_next_3D[i][j - 1][k] );
                        else

                            v_next_3D[i][j][k] = v_first_step_3D[i][j][k] - dt / hy * ( p_next_3D[i][j][k] - p_next_3D[i][j - 1][k] );
                                                                
                    }



                    if ( i == 0 ){
                        v_next_3D[i][j][k] = v_first_step_3D[i][j][k];
                        u_next_3D[i][j][k] = u_first_step_3D[i][j][k];
                        w_next_3D[i][j][k] = w_first_step_3D[i][j][k];
                    }
                    else{
                        if ( i == 1)
                            u_next_3D[i][j][k] = u_first_step_3D[params2d.cells_number_x - 1][j][k] - dt / hx * ( p_next_3D[1][j][k] - p_next_3D[params2d.cells_number_x - 2][j][k] );
                        else
                        if ( i == params2d.cells_number_x - 1 )
                            u_next_3D[i][j][k] = u_first_step_3D[i][j][k] - dt / hx * ( p_next_3D[1][j][k] - p_next_3D[i - 1][j][k] );
                        else
                            u_next_3D[i][j][k] = u_first_step_3D[i][j][k] - dt / hx * ( p_next_3D[i][j][k] - p_next_3D[i - 1][j][k] );
                    }



                    if ( k == 0 ){
                        v_next_3D[i][j][k] = v_first_step_3D[i][j][k];
                        u_next_3D[i][j][k] = u_first_step_3D[i][j][k];
                        w_next_3D[i][j][k] = w_first_step_3D[i][j][k];
                    }
                    else{

                        if ( k == 1)
                            w_next_3D[i][j][k] = w_first_step_3D[i][j][params2d.cells_number_z - 1] - dt / hx * ( p_next_3D[i][j][1] - p_next_3D[i][j][params2d.cells_number_z - 2] );
                        else
                        if ( k == params2d.cells_number_z - 1 )
                            w_next_3D[i][j][k] = w_first_step_3D[i][j][k] - dt / hx * ( p_next_3D[i][j][1] - p_next_3D[i][j][k - 1] );
                        else
                            w_next_3D[i][j][k] = w_first_step_3D[i][j][k] - dt / hz * ( p_next_3D[i][j][k] - p_next_3D[i][j][k - 1] );
                    }
                    u_next_3D[params2d.cells_number_x][j][k] = u_first_step_3D[params2d.cells_number_x][j][k];
                    v_next_3D[params2d.cells_number_x - 1][j][k] = v_first_step_3D[params2d.cells_number_x - 1][j][k];
                    w_next_3D[params2d.cells_number_x - 1][j][k] = w_first_step_3D[params2d.cells_number_x - 1][j][k];

                    v_next_3D[i][params2d.cells_number_y][k] = v_first_step_3D[i][params2d.cells_number_y][k];
                    u_next_3D[i][params2d.cells_number_y - 1][k] = u_first_step_3D[i][params2d.cells_number_y - 1][k];
                    w_next_3D[i][params2d.cells_number_y - 1][k] = w_first_step_3D[i][params2d.cells_number_y - 1][k];

                    w_next_3D[i][j][params2d.cells_number_z] = w_first_step_3D[i][j][params2d.cells_number_z];
                    u_next_3D[i][j][params2d.cells_number_z - 1] = u_first_step_3D[i][j][params2d.cells_number_z - 1];
                    v_next_3D[i][j][params2d.cells_number_z - 1] = v_first_step_3D[i][j][params2d.cells_number_z - 1];



                }
	    }
        }   
    }                       
}