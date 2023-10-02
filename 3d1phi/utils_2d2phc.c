// utils_2d.cc
// ������� �������, ����������� ��� ��������� ������� ������� ��������� ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// (c) ����� �����, 2018
// ������: 27 ������� 2018 �.

#include "utils_2d2phc.h"
#include "io_2d2phc.h"



// ������������� �������-�������
// params2d - ��������� � ����������� ��������� ������ (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out)
void init_solution_3d( struct Parameters2d *params2d, struct ParametersCommon* paramsc, int file_num, struct TimeMoment *time_mom, double ***u, double ***v, double ***w ) {
    if (!paramsc->isContinue){
		//printf(" %d %d \n", params2d->x_begin[0][0], params2d->x_end[0][0]);
	printf("\ninit new solution\n");
     for ( int i_block = 0; i_block < params2d->x_blocks_number; i_block++ ) { // ���� �� ������ ����� ��� x
        for ( int j_block = 0; j_block < params2d->y_blocks_number; j_block++ ) { // ���� �� ������ ����� ��� y
            for ( int k_block = 0; k_block < params2d->z_blocks_number; k_block++ ) { // ���� �� ������ ����� ��� z
                for ( int i_cell = params2d->x_begin[i_block][j_block][k_block]; i_cell <= params2d->x_end[i_block][j_block][k_block]; i_cell++ ) { // ���� �� ������� ������ ����� �� x  
                    for ( int j_cell = params2d->y_begin[i_block][j_block][k_block]; j_cell <= params2d->y_end[i_block][j_block][k_block]; j_cell++ ) { // ���� �� ������� ������ ����� �� y 
                        for ( int k_cell = params2d->z_begin[i_block][j_block][k_block]; k_cell <= params2d->z_end[i_block][j_block][k_block]; k_cell++ ) { // ���� �� ������� ������ ����� �� y
			                u[i_cell][j_cell][k_cell] = params2d->block_values_u[i_block][j_block][k_block];
			                v[i_cell][j_cell][k_cell] = params2d->block_values_v[i_block][j_block][k_block];
                            w[i_cell][j_cell][k_cell] = params2d->block_values_w[i_block][j_block][k_block];
			                u[i_cell + 1][j_cell][k_cell] = params2d->block_values_u[i_block][j_block][k_block];
			                v[i_cell][j_cell + 1][k_cell] = params2d->block_values_v[i_block][j_block][k_block];
                            w[i_cell][j_cell][k_cell + 1] = params2d->block_values_w[i_block][j_block][k_block];
                            //if ( j_cell < 20 && j_cell > 10 && k_cell < 20 && k_cell > 10 )
                            //u[38][j_cell][k_cell] = 1.0;
                            //printf(" %lf\n", u[i_cell][j_cell][k_cell + 1]);
                        } // k_cell
                    } // j_cell
                } // i_cell
            } // k_block
        } // j_block
     } // i_block
    }else{
    printf("\nread solution\n");
    read_solution(params2d, paramsc, file_num, u, v, w);


    }
    

}



// ������ ���� �������������� �� �������
// params - ��������� � ����������� ��������������� ������������
// *x - ���������� ����� �����
// **v_ncons - ������� ����������� ���������� � ������� �����
// ���������� ��� �������������� �� �������
double get_time_step_3d( const struct ParametersCommon* params, const struct Parameters2d* params2d, const double* x, const double* y, const double* z, double ***u, double ***v, double ***w ) 
{  
    int i, j, k;
    double h_x, h_y, h_z;
    double new_step = INFINITY1; // �������������� ��� �� �������

    double curr_step;



    if ( params->constant_time_step ) 
    {
        // ������ � ���������� ����� �� �������
        return params->dt;
    }
    else 
    {
        // ������ � ������������ ������� ���� �� �������
        for ( i = 0; i < params2d->cells_number_x; i++ )   
        {      
           for ( j = 0; j < params2d->cells_number_y; j++ )   
           {
               for (k = 0; k < params2d->cells_number_z; k++ )
               {


                   h_x = x[i+1] - x[i]; // ������� ������ ������ �� x
                   h_y = y[j+1] - y[j]; // ������� ������ ������ �� y  //(S) ��� ������ �� j ������
                   h_z = z[k+1] - z[k]; // ������� ������ ������ �� x


                   // ������������� ���� �� ������� 
                   curr_step = 6 * params->mu / ( u[i][j][k] * u[i][j][k] + v[i][j][k] * v[i][j][k] + w[i][j][k] * w[i][j][k] );
                   //printf("step %lf %lf %lf\n", u[i][j][k] * u[i][j][k], v[i][j][k] * v[i][j][k], w[i][j][k] * w[i][j][k]  );
                   if ( curr_step < new_step )
                       new_step = curr_step;
                   curr_step = ( 0.5 / ( 1 / h_x / h_x + 1 / h_y / h_y + 1 / h_z / h_z ) ) / params->mu;
			  //printf("step %lf\n",curr_step  );
                   if ( curr_step < new_step )
                       new_step = curr_step;
               }
            };
        //return params->cfl * new_step;
       };
		//printf("step %lf\n",params->cfl * new_step);
       return params->cfl * new_step;
	   
    };
}
