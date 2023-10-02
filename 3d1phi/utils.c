// utils.cc
// Рабочие функции, специфичные для рассматриваемого класса систем уравнений типа Baer-Nunziato и Saurel-Abgrall
// вне зависимости от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 26 августа 2018 г.

#include "utils.h"






// Приведение параметров задачи к безразмерному виду
// paramsc - структура с основными параметрами вычислительного эксперимента, часть полей проходит процедуру приведения к безразмерному виду (in/out)
// params1d - структура с параметрами одномерной задачи, часть полей проходит процедуру приведения к безразмерному виду (in/out)
//void dimensionalization( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d ) {
//	
//    double density_scale = paramsc->mass_scale / pow( paramsc->length_scale, 3.0 ); // характерный масштаб плотности
//    double velocity_scale = paramsc->length_scale / paramsc->time_scale; // характерный масштаб скорости
//    double pressure_scale = paramsc->mass_scale / paramsc->length_scale / pow( paramsc->time_scale, 2.0 ); // характерный масштаб давления
//
//    double energy_scale = pressure_scale / density_scale;
//
//    // шаг по времени
//    if ( paramsc->constant_time_step )
//        paramsc->dt /= paramsc->time_scale;
//
//    // УРС
//    paramsc->p01 /= pressure_scale; // константа в уравнении состояния дисперсной фазы
//    paramsc->p02 /= pressure_scale; // константа в уравнении состояния газовой фазы
//    // коэффициент динамической вязкости в газе
//    paramsc->mu2 /= ( paramsc->mass_scale / ( paramsc->length_scale * paramsc->time_scale ) );
//
//    // вывод результатов
//    paramsc->stop_time /= paramsc->time_scale;
//
//    // фоновые параметры для особого случая отсутствия дисперсной фазы по одну из сторон от разрыва
//    paramsc->background_density /= density_scale; // плотность дисперсной фазы
//    paramsc->background_velocity /= velocity_scale; // скорость дисперсной фазы
//    paramsc->background_pressure /= pressure_scale; // давление дисперсной фазы
//
//    // "физика"
//    paramsc->particle_diameter /= paramsc->length_scale; // диаметр частицы дисперсной фазы
//    paramsc->compaction_viscosity /= (pressure_scale * paramsc->time_scale);
//    paramsc->interface_drag_coef /= (density_scale / paramsc->time_scale);
//	
//    paramsc->tau_parameter /= (energy_scale / paramsc->mass_scale);
//
//    if (paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC){
//
//        // границы расчетной области
//        params1d->left_boundary_x /= paramsc->length_scale;
//        params1d->right_boundary_x /= paramsc->length_scale;
//
//        // начальные условия
//        for ( int i_block = 0; i_block < params1d->ic_blocks_number; i_block++ ) {
//            params1d->block_values[i_block][R_DISP] /= density_scale; // плотность дисперсной фазы
//            params1d->block_values[i_block][V_DISP] /= velocity_scale; // скорость дисперсной фазы
//            params1d->block_values[i_block][P_DISP] /= pressure_scale; // давление дисперсной фазы
//            params1d->block_values[i_block][R_GAS] /= density_scale; // плотность газовой фазы
//            params1d->block_values[i_block][V_GAS] /= velocity_scale; // скорость газовой фазы
//            params1d->block_values[i_block][P_GAS] /= pressure_scale; // давление газовой фазы
//        }
//    
//        // параметры вдува
//        if ( params1d->left_bc == INFLOW || params1d->right_bc == INFLOW || params1d->left_bc == COMPLEX_INFLOW || params1d->right_bc == COMPLEX_INFLOW ) {
//            params1d->inflow_values[R_DISP] /= density_scale; // плотность дисперсной фазы
//            params1d->inflow_values[V_DISP] /= velocity_scale; // скорость дисперсной фазы
//            params1d->inflow_values[P_DISP] /= pressure_scale; // давление дисперсной фазы
//            params1d->inflow_values[R_GAS] /= density_scale; // плотность газовой фазы
//            params1d->inflow_values[V_GAS] /= velocity_scale; // скорость газовой фазы
//            params1d->inflow_values[P_GAS] /= pressure_scale; // давление газовой фазы
//        }
//
//        // параметры сложного вдува
//        if ( params1d->left_bc == COMPLEX_INFLOW || params1d->right_bc == COMPLEX_INFLOW ) {
//            params1d->change_inflow_time /= paramsc->time_scale;
//            params1d->changed_inflow_values[R_DISP] /= density_scale; // плотность дисперсной фазы
//            params1d->changed_inflow_values[V_DISP] /= velocity_scale; // скорость дисперсной фазы
//            params1d->changed_inflow_values[P_DISP] /= pressure_scale; // давление дисперсной фазы
//            params1d->changed_inflow_values[R_GAS] /= density_scale; // плотность газовой фазы
//            params1d->changed_inflow_values[V_GAS] /= velocity_scale; // скорость газовой фазы
//            params1d->changed_inflow_values[P_GAS] /= pressure_scale; // давление газовой фазы
//        }
//    }
//    else if (paramsc->program_name == TWOD2PHC){
//
//        double velocity_scale1 = paramsc->length_scale / paramsc->time_scale;    /* характерный масштаб скорости */
//        double velocity_scale2 = paramsc->length_scale / paramsc->time_scale;   // (S) Ввёл 2-ую функцию для скорости просто для различия (которое возможно и не нужно)
//
//        /* границы расчетной области */
//        params2d->left_boundary_x /= paramsc->length_scale;
//        params2d->right_boundary_x /= paramsc->length_scale;
//
//        /* начальные условия */
//        for ( int i_block = 0; i_block < params2d->x_blocks_number; i_block++ ) 
//        {
//            for ( int j_block = 0; j_block < params2d->y_blocks_number; j_block++ ) 
//            {
//                params2d->block_values[i_block][j_block][R_DISP_2D] /= density_scale;     /* плотность дисперсной фазы */
//                params2d->block_values[i_block][j_block][V_DISP_2D] /= velocity_scale1;    /* скорость дисперсной фазы */
//                params2d->block_values[i_block][j_block][U_DISP_2D] /= velocity_scale2;    /* скорость дисперсной фазы */
//                params2d->block_values[i_block][j_block][P_DISP_2D] /= pressure_scale;    /* давление дисперсной фазы */
//                params2d->block_values[i_block][j_block][R_GAS_2D] /= density_scale;     /* плотность газовой фазы */
//                params2d->block_values[i_block][j_block][V_GAS_2D] /= velocity_scale1;    /* скорость газовой фазы */
//                params2d->block_values[i_block][j_block][U_GAS_2D] /= velocity_scale2;    /* скорость газовой фазы */
//                params2d->block_values[i_block][j_block][P_GAS_2D] /= pressure_scale;    /* давление газовой фазы */
//            }
//        }
//
//        if ( params2d->left_bc == INFLOW ) 
//        {
//            params2d->inflow_values_left_bc[R_DISP_2D] /= density_scale;     
//            params2d->inflow_values_left_bc[V_DISP_2D] /= velocity_scale1;   
//            params2d->inflow_values_left_bc[U_DISP_2D] /= velocity_scale2;    
//            params2d->inflow_values_left_bc[P_DISP_2D] /= pressure_scale;    
//            params2d->inflow_values_left_bc[R_GAS_2D] /= density_scale;     
//            params2d->inflow_values_left_bc[V_GAS_2D] /= velocity_scale1;   
//            params2d->inflow_values_left_bc[U_GAS_2D] /= velocity_scale2;    
//            params2d->inflow_values_left_bc[P_GAS_2D] /= pressure_scale;    
//        }
//        if ( params2d->right_bc == INFLOW ) 
//        {
//            params2d->inflow_values_right_bc[R_DISP_2D] /= density_scale;     
//            params2d->inflow_values_right_bc[V_DISP_2D] /= velocity_scale1;   
//            params2d->inflow_values_right_bc[U_DISP_2D] /= velocity_scale2;    
//            params2d->inflow_values_right_bc[P_DISP_2D] /= pressure_scale;    
//            params2d->inflow_values_right_bc[R_GAS_2D] /= density_scale;     
//            params2d->inflow_values_right_bc[V_GAS_2D] /= velocity_scale1;   
//            params2d->inflow_values_right_bc[U_GAS_2D] /= velocity_scale2;    
//            params2d->inflow_values_right_bc[P_GAS_2D] /= pressure_scale;    
//        }
//        if ( params2d->down_bc == INFLOW ) 
//        {
//            params2d->inflow_values_down_bc[R_DISP_2D] /= density_scale;     
//            params2d->inflow_values_down_bc[V_DISP_2D] /= velocity_scale1;   
//            params2d->inflow_values_down_bc[U_DISP_2D] /= velocity_scale2;    
//            params2d->inflow_values_down_bc[P_DISP_2D] /= pressure_scale;    
//            params2d->inflow_values_down_bc[R_GAS_2D] /= density_scale;     
//            params2d->inflow_values_down_bc[V_GAS_2D] /= velocity_scale1;   
//            params2d->inflow_values_down_bc[U_GAS_2D] /= velocity_scale2;    
//            params2d->inflow_values_down_bc[P_GAS_2D] /= pressure_scale;    
//        }
//        if ( params2d->up_bc == INFLOW ) 
//        {
//            params2d->inflow_values_up_bc[R_DISP_2D] /= density_scale;     
//            params2d->inflow_values_up_bc[V_DISP_2D] /= velocity_scale1;   
//            params2d->inflow_values_up_bc[U_DISP_2D] /= velocity_scale2;    
//            params2d->inflow_values_up_bc[P_DISP_2D] /= pressure_scale;    
//            params2d->inflow_values_up_bc[R_GAS_2D] /= density_scale;     
//            params2d->inflow_values_up_bc[V_GAS_2D] /= velocity_scale1;   
//            params2d->inflow_values_up_bc[U_GAS_2D] /= velocity_scale2;    
//            params2d->inflow_values_up_bc[P_GAS_2D] /= pressure_scale;    
//        }
//        if ( params2d->left_bc == COMPLEX_INFLOW ) 
//        {
//            params2d->inflow_values_left_bc[R_DISP_2D] /= density_scale;     
//            params2d->inflow_values_left_bc[V_DISP_2D] /= velocity_scale1;   
//            params2d->inflow_values_left_bc[U_DISP_2D] /= velocity_scale2;    
//            params2d->inflow_values_left_bc[P_DISP_2D] /= pressure_scale;    
//            params2d->inflow_values_left_bc[R_GAS_2D] /= density_scale;     
//            params2d->inflow_values_left_bc[V_GAS_2D] /= velocity_scale1;   
//            params2d->inflow_values_left_bc[U_GAS_2D] /= velocity_scale2;    
//            params2d->inflow_values_left_bc[P_GAS_2D] /= pressure_scale;    
//        }
//        if ( params2d->right_bc == COMPLEX_INFLOW ) 
//        {
//            params2d->inflow_values_right_bc[R_DISP_2D] /= density_scale;     
//            params2d->inflow_values_right_bc[V_DISP_2D] /= velocity_scale1;   
//            params2d->inflow_values_right_bc[U_DISP_2D] /= velocity_scale2;    
//            params2d->inflow_values_right_bc[P_DISP_2D] /= pressure_scale;    
//            params2d->inflow_values_right_bc[R_GAS_2D] /= density_scale;     
//            params2d->inflow_values_right_bc[V_GAS_2D] /= velocity_scale1;   
//            params2d->inflow_values_right_bc[U_GAS_2D] /= velocity_scale2;    
//            params2d->inflow_values_right_bc[P_GAS_2D] /= pressure_scale;    
//        }
//
//    }
//}

//void current_block_number(struct ParametersCommon *params, struct Parameters1d *params1d, struct Parameters2d *params2d, int step_number_x, int step_number_y, int step_number_z, int *number_of_block){
//
//    number_of_block[0] = 0;
//
//    if (params->program_name == ONED2PHC || params->program_name == ONED3PHC){
//        for (int i = 0; i < params1d->ic_blocks_number; i++){
//	    if (step_number_x <= params1d->cell_end[i] && step_number_x >= params1d->cell_begin[i])
//	        number_of_block[0] = i;
//        }
//    }
//    else if (params->program_name == TWOD2PHC){
//        number_of_block[1] = 0;
//        for (int i = 0; i < params2d->x_blocks_number; i++){
//            for (int j = 0; j < params2d->y_blocks_number; j++){
//                for (int k = 0; k < params2d->z_blocks_number; k++){
//                if (step_number_x <= params2d->x_end[i][j][k] && step_number_x >= params2d->x_begin[i][j][k] && step_number_y <= params2d->y_end[i][j][k] && step_number_z >= params2d->z_begin[i][j][k] && step_number_z >= params2d->z_begin[i][j][k]){
//                    number_of_block[0] = i;
//                    number_of_block[1] = j;
//                    number_of_block[2] = k;
//                }
//            }
//        }
//    }
//    else{
//        printf("Wrong program name is used in function current_block_number_2d\n");
//        exit(EXIT_FAILURE);
//    }
//   
//}

								//if ( cell_type[i + 1][j] == OUTER ){								
								//		u_first_step_2D[i + 1][j] = 0.0;
								//		if ( cell_type[i][j + 1] == OUTER ){
								//				v_first_step_2D[i][j + 1] = 0.0;
								//		}

								//		else
								//			v_first_step_2D[i][j + 1] = v_prev_2D[i][j + 1] / 2.0 + v_prev_2D[i - 1][j + 1] / 6.0;																																																
								//}
								//else
								//if ( cell_type[i - 1][j] == OUTER ){

								//		if ( cell_type[i][j + 1] == OUTER ){
								//				v_first_step_2D[i][j + 1] = 0.0;
								//				u_first_step_2D[i + 1][j] = u_prev_2D[i + 1][j] / 2.0 + u_prev_2D[i + 1][j - 1] / 6.0;
								//		}
								//		else
								//				if ( cell_type[i][j - 1] == OUTER ){
								//						v_first_step_2D[i][j + 1] = v_prev_2D[i][j + 1] / 2.0 + v_prev_2D[i + 1][j + 1] / 6.0;
								//						u_first_step_2D[i + 1][j] = u_prev_2D[i + 1][j] / 2.0 + u_prev_2D[i + 1][j + 1] / 6.0;
								//				}
								//				else{
								//						v_first_step_2D[i][j + 1] = v_prev_2D[i][j + 1] / 2.0 + v_prev_2D[i + 1][j + 1] / 6.0;
								//						u_first_step_2D[i + 1][j] = u_prev_2D[i + 1][j] + dt * ( ( u_prev_2D[i][j] * u_prev_2D[i + 1][j]  - u_prev_2D[i + 1][j] * u_prev_2D[i + 2][j] ) / hx + 
								//								( ( u_prev_2D[i + 1][j] + u_prev_2D[i + 1][j - 1] ) / 2.0 * ( v_prev_2D[i + 1][j] + v_prev_2D[i][j] ) / 2.0 - ( u_prev_2D[i + 1][j] + u_prev_2D[i + 1][j + 1] ) / 2.0 * ( v_prev_2D[i + 1][j + 1] + v_prev_2D[i][j + 1] ) / 2.0 ) / hy -
								//								paramsc.mu / hy * ( ( ( v_prev_2D[i + 1][j + 1] - v_prev_2D[i][j + 1] ) / hx - ( u_prev_2D[i + 1][j + 1] - u_prev_2D[i + 1][j] ) / hy ) -
								//								( ( v_prev_2D[i + 1][j] - v_prev_2D[i][j] ) / hx - ( u_prev_2D[i + 1][j] - u_prev_2D[i + 1][j - 1] ) / hy ) ) );														
								//				}
								//}
								//else
								//if ( cell_type[i][j + 1] == OUTER ){								
								//		v_first_step_2D[i][j + 1] = 0.0;
								//		if ( cell_type[i + 1][j] == OUTER ){
								//				u_first_step_2D[i + 1][j] = 0.0;
								//		}

								//		else
								//				u_first_step_2D[i + 1][j] = u_prev_2D[i + 1][j] / 2.0 + u_prev_2D[i + 1][j - 1] / 6.0;										
								//}
								//else
								//if ( cell_type[i][j - 1] == OUTER ){								
								//		if ( cell_type[i + 1][j] == OUTER ){
								//				u_first_step_2D[i + 1][j] = 0.0;
								//				v_first_step_2D[i][j + 1] = v_prev_2D[i][j + 1] / 2.0 + v_prev_2D[i - 1][j + 1] / 6.0;
								//		}

								//		else
								//				if ( cell_type[i - 1][j] == OUTER ){
								//						u_first_step_2D[i + 1][j] = u_prev_2D[i + 1][j] / 2.0 + u_prev_2D[i + 1][j + 1] / 6.0;
								//						v_first_step_2D[i][j + 1] = v_prev_2D[i][j + 1] / 2.0 + v_prev_2D[i + 1][j + 1] / 6.0;
								//				}
								//				else{
								//						u_first_step_2D[i + 1][j] = u_prev_2D[i + 1][j] / 2.0 + u_prev_2D[i + 1][j + 1] / 6.0;
								//						v_first_step_2D[i][j + 1] = v_prev_2D[i][j + 1] + dt * ( ( v_prev_2D[i][j + 1] * v_prev_2D[i][j]  - v_prev_2D[i][j + 1] * v_prev_2D[i][j + 2] ) / hy + 
								//								( ( v_prev_2D[i][j + 1] + v_prev_2D[i - 1][j + 1] ) / 2.0 * ( u_prev_2D[i][j] + u_prev_2D[i][j + 1] ) / 2.0 - ( v_prev_2D[i + 1][j + 1] + v_prev_2D[i][j + 1] ) / 2.0 * ( u_prev_2D[i + 1][j] + u_prev_2D[i + 1][j + 1] ) / 2.0 ) / hx +
								//								paramsc.mu / hx * ( ( ( v_prev_2D[i + 1][j + 1] - v_prev_2D[i][j + 1] ) / hx - ( u_prev_2D[i + 1][j + 1] - u_prev_2D[i + 1][j] ) / hy ) -
								//								( ( v_prev_2D[i][j + 1] - v_prev_2D[i - 1][j] ) / hx - ( u_prev_2D[i][j + 1] - u_prev_2D[i][j] ) / hy ) ) );														
								//				}										
								//}
								//else{