/******************************************************************************
 * Copyright (c) 1998 Lawrence Livermore National Security, LLC and other
 * HYPRE Project Developers. See the top-level COPYRIGHT file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 ******************************************************************************/



#include "mem.h"
#include "math_utils.h"
#include "io_2d2phc.h"
#include "utils.h"
#include "utils_2d2phc.h"
#include "io.h"
#include "first_step.h"
#include "first_step_rhs.h"
#include "second_step.h"
#include "second_step_rhs.h"
#include "third_step.h"
#include "poisson.h"
#include "ex5.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "ex.h"
#include "struct.h"
#include <time.h>

//#include "ex5.h"

#ifdef HYPRE_EXVIS
#include "vis.c"
#endif

int hypre_FlexGMRESModifyPCAMGExample(void *precond_data, int iterations,
                                      double rel_residual_norm);

#define my_min(a,b)  (((a)<(b)) ? (a) : (b))

int main ( int argc, char *argv[] )
{
    int myid, num_procs;

    printf( "\n3D one-phase incompressible solver \n" );
    printf( "\nPreparing:\n" );

    // ��������� ���������� ��������� ������
    // if ( argc != 3 ) {
    // argv[1] - ���� � ����������, ��� ��������� ���� � ����������� ������
    // argv[2] - ���� � ����������, ��� ����� �������� �������� �����
    //     printf( "\nmain -> wrong command line arguments number\n" );
    //     exit( EXIT_FAILURE );
    // }


    printf( "\n Command line arguments read success %d %d\n", myid, num_procs);
    struct ParametersCommon paramsc; //
    struct Parameters2d params2d; //
    double dt; //
    double time_gap; //
    int steps_gap; //
    int files_counter = 0; //
    struct TimeMoment time_mom = { 0, 0.0 }; //
    struct DebugInfo debug_info; // �
    char output_filename[MAX_STRING_SIZE] = " "; //
    char tmp_str[MAX_STRING_SIZE]; //
    paramsc.isContinue = false;
    int file_num = 0;
    if (atoi(argv[3]) == 1){
    printf("restart\n");
	paramsc.isContinue = true;
	file_num = atoi(argv[4]);
    }else{
        paramsc.isContinue = false;
        file_num = 0.0;
    }
    printf("\n%s %d\n",argv[3], file_num);
    fill_parameters_struct_2d( &paramsc, &params2d, argv[1], argv[2], 0 );
    printf( "\n> File parameters.dat is processed successfully. Correspondent structure is filled. \n" );


    double ***u_prev_3D; // ������� �-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x+1, params2d.cells_number_y, params2d.cells_number_z, &u_prev_3D );
    double ***v_prev_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y+1, params2d.cells_number_z, &v_prev_3D );
    double ***w_prev_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z+1, &w_prev_3D );
    double ***u_inter_3D; // ������� �-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x+1, params2d.cells_number_y, params2d.cells_number_z, &u_inter_3D );
    double ***v_inter_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y+1, params2d.cells_number_z, &v_inter_3D );
    double ***w_inter_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z+1, &w_inter_3D );
    double *p_prev_3D; // �������� � ��������� ������� �� n-�� ����
    if (BOUN_TYPE == 1)
        get_memory_for_1D_double_array( (params2d.cells_number_x-2) * (params2d.cells_number_y-2) * (params2d.cells_number_z-2), &p_prev_3D );
    if (BOUN_TYPE == 0)
        get_memory_for_1D_double_array( (params2d.cells_number_x) * (params2d.cells_number_y) * (params2d.cells_number_z), &p_prev_3D );
    double ***u_first_step_3D; // ������� �-���������� �������� � ��������� ������� ����� ���������� ������� ����� �����������
    get_memory_for_3D_double_array( params2d.cells_number_x+1, params2d.cells_number_y, params2d.cells_number_z, &u_first_step_3D );
    double ***v_first_step_3D; // ������� y-���������� �������� � ��������� ������� ����� ���������� ������� ����� �����������
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y+1, params2d.cells_number_z, &v_first_step_3D );
    double ***w_first_step_3D; // ������� y-���������� �������� � ��������� ������� ����� ���������� ������� ����� �����������
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z+1, &w_first_step_3D );
    double ***p_next_3D; // �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z, &p_next_3D );
    double ***u_next_3D; // ������� �-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x+1, params2d.cells_number_y, params2d.cells_number_z, &u_next_3D );
    double ***v_next_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y+1, params2d.cells_number_z, &v_next_3D );
    double ***w_next_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z+1, &w_next_3D );

    double ***nonlin_x_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z, &nonlin_x_3D );
    double ***nonlin_y_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z, &nonlin_y_3D );
    double ***nonlin_z_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z, &nonlin_z_3D );
    init_nonlin(nonlin_x_3D, nonlin_y_3D, nonlin_z_3D, params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z);

    double ***force_x_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z, &force_x_3D );
    double ***force_y_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z, &force_y_3D );
    double ***force_z_3D; // ������� y-���������� �������� � ��������� ������� �� n-�� ����
    get_memory_for_3D_double_array( params2d.cells_number_x, params2d.cells_number_y, params2d.cells_number_z, &force_z_3D );
    double N_fs = ( params2d.cells_number_x + 1 ) * ( params2d.cells_number_y) * ( params2d.cells_number_z ) +
   			( params2d.cells_number_x ) * ( params2d.cells_number_y + 1) * ( params2d.cells_number_z ) +
			( params2d.cells_number_x ) * ( params2d.cells_number_y ) * ( params2d.cells_number_z + 1 );

    double **matrix; //������� ��� ������� ��������� ��������
    int cell_num = params2d.cells_number_x * params2d.cells_number_y * params2d.cells_number_z;
    get_memory_for_2D_double_array( cell_num, 7, &matrix );
    double **u_mean;
    get_memory_for_2D_double_array( params2d.cells_number_x, params2d.cells_number_y, &u_mean );
    double **v_mean;
    get_memory_for_2D_double_array( params2d.cells_number_x, params2d.cells_number_y, &v_mean );
    double **w_mean;
    get_memory_for_2D_double_array( params2d.cells_number_x, params2d.cells_number_y, &w_mean );

    //double *b; //������ ����� ��� ��������� ��������
    //get_memory_for_1D_double_array( cell_num, &b );
    double *solution; //������ ����� ��� ��������� ��������
    get_memory_for_1D_double_array( cell_num, &solution );
    double *x_grid; // ������ ��������� ����� ����� x
    get_memory_for_1D_double_array( params2d.cells_number_x + 1, &x_grid );
    double *y_grid; // ������ ��������� ����� ����� �� y
    get_memory_for_1D_double_array( params2d.cells_number_y + 1, &y_grid );
    double *z_grid; // ������ ��������� ����� ����� �� y
    get_memory_for_1D_double_array( params2d.cells_number_z + 1, &z_grid );

    printf( "\n> All the necessary memory is allocated successfully.\n" );

    for ( int i = 0; i < params2d.cells_number_x + 1; i++ )
        x_grid[i] = params2d.left_boundary_x + ( params2d.right_boundary_x - params2d.left_boundary_x ) / ( params2d.cells_number_x ) * i;
    for ( int j = 0; j < params2d.cells_number_y + 1; j++ )
        y_grid[j] = params2d.down_boundary_y + ( params2d.up_boundary_y - params2d.down_boundary_y ) / ( params2d.cells_number_y ) * j;
    for ( int k = 0; k < params2d.cells_number_z + 1; k++ )
        z_grid[k] = params2d.rear_boundary_z + ( params2d.front_boundary_z - params2d.rear_boundary_z ) / ( params2d.cells_number_z ) * k;

    double hx = x_grid[1] - x_grid[0];
    double hy = y_grid[1] - y_grid[0];
    double hz = z_grid[1] - z_grid[0];
    init_solution_3d( &params2d, &paramsc, file_num, /* argv[1],*/ &time_mom, u_prev_3D, v_prev_3D, w_prev_3D);
    files_counter = file_num;

    for ( int i = 0; i < params2d.cells_number_x; i++){
        for (int j = 0; j < params2d.cells_number_y; j++){
    	    u_mean[i][j] = 0.0;
    	    v_mean[i][j] = 0.0;
    	    w_mean[i][j] = 0.0;
            for (int k = 0; k < params2d.cells_number_z; k++){
        	p_next_3D[i][j][k] = 0.0;
                u_inter_3D[i][j][k] = u_prev_3D[i][j][k];
                v_inter_3D[i][j][k] = v_prev_3D[i][j][k];
                w_inter_3D[i][j][k] = w_prev_3D[i][j][k];
            }
        }
    }
    for (int j = 0; j < params2d.cells_number_y; j++){
        for (int k = 0; k < params2d.cells_number_z; k++){
            u_inter_3D[params2d.cells_number_x][j][k] = u_prev_3D[params2d.cells_number_x][j][k];
        }
    }
    for (int i = 0; i < params2d.cells_number_x; i++){
        for (int k = 0; k < params2d.cells_number_z; k++){
            v_inter_3D[i][params2d.cells_number_y][k] = v_prev_3D[i][params2d.cells_number_y][k];
        }
    }
    for (int i = 0; i < params2d.cells_number_x; i++){
        for (int j = 0; j < params2d.cells_number_y; j++){
            w_inter_3D[i][j][params2d.cells_number_z] = w_prev_3D[i][j][params2d.cells_number_z];
        }
    }
    printf("\n omega = %e", paramsc. omega);
    printf("\n force = %e", paramsc.f0);
    std::srand(time(0));
    int qvx, qvy, qvz;
    double Bx, By;
    FILE *kinetic;
	char string[MAX_STRING_SIZE];
	strcpy( string, argv[2] );
	strcat(string, "/energy.dat");
	kinetic = fopen( string, "w" );

	if ( NULL == kinetic )
	{
		printf( "\nwrite_solution -> can't open file %s for writing\n", output_filename );
	}
	fprintf(kinetic, "time energy\n");
	fclose(kinetic);

	FILE *vel;
	char string2[MAX_STRING_SIZE];
	strcpy( string2, argv[2] );
	strcat(string2, "/velocity.dat");
	vel = fopen( string2, "w" );

	if ( NULL == vel )
	{
		printf( "\nwrite_solution -> can't open file %s for writing\n", output_filename );
	}
	fprintf(vel, "time U\n");
	fclose(vel);

    for (int l = 0; l < cell_num; l ++ )
        solution[l] = 0.0;
    printf( "\n> The solution is initiated successfully.\n" );
    if (!paramsc.isContinue){
    	write_solution_3d( &paramsc, &params2d, argv[2], x_grid, y_grid, z_grid, p_next_3D, u_prev_3D, v_prev_3D, w_prev_3D, 100, output_filename, M2D );
    	files_counter++;
    }
    printf( "\nCalculation:\n" );
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    dt = get_time_step_3d( &paramsc, &params2d, x_grid, y_grid, z_grid, u_prev_3D, v_prev_3D, w_prev_3D );
    init_force(paramsc, params2d, force_x_3D, force_y_3D, force_z_3D, N_fs, x_grid, y_grid, z_grid, dt);
    HYPRE_IJVector x_fs;
    HYPRE_ParVector par_x_fs;
    HYPRE_IJVector b_fs;
    HYPRE_ParVector par_b_fs;

    HYPRE_IJMatrix A_fs;
    HYPRE_ParCSRMatrix parcsr_A_fs;
    HYPRE_Solver solver_fs, precond_fs;
    HYPRE_Init();
    calc_first_step(paramsc, params2d, p_next_3D, u_inter_3D, v_inter_3D, w_inter_3D, u_prev_3D, v_prev_3D, w_prev_3D, u_first_step_3D, v_first_step_3D,
                        w_first_step_3D, x_grid, y_grid, z_grid, dt, time_mom.curr_t, myid, num_procs, &x_fs, &par_x_fs,
						&b_fs, &par_b_fs, &A_fs, &parcsr_A_fs);


	int solver_id = 0;
;
   	int cx = params2d.cells_number_x;
   	int cy = params2d.cells_number_y;
   	int cz = params2d.cells_number_z;
    int ilower, iupper;
    int local_size, extra;
    local_size = N_fs / num_procs;
    extra = N_fs - local_size * num_procs;

    ilower = local_size * myid;
    ilower += my_min(myid, extra);

    iupper = local_size * (myid + 1);
    iupper += my_min(myid + 1, extra);
    iupper = iupper - 1;

    /* How many rows do I have? */
    local_size = iupper - ilower + 1;
	int num_iterations;
	double final_res_norm;
	/* Create solver */
	HYPRE_BoomerAMGCreate(&solver_fs);
	HYPRE_BoomerAMGSetMaxIter(solver_fs,20);
	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_BoomerAMGSetPrintLevel(solver_fs, 0);  /* print solve info + parameters */
	//HYPRE_BoomerAMGSetOldDefault(solver); /* Falgout coarsening with modified classical interpolaiton */
	HYPRE_BoomerAMGSetRelaxType(solver_fs, 1);   /* G-S/Jacobi hybrid relaxation */
	HYPRE_BoomerAMGSetRelaxOrder(solver_fs, 2);   /* uses C/F relaxation */
	HYPRE_BoomerAMGSetNumSweeps(solver_fs, 1);   /* Sweeeps on each level */
	HYPRE_BoomerAMGSetMaxLevels(solver_fs, 1);  /* maximum number of levels */
	HYPRE_BoomerAMGSetTol(solver_fs, 1e-2);      /* conv. tolerance */

    HYPRE_IJVector x_ss;
    HYPRE_ParVector par_x_ss;
    HYPRE_IJVector b_ss;
    HYPRE_ParVector par_b_ss;

    HYPRE_IJMatrix A_ss;
    HYPRE_ParCSRMatrix parcsr_A_ss;
    HYPRE_Solver solver_ss, precond_ss;
    calc_second_step(paramsc, params2d, u_first_step_3D, v_first_step_3D, w_first_step_3D, p_prev_3D, x_grid, y_grid, z_grid, dt, myid, num_procs, &A_ss, &parcsr_A_ss );

    HYPRE_BoomerAMGCreate(&solver_ss);
    HYPRE_BoomerAMGSetMaxIter(solver_ss,20);

    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_BoomerAMGSetPrintLevel(solver_ss, 0);  /* print solve info + parameters */
    //HYPRE_BoomerAMGSetOldDefault(solver_ss); /* Falgout coarsening with modified classical interpolaiton */
    HYPRE_BoomerAMGSetRelaxType(solver_ss, 1);   /* G-S/Jacobi hybrid relaxation */
    HYPRE_BoomerAMGSetRelaxOrder(solver_ss, 3);   /* uses C/F relaxation */
    HYPRE_BoomerAMGSetNumSweeps(solver_ss, 1);   /* Sweeeps on each level */
    HYPRE_BoomerAMGSetMaxLevels(solver_ss, 10);  /* maximum number of levels */
    HYPRE_BoomerAMGSetTol(solver_ss, 1e-3);      /* conv. tolerance */

    /* Now setup and solve! */
    HYPRE_BoomerAMGSetup(solver_ss, parcsr_A_ss, par_b_ss, par_x_ss);


    if ( paramsc.is_output_on_time ){
        time_gap = paramsc.stop_time / paramsc.output_number;
        time_mom.curr_t = file_num *  time_gap;
    }
    else{
        steps_gap = paramsc.stop_steps_number / paramsc.output_number;
        time_mom.curr_t = file_num *  steps_gap;
    }
    while ( ( ( paramsc.stop_steps_number - time_mom.steps_num > 0 ) && ( !paramsc.is_output_on_time ) ) ||
            ( ( paramsc.stop_time - time_mom.curr_t > 0 ) && ( paramsc.is_output_on_time ) ) )
    {
        dt = get_time_step_3d( &paramsc, &params2d, x_grid, y_grid, z_grid, u_prev_3D, v_prev_3D, w_prev_3D );
        printf("\nbegin first step\n");

        calc_first_step_rhs(paramsc, params2d, p_next_3D, u_inter_3D, v_inter_3D, w_inter_3D, u_prev_3D, v_prev_3D, w_prev_3D, u_first_step_3D, v_first_step_3D,
                            w_first_step_3D, x_grid, y_grid, z_grid, dt, time_mom.curr_t, myid, num_procs, &x_fs, &par_x_fs,
							&b_fs, &par_b_fs, nonlin_x_3D, nonlin_y_3D, nonlin_z_3D, force_x_3D, force_y_3D, force_z_3D);


        printf("\n parameters set\n");
        /* Now setup and solve! */
        HYPRE_BoomerAMGSetup(solver_fs, parcsr_A_fs, par_b_fs, par_x_fs);
        printf("\n solver setup\n");
        HYPRE_BoomerAMGSolve(solver_fs, parcsr_A_fs, par_b_fs, par_x_fs);
        printf("\n solver solved\n");
        /* Run info - needed logging turned on */
        HYPRE_BoomerAMGGetNumIterations(solver_fs, &num_iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver_fs, &final_res_norm);
        if (myid == 0)
        {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }

        /* Destroy solver */
        fin_first_step(paramsc, params2d, p_next_3D, u_first_step_3D, v_first_step_3D,
                w_first_step_3D, myid, num_procs, &x_fs);
        HYPRE_IJVectorDestroy(b_fs);
        HYPRE_IJVectorDestroy(x_fs);
        printf("\nfirst step -- done\n");
        calc_second_step_rhs(paramsc, params2d, u_first_step_3D, v_first_step_3D, w_first_step_3D, p_prev_3D, x_grid, y_grid, z_grid, dt, myid, num_procs, &x_ss, &par_x_ss,
				&b_ss, &par_b_ss );


        HYPRE_BoomerAMGSolve(solver_ss, parcsr_A_ss, par_b_ss, par_x_ss);
        HYPRE_BoomerAMGGetNumIterations(solver_ss, &num_iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver_ss, &final_res_norm);
        fin_second_step(paramsc, params2d, p_prev_3D, myid, num_procs, &x_ss);
        HYPRE_IJVectorDestroy(b_ss);
        HYPRE_IJVectorDestroy(x_ss);
        for ( int i = 0; i < params2d.cells_number_x; i++ ){
            for ( int j = 0; j < params2d.cells_number_y; j++ ){
                for ( int k = 0; k < params2d.cells_number_z; k++ ){
                    //p_next_3D[i][j][k] = solution[k * params2d.cells_number_x * params2d.cells_number_y + j * params2d.cells_number_x + i];
                    if (BOUN_TYPE == 1) {
                        if (i != 0 && j != 0 && k != 0 && i != (params2d.cells_number_x - 1) &&
                            j != (params2d.cells_number_x - 1) && k != (params2d.cells_number_x - 1))
                            p_next_3D[i][j][k] = p_prev_3D[
                                    (k - 1) * (params2d.cells_number_x - 2) * (params2d.cells_number_y - 2) +
                                    (j - 1) * (params2d.cells_number_x - 2) + (i - 1)];
                        else
                            p_next_3D[i][j][k] = 0.0;
                    }
                    if (BOUN_TYPE == 0)
                        p_next_3D[i][j][k] = p_prev_3D[
                                (k) * (params2d.cells_number_x) * (params2d.cells_number_y) +
                                (j) * (params2d.cells_number_x) + (i)];
                    //p_next_3D[i][j][k] = p_next_3D[i][j][k] + ( pow((y[j] + y[j+1]) / 2, 2) + pow((x[i] + x[i+1]) / 2, 2) ) * paramsc.omega * paramsc.omega / 2;
                    //printf("%lf\n", p_next_3D[i][j][k]);
                }
            }
        }
        if (BOUN_TYPE == 1) {
            for (int i = 0; i < params2d.cells_number_x; i++) {
                for (int j = 0; j < params2d.cells_number_y; j++) {
                    for (int k = 0; k < params2d.cells_number_z; k++) {
                        if (i == 0)
                            p_next_3D[i][j][k] = p_next_3D[1][1][1];
                        else if (i == params2d.cells_number_x - 1)
                            p_next_3D[i][j][k] = p_next_3D[1][1][1];
                        else if (j == 0)
                            p_next_3D[i][j][k] = p_next_3D[1][1][1];
                        else if (j == params2d.cells_number_y - 1)
                            p_next_3D[i][j][k] = p_next_3D[1][1][1];
                        else if (k == 0)
                            p_next_3D[i][j][k] = p_next_3D[1][1][1];
                        else if (k == params2d.cells_number_z - 1)
                            p_next_3D[i][j][k] = p_next_3D[1][1][1];

                    }
                }
            }
        }

        calc_third_step(paramsc, params2d, u_next_3D, v_next_3D, w_next_3D, u_first_step_3D, v_first_step_3D, w_first_step_3D, p_next_3D, hx, hy, hz, dt);


double u_max = 0;
int i_max, j_max, k_max;
double kin_prev = 0.0;
        for ( int i = 0; i < params2d.cells_number_x; i++){
            for (int j = 0; j < params2d.cells_number_y; j++){
                for (int k = 0; k < params2d.cells_number_z; k++){
                kin_prev += u_inter_3D[i][j][k] * u_inter_3D[i][j][k] / 2.0 + v_inter_3D[i][j][k] * v_inter_3D[i][j][k] / 2.0 + w_inter_3D[i][j][k] * w_inter_3D[i][j][k] / 2.0;
        	    u_mean[i][j] += u_next_3D[i][j][k] * hz / 2.0 / 3.1415;
        	    v_mean[i][j] += v_next_3D[i][j][k] * hz / 2.0 / 3.1415;
        	    w_mean[i][j] += w_next_3D[i][j][k] * hz / 2.0 / 3.1415;
                    u_inter_3D[i][j][k] = u_prev_3D[i][j][k];
                    v_inter_3D[i][j][k] = v_prev_3D[i][j][k];
                    w_inter_3D[i][j][k] = w_prev_3D[i][j][k];

                }
            }
        }
        
        double kin_z = 0.0;
        for (int i = 0; i < params2d.cells_number_x; i++){
    	    for (int j = 0; j < params2d.cells_number_y; j++){
    		kin_z += ( u_mean[i][j] * u_mean[i][j] / 2.0 + v_mean[i][j] * v_mean[i][j] / 2.0 + w_mean[i][j] * w_mean[i][j] / 2.0 ) * 2.0 * 3.1415 * hx * hy;
    		u_mean[i][j] = 0.0;
    		v_mean[i][j] = 0.0;
    		w_mean[i][j] = 0.0;
    	    }
    	}
        
        for (int j = 0; j < params2d.cells_number_y; j++){
            for (int k = 0; k < params2d.cells_number_z; k++){
                u_inter_3D[params2d.cells_number_x][j][k] = u_prev_3D[params2d.cells_number_x][j][k];
            }
        }
        for (int i = 0; i < params2d.cells_number_x; i++){
            for (int k = 0; k < params2d.cells_number_z; k++){
                v_inter_3D[i][params2d.cells_number_y][k] = v_prev_3D[i][params2d.cells_number_y][k];
            }
        }
        for (int i = 0; i < params2d.cells_number_x; i++){
            for (int j = 0; j < params2d.cells_number_y; j++){
                w_inter_3D[i][j][params2d.cells_number_z] = w_prev_3D[i][j][params2d.cells_number_z];
            }
        }
double kin = 0.0;
//double kin_prev = 0.0;
double work_int = 0.0;
double vort_int = 0.0;
double visc_int = 0.0;
        for ( int i = 0; i < params2d.cells_number_x; i++){
            for (int j = 0; j < params2d.cells_number_y; j++){
                for (int k = 0; k < params2d.cells_number_z; k++){
                    u_prev_3D[i][j][k] = u_next_3D[i][j][k];
                    v_prev_3D[i][j][k] = v_next_3D[i][j][k];
                    w_prev_3D[i][j][k] = w_next_3D[i][j][k];
                    double dwdy, dvdz, dudz, dwdx, dvdx, dudy;
                    //if (i > 0 && j > 0 && k > 0 && i < params2d.cells_number_x-1 && j < params2d.cells_number_y-1 && k < params2d.cells_number_z-1){
                    /*
                    dwdy = ( w_prev_3D[i][j+1][k] - w_prev_3D[i][j-1][k] ) / hy / 2.0;
                    dvdz = ( v_prev_3D[i][j][k+1] - v_prev_3D[i][j][k-1] ) / hz / 2.0;
                    dudz = ( u_prev_3D[i][j][k+1] - u_prev_3D[i][j][k-1] ) / hz / 2.0;
                    dwdx = ( w_prev_3D[i+1][j][k] - w_prev_3D[i-1][j][k] ) / hx / 2.0;
                    dvdx = ( v_prev_3D[i+1][j][k] - v_prev_3D[i-1][j][k] ) / hx / 2.0;
                    dudy = ( u_prev_3D[i][j+1][k] - u_prev_3D[i][j-1][k] ) / hy / 2.0;
                    */
                    kin += u_prev_3D[i][j][k] * u_prev_3D[i][j][k] / 2.0 + v_prev_3D[i][j][k] * v_prev_3D[i][j][k] / 2.0 + w_prev_3D[i][j][k] * w_prev_3D[i][j][k] / 2.0;
                    //kin_prev += u_inter_3D[i][j][k] * u_inter_3D[i][j][k] / 2.0 + v_inter_3D[i][j][k] * v_inter_3D[i][j][k] / 2.0 + w_inter_3D[i][j][k] * w_inter_3D[i][j][k] / 2.0;
                    work_int += u_inter_3D[i][j][k] * force_x_3D[i][j][k] + v_inter_3D[i][j][k] * force_y_3D[i][j][k] + w_inter_3D[i][j][k] * force_z_3D[i][j][k];
                    //vort_int += ( dwdy - dvdz ) * ( dwdy - dvdz ) + ( dudz - dwdx ) * ( dudz - dwdx ) + ( dvdx - dudy ) * ( dvdx - dudy );
                    if (i > 0 && j > 0 && k > 0 && i < params2d.cells_number_x-1 && j < params2d.cells_number_y-1 && k < params2d.cells_number_z-1){
                	visc_int += u_inter_3D[i][j][k] * calc_x_visc_term( paramsc, u_inter_3D, v_inter_3D, w_inter_3D, i, j, k, hx, hy, hz ) +
                	    v_inter_3D[i][j][k] * calc_y_visc_term( paramsc, u_inter_3D, v_inter_3D, w_inter_3D, i, j, k, hx, hy, hz ) +
                	    w_inter_3D[i][j][k] * calc_z_visc_term( paramsc, u_inter_3D, v_inter_3D, w_inter_3D, i, j, k, hx, hy, hz );
                    }
                    //}
                    if ( u_prev_3D[i][j][k] > u_max){
                        u_max = u_prev_3D[i][j][k];
                        i_max = i;
                        j_max = j;
                        k_max = k;
                    }
                }
            }
        }
        for (int j = 0; j < params2d.cells_number_y; j++){
            for (int k = 0; k < params2d.cells_number_z; k++){
                u_prev_3D[params2d.cells_number_x][j][k] = u_next_3D[params2d.cells_number_x][j][k];
                if ( u_prev_3D[params2d.cells_number_x][j][k] > u_max){
                    u_max = u_prev_3D[params2d.cells_number_x][j][k];
                    i_max = params2d.cells_number_x;
                    j_max = j;
                    k_max = k;
                }
            }
        }
        for (int i = 0; i < params2d.cells_number_x; i++){
            for (int k = 0; k < params2d.cells_number_z; k++){
                v_prev_3D[i][params2d.cells_number_y][k] = v_next_3D[i][params2d.cells_number_y][k];
            }
        }
        for (int i = 0; i < params2d.cells_number_x; i++){
            for (int j = 0; j < params2d.cells_number_y; j++){
                w_prev_3D[i][j][params2d.cells_number_z] = w_next_3D[i][j][params2d.cells_number_z];
            }
        }
        printf("%e %d %d %d\n", u_max, i_max, j_max, k_max);
        int i = i_max;
        int j = j_max;
        int k = k_max;
        kin = kin * hy * hx * hz;
        kin_prev = kin_prev * hy * hx * hz;
        work_int = work_int * hx * hy * hz;
        vort_int = vort_int * hx * hy * hz;
        visc_int = visc_int * hx * hy * hz;
        double Energy_balance = ( kin - kin_prev ) / dt / 2.0  - work_int / dt -visc_int;
        printf("%e %e %e\n", u_prev_3D[i][j][k], u_first_step_3D[i][j][k], kin);

        kinetic = fopen( string, "a" );
        fprintf(kinetic,"%e %e %e %e %e %e %e \n",time_mom.curr_t, kin, kin_z, Energy_balance, (kin - kin_prev) / dt / 2.0 , work_int / dt, visc_int);
        fclose(kinetic);

        double vel_aver;
        vel_aver = 0.0;

        vel = fopen( string2, "a" );
        fprintf(vel,"%e %e\n",time_mom.curr_t, vel_aver);
        fclose(vel);
        time_mom.curr_t += dt;
        time_mom.steps_num++;
        // ������ ��������������� ��������� � ����� ����
        printf( "\nStep %d, time %f: dt = %.2e\n", time_mom.steps_num - 1, time_mom.curr_t - dt, dt );

        if ( paramsc.is_output_on_time  ) {
            // ����� �� �������� �������
            if ( time_mom.curr_t >= ( files_counter + 1 ) * time_gap ) {
                write_solution_3d( &paramsc, &params2d, argv[2], x_grid, y_grid, z_grid, p_next_3D, u_prev_3D, v_prev_3D, w_prev_3D, files_counter, output_filename, M2D );
                files_counter++;
            }


        }
        else {
            // ����� �� ���������� ����� �� �������
            if ( time_mom.steps_num >= ( files_counter  ) * steps_gap ) {
                write_solution_3d( &paramsc, &params2d, argv[2], x_grid, y_grid, z_grid, p_next_3D, u_prev_3D, v_prev_3D, w_prev_3D, files_counter, output_filename, M2D );
                files_counter++;
            }


        }
    }


    free(x_grid);
    free(y_grid);
    free(z_grid);
    //free(b);
    free(solution);
    for ( int i = 0; i < params2d.cells_number_x * params2d.cells_number_y; i++ )
        free (matrix[i]);
    free(matrix);
    MPI_Finalize();
    fclose(kinetic);


    printf( "\nNormal finish: total time steps = %d, total time = %f\n\n", time_mom.steps_num, time_mom.curr_t );
    //printf("%lf\n", u_first_step_3D[1][1][1]);

}

/*--------------------------------------------------------------------------
   hypre_FlexGMRESModifyPCAMGExample -

    This is an example (not recommended)
   of how we can modify things about AMG that
   affect the solve phase based on how FlexGMRES is doing...For
   another preconditioner it may make sense to modify the tolerance..

 *--------------------------------------------------------------------------*/

int hypre_FlexGMRESModifyPCAMGExample(void *precond_data, int iterations,
                                      double rel_residual_norm)
{


   if (rel_residual_norm > .1)
   {
      HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver)precond_data, 10);
   }
   else
   {
      HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver)precond_data, 1);
   }


   return 0;
}
