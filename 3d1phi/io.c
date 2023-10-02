// io.cc
// ����� ������� �����/������, �� ��������� �� ����������� ������
// (c) ����� �����, 2018
// ������: 26 ������� 2018 �.

#include "io.h"

// ���������� ��������� � ��������� ����������� ��������������� ������������
// params - ���������� ����������������� ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
void fill_parameters_common( FILE *params, struct ParametersCommon* paramsc ) {

    int dtmp; // ��� ���������� ������������� ����������
    char string[MAX_STRING_SIZE]; // ��� ���������� ��������� ���������� �� �����

    fscanf( params, "%s", string ); // ��������� �������
    fscanf( params, "%s", string ); // ���������, ������������, ��� ����������� ����� �����, �� ��������� �� ����������� ������

    
    // ��� �� �������
    fscanf( params, "%s", string ); // ��������� �������
    fscanf( params, "%s %d", string, &dtmp ); // ��� - ���������� ��� ���
    switch ( dtmp ) {
        case 0:
            paramsc->constant_time_step = false;
            break;
        case 1:
            paramsc->constant_time_step = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> constant_time_step should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->cfl) ); // ����� cfl
    if ( !paramsc->constant_time_step ) {
        if ( paramsc->cfl <= 0.0 ) {
            printf( "\nfill_parameters_struct -> cfl number should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf( params, "%s %lf", string, &(paramsc->dt) ); // ���������� ��� �� �������
    if ( paramsc->constant_time_step ) {
        if ( paramsc->dt <= 0.0 ) {
            printf( "\nfill_parameters_struct -> dt should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }


    fscanf( params, "%s %d", string, &(paramsc->approximation_order) ); // ������� ������������� ������
    if ( paramsc->approximation_order != 1 && paramsc->approximation_order != 2 ) {
        printf( "\nfill_parameters_struct -> approximation_order should be 1 or 2\n\n" );
        exit( EXIT_FAILURE );
    }
            
    
    fscanf( params, "%s %lf", string, &(paramsc->mu) ); // ����������� ������������ �������� � ����
    if ( paramsc->mu < 0.0 ) {
        printf( "\nfill_parameters_struct -> mu_g should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
    
    // ����� �������
    fscanf( params, "%s", string ); // ��������� �������
    fscanf( params, "%s %d", string, &dtmp );
    switch ( dtmp ) {
        case 1:
            paramsc->is_debug = true; // ������ � ���������� ��������� �������� � ������������ ������� ������ ���������
            break;
        case 0:
            paramsc->is_debug = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_debug should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }

    // �����
    fscanf( params, "%s", string ); // ��������� �������
    fscanf( params, "%s %d", string, &dtmp );
    switch ( dtmp ) {
        case 1:
            paramsc->is_output_on_time = true; // �������� ������ ����������� - ���������� ��������� ������� �������
            break;
        case 0:
            paramsc->is_output_on_time = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_output_on_time should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->stop_time) ); // ������ ������� ��������� �������
    printf("%lf\n",paramsc->stop_time);
    if ( paramsc->is_output_on_time ) {
        if ( paramsc->stop_time <= 0.0 ) {
            printf( "\nfill_parameters_struct -> stop_time should be a positive value.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf( params, "%s %d", string, &(paramsc->stop_steps_number) ); // ��������� ���������� ����� �� �������
    if ( !paramsc->is_output_on_time ) {
        if ( paramsc->stop_steps_number <= 0.0 ) {
            printf( "\nfill_parameters_struct -> stop_steps_number should be a positive value.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf( params, "%s %d", string,  &(paramsc->output_number) ); // ���������� ������ � �������������� ������������
    if ( paramsc->output_number <= 0 ) {
        printf( "\nfill_parameters_struct -> output_number should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %d", string,  &dtmp ); // �������� �� ��������� ��� ������������ � Tecplot?
    switch ( dtmp ) {
        case 1:
            paramsc->is_tecplot_header = true;
            break;
        case 0:
            paramsc->is_tecplot_header = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_tecplot_header should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }

    fscanf( params, "%s %lf", string, &(paramsc->f0) );
    fscanf( params, "%s %lf", string, &(paramsc->omega) );
    // ����� ���������
    fscanf( params, "%s", string ); // ��������� �������
    fscanf( params, "%s %lf", string, &(paramsc->eps_general) ); // ���������� ����� ����� ��� ��������� ������������ �����,
                                                                                    // ������������� � �������������� ��������, ��������
                                                                                    // ���������� ������������ ��������� (���� ���� �� ��������� �����)
    if ( paramsc->eps_general <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_general should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->eps_ludcmp) ); // ����� ����� � ������� ludcmp (math_utilc.cc),
                                                                                   // ��������� ��� ���� � ������������ ������
    if ( paramsc->eps_ludcmp <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_ludcmp should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->eps_decouple) ); // ����������� ���������� ������� �� �������� ���� ����� � ������,
                                                                                     // ��� ������� ����������� ���������������� �����, ����� - �����������
                                                                                     // ������������ ����� ��������� ��� ���
    if ( paramsc->eps_decouple <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_decouple should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->eps_thin_layer) ); // �������� ������� ������� ���������� �������������� ���������
                                                                                       // "������� ����" ��� ��������
    if ( paramsc->eps_thin_layer <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_thin_layer should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->eps_contact) ); // ����� �������� ��� ������� �������� � ������ �������������
                                                                                    // ����������� �������
    if ( paramsc->eps_contact <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_contact should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->eps_disp_abs) ); // ����� �������� �������� ���� ����������
                                                                                     // ����, ������� � ������� ���������, ��� ���
                                                                                     // �����������
    if ( paramsc->eps_disp_abs <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_disp_abs should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf( params, "%s %lf", string, &(paramsc->eps_cut_out) ); // ����� �������� �������� ���� ����������
                                                                                    // ����, ������� � ������� ��� ������
                                                                                    // ����������� ��������� ���������� ����
                                                                                    // �������� ������� ������� ����������
    if ( paramsc->eps_cut_out <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_cut_out should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }

}

// ������ ����� � ����������� ��� ����������� ����������� �������
// files_directory - ����������, � ������� ���������� ������ (in)
// time_moment - ������ � ����������� � ������� ������� ��� ����� �����, ��� �������
// � ��������� ��� ��� ������� ���� � �������������� ������������ (in)
// filename - ��� ���������� ����������� ����� (in)
void write_restart_info( char *files_directory, char *time_moment, char *filename ) {

    FILE *restart_info_file; // ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������
    char string[MAX_STRING_SIZE]; // ��� ������������ ����� �����

    // ����������, ����������� ��� ��������, ��������� � ����� restart_info.dat
    strcpy( string, files_directory );
    strcat( string, "\\restart_info.dat" );
    if ( (restart_info_file = fopen(  string, "wt" ) ) == NULL ) {
        printf( "\nwrite_restart_info -> can't open file %s for writing\n\n", string );
        exit( EXIT_FAILURE );
    }

    // ������ ���������� � ������� ������� � �� ����� ���������� ����� � �������������� ������������
    fprintf( restart_info_file, "%s %s", time_moment, filename );

    fclose( restart_info_file );

}	



/* ���������� ����� � �������� ��� ���������� ����������� ����������� ����������� ������� 

   params - ��������� � ����������� ��������������� ������������ (in)
   restart_info_file - ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������ (in)

   **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out) */
void read_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, FILE *restart_info_file, double **initial_solution ) {

    // char string[MAX_STRING_SIZE];   /* ��� ���������� ��������� ���������� �� ����� */

}

// ������� ������� ���� � ����������� � �������� �, � ������ ������, ��������� �����������
// ����� � ���������� ����������� ���������������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in) 
// files_directory - ����������, � ������� ���������� ������ (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// file_to_read - ���������� ����� � ��������������� ��� �������� (out)
void try_to_restart( struct ParametersCommon *paramsc, char *files_directory, struct TimeMoment *time_mom, FILE *file_to_read ) {


    printf( "Data from is not loading.\n");

}
