// io_2d.cc
// ������� �����/������ ��� ���������� ����
// (c) ����� �����, 2018
// ������: 27 ������� 2018 �.

#include "io_2d2phc.h"
#include "io.h"

// ���������� ��������� � ����������� ��������� ������
// parameters - ���������� ����������������� ����� ������
// params2d - ��������� � ����������� ��������� ������
void fill_parameters_2d( FILE *parameters, struct Parameters2d* params2d ) {
    
    int dtmp; // ��� ���������� ������������� ����������
    char string[MAX_STRING_SIZE]; // ��� ���������� ��������� ���������� �� �����

    fscanf( parameters, "%s", string ); // ���������, �������������, ��� ����������� ����� ����������, ����������� ���������� ������

    // ��������� �����
    fscanf( parameters, "%s", string ); // ��������� �������                               (S) ---GRID---
    fscanf( parameters, "%s %lf", string, &(params2d->left_boundary_x) ); // ���������� ������ ����� ��������� �������
    fscanf( parameters, "%s %lf", string, &(params2d->right_boundary_x) ); // ���������� ������� ����� ��������� �������
    fscanf( parameters, "%s %lf", string, &(params2d->down_boundary_y) ); // ���������� ������� ����� ��������� �������
    fscanf( parameters, "%s %lf", string, &(params2d->up_boundary_y) ); // ���������� �������� ����� ��������� �������
    fscanf( parameters, "%s %lf", string, &(params2d->rear_boundary_z) ); // ���������� ������� ����� ��������� �������
    fscanf( parameters, "%s %lf", string, &(params2d->front_boundary_z) ); // ���������� �������� ����� ��������� �������
    fscanf( parameters, "%s %d", string, &(params2d->cells_number_x) ); // ���������� ����� x
    fscanf( parameters, "%s %d", string, &(params2d->cells_number_y) ); // ���������� ����� �� y
    fscanf( parameters, "%s %d", string, &(params2d->cells_number_z) ); // ���������� ����� �� z
	//printf("%d %d\n",params2d->cells_number_y,params2d->cells_number_x);

    if ( params2d->cells_number_x <= 1 ) 
    {
        printf( "\nread_parameters -> cells_number_x should be a positive value %d\n\n", params2d->cells_number_x );
        exit( EXIT_FAILURE );
    }

    if ( params2d->cells_number_y <= 1 ) 
    {
        printf( "\nread_parameters -> cells_number_y should be a positive value %d\n\n", params2d->cells_number_y );
        exit( EXIT_FAILURE );
    }

    
    // ��������� �������
    fscanf( parameters, "%s", string ); // ��������� �������
    
    // ������������� ����� � ��������� �������
    fscanf( parameters, "%s %d", string,  &(params2d->x_blocks_number) ); // ���������� ������ �� x
    fscanf( parameters, "%s %d", string,  &(params2d->y_blocks_number) ); // ���������� ������ �� y
    fscanf( parameters, "%s %d", string,  &(params2d->z_blocks_number) ); // ���������� ������ �� z

    if ( params2d->x_blocks_number <= 0 ) 
    {
        printf( "\nread_parameters -> x_blocks_number should be positive.\n\n" );
        exit( EXIT_FAILURE );
    }
    if ( params2d->y_blocks_number <= 0 ) 
    {
        printf( "\nread_parameters -> y_blocks_number should be positive.\n\n" );
        exit( EXIT_FAILURE );
    }

    for ( int i_block = 0; i_block < params2d->x_blocks_number; i_block++ ){ // ���� �� ������ ��������� �������
        for ( int j_block = 0; j_block < params2d->y_blocks_number; j_block++ ){ 
            for ( int k_block = 0; k_block < params2d->z_blocks_number; k_block++ ){ 
         
                fscanf( parameters, "%s", string); // ��������� - ����� �����
                fscanf( parameters, "%s %d", string, &(params2d->x_begin[i_block][j_block][k_block]) ); // ����� ������-������ �� x
                if ( params2d->x_begin[i_block][j_block][k_block] < 0 ) 
                { 
                    printf( "\nread_parameters -> params->x_begin[%d][%d][%d] should be non-negative .\n\n", i_block, j_block, k_block );
                    exit( EXIT_FAILURE );
                }
                fscanf( parameters, "%s %d", string,  &(params2d->x_end[i_block][j_block][k_block]) ); // ����� ������-����� �� x
                if ( params2d->x_end[i_block][j_block][k_block] < 0 ) 
                {
                    printf( "\nread_parameters -> params->x_end[%d][%d][%d] should be non-negative .\n\n", i_block, j_block, k_block );
                    exit( EXIT_FAILURE );
                }

                fscanf( parameters, "%s %d", string,  &(params2d->y_begin[i_block][j_block][k_block]) ); // ����� ������-������ �� y
                if ( params2d->y_begin[i_block][j_block][k_block] < 0 ) 
                {
                    printf( "\nread_parameters -> params->y_begin[%d][%d][%d] should be non-negative .\n\n", i_block, j_block, k_block );
                    exit( EXIT_FAILURE );
                }
                fscanf( parameters, "%s %d", string,  &(params2d->y_end[i_block][j_block][k_block]) ); // ����� ������-����� �� y
                if ( params2d->y_end[i_block][j_block][k_block] < 0 ) 
                {
                    printf( "\nread_parameters -> params->y_end[%d][%d][%d] should be non-negative .\n\n", i_block, j_block, k_block );
                    exit( EXIT_FAILURE );
                }

                fscanf( parameters, "%s %d", string,  &(params2d->z_begin[i_block][j_block][k_block]) ); // ����� ������-������ �� y
                if ( params2d->z_begin[i_block][j_block][k_block] < 0 ) 
                {
                    printf( "\nread_parameters -> params->y_begin[%d][%d][%d] should be non-negative .\n\n", i_block, j_block, k_block );
                    exit( EXIT_FAILURE );
                }
                fscanf( parameters, "%s %d", string,  &(params2d->z_end[i_block][j_block][k_block]) ); // ����� ������-����� �� y
                if ( params2d->z_end[i_block][j_block][k_block] < 0 ) 
                {
                    printf( "\nread_parameters -> params->y_end[%d][%d][%d] should be non-negative .\n\n", i_block, j_block, k_block );
                    exit( EXIT_FAILURE );
                }

                fscanf( parameters, "%s %d", string,  &(params2d->block_type[i_block][j_block][k_block]) );
                if ( params2d->block_type[i_block][j_block][k_block] != INNER && params2d->block_type[i_block][j_block][k_block] != OUTER ) 
                {
                    printf( "\nread_parameters -> params2d->block_type[%d][%d] should be 1 or 2 .\n\n", i_block, j_block );
                    exit( EXIT_FAILURE );
                }

                // ������ ����������� ���������� � �����
                fscanf( parameters, "%s %lf", string,  &(params2d->block_values_p[i_block][j_block][k_block]) );
                fscanf( parameters, "%s %lf", string,  &(params2d->block_values_u[i_block][j_block][k_block]) );
                fscanf( parameters, "%s %lf", string,  &(params2d->block_values_v[i_block][j_block][k_block]) );
                fscanf( parameters, "%s %lf", string,  &(params2d->block_values_w[i_block][j_block][k_block]) );

            }
        }   
    }
}

// ���������� ��������� � ����������� ��������� ������ � ���������� ���������� �������� �����, � ����� ��������� ���������� ��������� ������.
// �������� ������������ ������� ����������.
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params2d - ��������� � ����������� ��������� ������
// input_file_directory - ����������, � ������� ��������� ���� � ����������� ������
// output_file_directory - ����������, � ������� ������ ���� �������� ����� � ������������
// file_num - ������� ����� ����� Parameters1d.dat ��� ����������
void fill_parameters_struct_2d( struct ParametersCommon* paramsc, struct Parameters2d* params2d, const char* input_file_directory, const char* output_file_directory,
                               const int file_num ){
    
    FILE *params; // ���������� ����������������� ����� ������
    char string[MAX_STRING_SIZE]; // ��� ���������� ��������� ���������� �� �����
            
    // ��� ��������� ������ ��������� � ����� parameters.dat
    strcpy( string, input_file_directory );
    strcat( string, "//parameters" );
    if ( file_num != 0 ) { // ���� 0, �� ������ �� ���������� ��� ���������� ��������������� � ���������� ���������
        char tmp_str[MAX_STRING_SIZE]; // ��������� ���������� ��� �������� �������� ������
        sprintf( tmp_str, "%d", file_num );
        strcat( string, tmp_str );
    }
    strcat( string, ".dat" );
    if ( ( params = fopen( string, "rt" ) ) == NULL ) {
        printf( "\nfill_parameters_struct_2d -> can't open file %s for reading\n\n", string );
        exit( EXIT_FAILURE );
    }

    // ��������� �������� ����������, � ������� ��������� ���� � ����������� ������, � params
    strcpy( paramsc->input_file_directory, input_file_directory );
    
    // ��������� �������� ����������, ���� ����� �������� ���������� ��������, � params
    strcpy( paramsc->output_file_directory, output_file_directory );
    
    // ���������� ��������� �����
    fscanf( params, "%s", string);

    // ���������� ����� ���������� �������
    fill_parameters_common( params, paramsc );

    // ���������� ���������� ���������� ������
    fill_parameters_2d( params, params2d );

    fclose( params );
}

/* ������ ������� � ����

   params - ��������� � ����������� ��������������� ������������ (in)
   output_directory - ����������, ���� ����� ������� ���� � ������ �������� (in)
   cells_number - ����� ����� (in) 
   *xc - ������ ��������� ������� ����� (in)
   **v_ncons_x - ��������� ������ ����������� ���������� � ������� ���� ����� ����� �� x (in)
   **v_ncons_y - ��������� ������ ����������� ���������� � ������� ���� ����� ����� �� y (in)   // (S) ����� �� �������� ; �� ����� ���� ������  v_ncons
   file_number - ����� ����� � �������������� ������������ (in)
   output_filename - ��� �������� ����� � �������������� ������������ (f_out)
   n - �������� ������ �������� (in)    */
void write_solution_3d( struct ParametersCommon *paramsc, struct Parameters2d *params2d, char *output_directory, double *x, double *y, double *z, double ***p, double ***u, double ***v, double ***w, int file_number, char *output_filename, int n ) 
{
    FILE *f_out; // ���������� ����� � ������������
    int i, j, k; // �������
    char tmp_str[MAX_STRING_SIZE]; // ��������� ���������� ��� �������� �������� ������
        
    // ������������ ����� ��������� ����� � ��������� ��������
    strcpy( output_filename,  output_directory );
    strcat( output_filename, "//numerical_solution_" );
    // ���������� ����������� ������
    if ( file_number < TEN_TO_FIFTH )
        strcat( output_filename,  "0" );
    if ( file_number < TEN_TO_FOURTH )
        strcat( output_filename,  "0" );
    if ( file_number < TEN_TO_THIRD )
        strcat( output_filename,  "0" );
    if ( file_number < TEN_TO_SECOND )
        strcat( output_filename,"0" );
    if ( file_number < TEN )
        strcat( output_filename, "0" );
    sprintf( tmp_str, "%d", file_number );
    // ������ ������ �����
    strcat( output_filename, tmp_str );
    strcat( output_filename, ".plt" );
    
    // �������� ����� ��� ������ � �������� �������
    f_out = fopen( output_filename, "wb" );

    if ( NULL == f_out ) 
    {
        printf( "\nwrite_solution -> can't open file %s for writing\n", output_filename );
    }

    // ������ ��������� ��� Tecplot � ������ �������������
    if ( paramsc->is_tecplot_header )
        write_solution_tecplot_header_3d( paramsc, params2d, f_out, file_number );

    // ������ ������� �����
    for ( k = 0; k < params2d->cells_number_z + 1; k++ )
    {
    	for ( j = 0; j < params2d->cells_number_y + 1; j++ )
    	{
    		for ( i = 0; i < params2d->cells_number_x + 1; i++ )
    		{

                fwrite( &( x[i] ), sizeof( double ), 1, f_out );
            }
        }
    }

    // ������ ������� �����
    for ( k = 0; k < params2d->cells_number_z + 1; k++ )
    {
    	for ( j = 0; j < params2d->cells_number_y + 1; j++ )
    	{
    		for ( i = 0; i < params2d->cells_number_x + 1; i++ )
    		{

                fwrite( &( y[j] ), sizeof( double ), 1, f_out );
            }
        }
    }
    // ������ �������� �����
    for ( k = 0; k < params2d->cells_number_z + 1; k++ )
    {
    	for ( j = 0; j < params2d->cells_number_y + 1; j++ )
    	{
    		for ( i = 0; i < params2d->cells_number_x + 1; i++ )
    		{

                fwrite( &( z[k] ), sizeof( double ), 1, f_out );
            }
        }
    }
	double u_center;
	double v_center;
        double w_center;
    // ������ ���������� � ������
	for ( k = 0; k < params2d->cells_number_z; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x; i++ ) 
            {

                    fwrite( &( p[i][j][k] ), sizeof( double ), 1, f_out );
                }

            }
    }
	for ( k = 0; k < params2d->cells_number_z; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x; i++ ) 
            {

		    u_center = ( u[i+1][j][k] + u[i][j][k] ) / 2;
                    //printf( "\n%e\n", u_center);
		    fwrite( &( u_center ), sizeof( double ), 1, f_out );
                }
            }
	}
	for ( k = 0; k < params2d->cells_number_z; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x; i++ ) 
            {

				v_center = ( v[i][j+1][k] + v[i][j][k] ) / 2;
				fwrite( &( v_center ), sizeof( double ), 1, f_out );
                }
            }
	}
	for ( k = 0; k < params2d->cells_number_z; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x; i++ ) 
            {

				w_center = ( w[i][j][k+1] + w[i][j][k] ) / 2;
				fwrite( &( w_center ), sizeof( double ), 1, f_out );
                }
            }
	}

    fclose( f_out );
}

/* ������ � ���� � ������������ ������������ ���������� ��� Tecplot

   params - ��������� � ����������� ��������������� ������������ (in) 
   file_to_write - ���������� �����, ��������� ��� ������ (in)
   file_number - ����� ����� � �������������� ������������ (in) */
void write_solution_tecplot_header_3d( struct ParametersCommon *paramsc, struct Parameters2d *params2d, FILE *f_out, int file_number )
{
    char tmp_str[MAX_STRING_SIZE]; // ��������� ���������� ��� ������ ������
    double time_gap; // ����� ����� �������� ������ � �������������� ������������, ����� ����� ��� params->is_output_on_time = 1
    int steps_gap; // ���������� ����� �� ������� ����� �������� ������ � �������������� ������������, ����� ����� ��� params->is_output_on_time = 0
    char version[] = "#!TDV102";
    char title[] = "3D Incompressible solver ";
    char vars[] = "X Y Z P U V W ";
    char zone[MAX_STRING_SIZE];
    int integer;
    float real32;
    int i;

    // ����� ������
    fwrite( &version, sizeof( char ), 8, f_out );
    
    // ������������� �������
    integer = 1;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // ���������
    for ( unsigned i = 0; i < strlen( title ); i++ ) 
    {
        integer = title[i];
        fwrite( &integer, sizeof( int ), 1, f_out );
    }
    integer = 0;
    fwrite( &integer, sizeof(int), 1, f_out );

    // ���������� ����������
    integer = M2D+3;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // ����� ����������
    for ( unsigned i = 0; i < strlen( vars ); i++ ) 
    {
        integer = vars[i];
        if ( vars[i] != ' ' ) 
        {
            fwrite( &integer, sizeof( int ), 1, f_out );
        }
        else 
        {
            integer = 0;
            fwrite( &integer, sizeof( int ), 1, f_out );
        }
    }
    
    // ������ ����
    real32 = 299.0;
    fwrite( &real32, sizeof( float ), 1, f_out );

    // ������������ �������� ����
    strcpy( zone,  "ZONE T=" );
    if ( paramsc->is_output_on_time ) 
    {
        // ������ ������ �� ���������� ��������� ������� �������
        time_gap = paramsc->stop_time / paramsc->output_number;
        sprintf( tmp_str, "%f", time_gap * file_number );
        strcat( tmp_str, " time units" );
    }
    else 
    {
        // ������ ������ �� ���������� ��������� ���������� �����
        steps_gap = paramsc->stop_steps_number / paramsc->output_number;
        sprintf( tmp_str, "%d", steps_gap * file_number );
        strcat( tmp_str, " time steps" );
    }
    strcat( zone, tmp_str );

    // ������ �������� ����
    for ( unsigned i = 0; i < strlen( zone ); i++ ) 
    {
        integer = zone[i];
        fwrite( &integer, sizeof( int ), 1, f_out );
    }
    integer = 0;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // ���� ����: -1 - ���� ���������� ��� Tecplot
    integer = -1;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // ��� ����: 0 - �������������
    integer = 0;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // �������������� ������: 0 - �� ������, 1 - �� ������
    integer = 0;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // ��������� ����������: 1 - ����������, 0 - �� ����������
    integer = 1;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // 0 - � �����; 1 - � ������� �����
    integer = 0;
    // �������� ����
    fwrite( &integer, sizeof( int ), 1, f_out );
    // �������� ����
    fwrite( &integer, sizeof( int ), 1, f_out );
    // ��������� ���� 
    fwrite( &integer, sizeof( int ), 1, f_out );
    integer = 1;
    // ���������� � ������
    for ( i = 0; i < M2D; i++ ) 
    {
        fwrite( &integer, sizeof( int ), 1, f_out );
    }

    // "number of user defined face neighbor connections": 0
    integer = 0;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // imax: nx
    integer = params2d->cells_number_x + 1;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // jmax: ny
    integer = params2d->cells_number_y + 1;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // kmax: nz
    integer = params2d->cells_number_z + 1;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // kmax: 1
    //integer = 1;
    //fwrite( &integer, sizeof( int ), 1, f_out );

    // 0 - "no more auxiliar name/value pairs"
    integer = 0;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // ������ ����� ���������� � �������
    real32 = 357.0;
    fwrite( &real32, sizeof( float ), 1, f_out );

    // ������ ����
    real32 = 299.0;
    fwrite( &real32, sizeof( float ), 1, f_out );

    // ������ ����������
    integer = 2; // 2 - double
    for ( i = 0; i < M2D + 3; i++ ) 
    {
        fwrite( &integer, sizeof( int ), 1, f_out );
    }
	
    // "has variable sharing": 0 - no
    integer = 0;
    fwrite( &integer, sizeof( int ), 1, f_out );

    // "zone number to share connectivity list with":  -1 - no sharing
    integer = -1;
    fwrite( &integer, sizeof( int ), 1, f_out );

}



/* ���������� ����� � �������� ��� ���������� ����������� ����������� ����������� ������� 

   params - ��������� � ����������� ��������������� ������������ (in)
   restart_info_file - ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������ (in)

   **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (f_out) */
void read_solution( struct Parameters2d *params2d, struct ParametersCommon* paramsc, int file_number, double ***u, double ***v, double ***w )
{

    char tmp_str[MAX_STRING_SIZE]; // ��������� ���������� ��� ������ ������
    double time_gap; // ����� ����� �������� ������ � �������������� ������������, ����� ����� ��� params->is_output_on_time = 1
    int steps_gap; // ���������� ����� �� ������� ����� �������� ������ � �������������� ������������, ����� ����� ��� params->is_output_on_time = 0
    char version[] = "#!TDV102";
    char title[] = "3D Incompressible solver ";
    char vars[] = "X Y Z P U V W ";
    char zone[MAX_STRING_SIZE];
    int integer;
    float real32;
    FILE *f_out; // ���������� ����� � ������������
     int i, j, k; // �������
     char output_filename[MAX_STRING_SIZE];
     // ������������ ����� ��������� ����� � ��������� ��������
     strcpy( output_filename, "./numerical_solution_" );
     // ���������� ����������� ������
     if ( file_number < TEN_TO_FIFTH )
         strcat( output_filename,  "0" );
     if ( file_number < TEN_TO_FOURTH )
         strcat( output_filename,  "0" );
     if ( file_number < TEN_TO_THIRD )
         strcat( output_filename,  "0" );
     if ( file_number < TEN_TO_SECOND )
         strcat( output_filename,"0" );
     if ( file_number < TEN )
         strcat( output_filename, "0" );
     sprintf( tmp_str, "%d", file_number );
     // ������ ������ �����
     strcat( output_filename, tmp_str );
     strcat( output_filename, ".plt" );
    f_out = fopen( output_filename, "rb" );

    if ( NULL == f_out )
    {
        printf( "\nwrite_solution -> can't open file %s for read\n", output_filename );
    }
    else
	printf("\n opened %s for read\n", output_filename);

    // ����� ������
    fread( &version, sizeof( char ), 8, f_out );
    
    // ������������� �������
    integer = 1;
    fread( &integer, sizeof( int ), 1, f_out );

    // ���������
    for ( unsigned i = 0; i < strlen( title ); i++ ) 
    {
        integer = title[i];
        fread( &integer, sizeof( int ), 1, f_out );
    }
    integer = 0;
    fread( &integer, sizeof(int), 1, f_out );

    // ���������� ����������
    integer = M2D+3;
    fread( &integer, sizeof( int ), 1, f_out );

    // ����� ����������
    for ( unsigned i = 0; i < strlen( vars ); i++ ) 
    {
        integer = vars[i];
        if ( vars[i] != ' ' ) 
        {
            fread( &integer, sizeof( int ), 1, f_out );
        }
        else 
        {
            integer = 0;
            fread( &integer, sizeof( int ), 1, f_out );
        }
    }
    
    // ������ ����
    real32 = 299.0;
    fread( &real32, sizeof( float ), 1, f_out );
    printf("%f\n", real32);
    // ������������ �������� ����
    strcpy( zone,  "ZONE T=" );
    if ( paramsc->is_output_on_time ) 
    {
        // ������ ������ �� ���������� ��������� ������� �������
        time_gap = paramsc->stop_time / paramsc->output_number;
        sprintf( tmp_str, "%f", time_gap * file_number );
        strcat( tmp_str, " time units" );
    }
    else 
    {
        // ������ ������ �� ���������� ��������� ���������� �����
        steps_gap = paramsc->stop_steps_number / paramsc->output_number;
        sprintf( tmp_str, "%d", steps_gap * file_number );
        strcat( tmp_str, " time steps" );
    }
    strcat( zone, tmp_str );

    // ������ �������� ����
    for ( unsigned i = 0; i < strlen( zone ); i++ ) 
    {
        integer = zone[i];
        fread( &integer, sizeof( int ), 1, f_out );
    }
    integer = 0;
    fread( &integer, sizeof( int ), 1, f_out );

    // ���� ����: -1 - ���� ���������� ��� Tecplot
    integer = -1;
    fread( &integer, sizeof( int ), 1, f_out );

    // ��� ����: 0 - �������������
    integer = 0;
    fread( &integer, sizeof( int ), 1, f_out );

    // �������������� ������: 0 - �� ������, 1 - �� ������
    integer = 0;
    fread( &integer, sizeof( int ), 1, f_out );

    // ��������� ����������: 1 - ����������, 0 - �� ����������
    integer = 1;
    fread( &integer, sizeof( int ), 1, f_out );

    // 0 - � �����; 1 - � ������� �����
    integer = 0;
    // �������� ����
    fread( &integer, sizeof( int ), 1, f_out );
    // �������� ����
    fread( &integer, sizeof( int ), 1, f_out );
    // ��������� ���� 
    fread( &integer, sizeof( int ), 1, f_out );
    integer = 1;
    // ���������� � ������
    for ( i = 0; i < M2D; i++ ) 
    {
        fread( &integer, sizeof( int ), 1, f_out );
    }

    // "number of user defined face neighbor connections": 0
    integer = 0;
    fread( &integer, sizeof( int ), 1, f_out );

    // imax: nx
    integer = params2d->cells_number_x + 1;
    fread( &integer, sizeof( int ), 1, f_out );

    // jmax: ny
    integer = params2d->cells_number_y + 1;
    fread( &integer, sizeof( int ), 1, f_out );

    // kmax: nz
    integer = params2d->cells_number_z + 1;
    fread( &integer, sizeof( int ), 1, f_out );
    //fread( &integer, sizeof( int ), 1, f_out );
    printf("%d\n", integer);

    // kmax: 1
    //integer = 1;
    //fread( &integer, sizeof( int ), 1, f_out );

    // 0 - "no more auxiliar name/value pairs"
    integer = 0;
    fread( &integer, sizeof( int ), 1, f_out );

    // ������ ����� ���������� � �������
    real32 = 357.0;
    fread( &real32, sizeof( float ), 1, f_out );
    printf("%f\n", real32);

    // ������ ����
    real32 = 299.0;
    fread( &real32, sizeof( float ), 1, f_out );

    // ������ ����������
    integer = 2; // 2 - double
    for ( i = 0; i < M2D + 3; i++ ) 
    {
        fread( &integer, sizeof( int ), 1, f_out );
    }
	
    // "has variable sharing": 0 - no
    integer = 0;
    fread( &integer, sizeof( int ), 1, f_out );

    // "zone number to share connectivity list with":  -1 - no sharing
    integer = -1;
    fread( &integer, sizeof( int ), 1, f_out );

    double x, y, z;
    // ������ ������� �����
    for ( k = 0; k < params2d->cells_number_z + 1; k++ )
    {
    	for ( j = 0; j < params2d->cells_number_y + 1; j++ )
    	{
    		for ( i = 0; i < params2d->cells_number_x + 1; i++ )
    		{

                fread( &( x ), sizeof( double ), 1, f_out );
            }
        }
    }

    // ������ ������� �����
    for ( k = 0; k < params2d->cells_number_z + 1; k++ )
    {
    	for ( j = 0; j < params2d->cells_number_y + 1; j++ )
    	{
    		for ( i = 0; i < params2d->cells_number_x + 1; i++ )
    		{

                fread( &( y ), sizeof( double ), 1, f_out );
            }
        }
    }
    // ������ �������� �����
    for ( k = 0; k < params2d->cells_number_z + 1; k++ )
    {
    	for ( j = 0; j < params2d->cells_number_y + 1; j++ )
    	{
    		for ( i = 0; i < params2d->cells_number_x + 1; i++ )
    		{

                fread( &( z ), sizeof( double ), 1, f_out );
                //printf("%e\n", z);
            }
        }
    }
	double u_center;
	double v_center;
    double w_center;
    double p;
    // ������ ���������� � ������
	for ( k = 0; k < params2d->cells_number_z; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x; i++ )
            {

                    fread( &( p), sizeof( double ), 1, f_out );
                }

            }
    }
	for ( k = 0; k < params2d->cells_number_z; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x + 1; i++ )
            {
            	if ( i < params2d->cells_number_x)
            		fread( &( u[i][j][k] ), sizeof( double ), 1, f_out );
            	else
            		u[i][j][k] = u[i-1][j][k];
            }
        }
	}
	for ( k = 0; k < params2d->cells_number_z; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y + 1; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x; i++ )
            {
            	if ( j < params2d->cells_number_y)
            		fread( &( v[i][j][k] ), sizeof( double ), 1, f_out );
            	else
            		v[i][j][k] = v[i][j-1][k];
            }
        }
	}
	for ( k = 0; k < params2d->cells_number_z + 1; k++ )
	{
        for ( j = 0; j < params2d->cells_number_y; j++ )
        {
            for ( i = 0; i < params2d->cells_number_x; i++ )
            {
            	if ( k < params2d->cells_number_z)
            		fread( &( w[i][j][k] ), sizeof( double ), 1, f_out );
            	else
            		w[i][j][k] = w[i][j][k-1];
            	//printf("%e\n", w[i][j][k]);
            }
        }
	}

    fclose( f_out );

}
