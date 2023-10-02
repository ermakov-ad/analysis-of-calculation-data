// memory.cc
// ������� ��� ������ � �������
// (c) ����� �����, 2013
// ������: 12 ������� 2013 �.

#include "mem.h"

// ��������� ������ ��� ���������� ������ ��������� ���� double
// size - ������ ������� (in) 
// **a - ��������� �� ������ (out)
void get_memory_for_1D_double_array( int size, double **a ) {

    *a = ( double * )malloc( size  * sizeof( double ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_1D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }

}

// ��������� ������ ��� ���������� ������ ��������� ���� int
// size - ������ ������� (in) 
// **a - ��������� �� ������ (out)
void get_memory_for_1D_int_array( int size, int **a ) {

    *a = ( int * )malloc( size * sizeof( int ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_1D_int_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }

}

// ��������� ������ ��� ��������� ������ ��������� ���� double
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// ***a - ��������� �� ������ (out)
void get_memory_for_2D_double_array( int size1, int size2, double ***a ) {
    
    *a = ( double ** )malloc( size1 * sizeof( double * ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_2D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( int i = 0; i < size1; i++ )
        get_memory_for_1D_double_array( size2, &((*a)[i]) );

}

// ��������� ������ ��� ��������� ������ ��������� ���� int
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// ***a - ��������� �� ������ (out)
void get_memory_for_2D_int_array( int size1, int size2, int ***a ) {
    
    *a = ( int ** )malloc( size1 * sizeof( int * ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_2D_int_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( int i = 0; i < size1; i++ )
        get_memory_for_1D_int_array( size2, &((*a)[i]) );

}

// ��������� ������ ��� ���������� ������ ��������� ���� double
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// size3 - ������ ����������� ������� (in)
// ****a - ��������� �� ������ (out)
void get_memory_for_3D_double_array( int size1, int size2, int size3, double ****a ) {
    
    *a = ( double *** )malloc( size1 * sizeof( double ** ) );
    if ( NULL == *a ) {
        printf( "get_memory_for_3D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( int j = 0; j < size1; j++ )
        get_memory_for_2D_double_array( size2, size3, &((*a)[j]) );
}

// ��������� ������ ��� ���������� ������ ��������� ���� double
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// size3 - ������ ����������� ������� (in)
// ****a - ��������� �� ������ (out)
void get_memory_for_3D_int_array( int size1, int size2, int size3, int ****a ) {
    
    *a = ( int *** )malloc( size1 * sizeof( int ** ) );
    if ( NULL == *a ) {
        printf( "get_memory_for_3D_int_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( int j = 0; j < size1; j++ )
        get_memory_for_2D_int_array( size2, size3, &((*a)[j]) );
}