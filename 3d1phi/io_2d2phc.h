// io_2d.h
// ������� �����/������ ��� ���������� ����
// (c) ����� �����, 2018
// ������: 27 ������� 2018 �.

#ifndef __IO_2D_H_
#define __IO_2D_H_

#include "ex5.h"

// ��������� ��� ���������� ����� ����� ����������� ������
#define TEN_TO_FIFTH    100000
#define TEN_TO_FOURTH   10000
#define TEN_TO_THIRD    1000
#define TEN_TO_SECOND   100
#define TEN             10

// ���������� ��������� � ����������� ��������� ������
// parameters - ���������� ����������������� ����� ������
// params2d - ��������� � ����������� ��������� ������
void fill_parameters_2d( FILE *parameters, struct Parameters2d* params2d );

// ���������� ��������� � ����������� ��������� ������ � ���������� ���������� �������� �����, � ����� ��������� ���������� ��������� ������.
// �������� ������������ ������� ����������.
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params2d - ��������� � ����������� ��������� ������
// input_file_directory - ����������, � ������� ��������� ���� � ����������� ������
// output_file_directory - ����������, � ������� ������ ���� �������� ����� � ������������
// file_num - ������� ����� ����� Parameters1d.dat ��� ����������
void fill_parameters_struct_2d( struct ParametersCommon* paramsc, struct Parameters2d* params2d, const char* input_file_directory, const char* output_file_directory,
                               const int file_num );

/* ������ ������� � ����

   params - ��������� � ����������� ��������������� ������������ (in)
   output_directory - ����������, ���� ����� ������� ���� � ������ �������� (in)
   cells_number - ����� ����� (in) 
   *xc - ������ ��������� ������� ����� (in)
   **v_ncons - ��������� ������ ����������� ���������� � ������� ���� ����� ����� (in)
   file_number - ����� ����� � �������������� ������������ (in)
   
   output_filename - ��� �������� ����� � �������������� ������������ (out) */
void write_solution_3d( struct ParametersCommon *paramsc, struct Parameters2d *params2d, char *output_directory, double *x, double *y, double *z, double ***p, double ***u, double ***v, double ***w, int file_number, char *output_filename, int n ) ;

/* ������ � ���� � ������������ ������������ ���������� ��� Tecplot

   params - ��������� � ����������� ��������������� ������������ (in) 
   file_to_write - ���������� �����, ��������� ��� ������ (in)
   file_number - ����� ����� � �������������� ������������ (in) */
void write_solution_tecplot_header_3d( struct ParametersCommon *paramsc, struct Parameters2d *params2d, FILE *out, int file_number );





/* ���������� ����� � �������� ��� ���������� ����������� ����������� ����������� ������� 

   params - ��������� � ����������� ��������������� ������������ (in)
   restart_info_file - ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������ (in)

   **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out) */
void read_solution( struct Parameters2d *params2d, struct ParametersCommon* paramsc, int file_number, double ***u, double ***v, double ***w ) ;

#endif /* __IO_2D_H_ */
