// io.h
// ����� ������� �����/������, �� ��������� �� ����������� ������
// (c) ����� �����, 2018
// ������: 26 ������� 2018 �.

#ifndef __IO_H_
#define __IO_H_

#include "struct.h"

// ���������� ��������� � ��������� ����������� ��������������� ������������
// params - ���������� ����������������� ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
void fill_parameters_common( FILE *params, struct ParametersCommon* paramsc );

// ������ ����� � ����������� ��� ����������� ����������� �������
// files_directory - ����������, � ������� ���������� ������ (in)
// time_moment - ������ � ����������� � ������� ������� ��� ����� �����, ��� �������
// � ��������� ��� ��� ������� ���� � �������������� ������������ (in)
// filename - ��� ���������� ����������� ����� (in)
void write_restart_info( char *files_directory, char *time_moment, char *filename );



/* ���������� ����� � �������� ��� ���������� ����������� ����������� ����������� ������� 

   params - ��������� � ����������� ��������������� ������������ (in)
   restart_info_file - ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������ (in)

   **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out) */
void read_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, FILE *restart_info_file, double **initial_solution );

// ������� ������� ���� � ����������� � �������� �, � ������ ������, ��������� �����������
// ����� � ���������� ����������� ���������������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in) 
// files_directory - ����������, � ������� ���������� ������ (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// file_to_read - ���������� ����� � ��������������� ��� �������� (out)
void try_to_restart( struct ParametersCommon *paramsc, char *files_directory, struct TimeMoment *time_mom, FILE *file_to_read );

#endif /* __IO_COMMON_H_ */