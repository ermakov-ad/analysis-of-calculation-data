// utils_2d.cc
// ������� �������, ����������� ��� ��������� ������� ������� ��������� ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// (c) ����� �����, 2013 - 2015
// ������: 16 ������� 2013 �.

#ifndef __UTILS_2D_H_
#define __UTILS_2D_H_

#include "ex5.h"
#include "struct.h"
#include "utils.h"




// ������������� �������-�������
// params2d - ��������� � ����������� ��������� ������ (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out)
void init_solution_3d( struct Parameters2d *params2d, struct ParametersCommon* paramsc, int file_num, struct TimeMoment *time_mom, double ***u, double ***v, double ***w );

// ������ ���� �������������� �� �������
// params - ��������� � ����������� ��������������� ������������
// *x - ���������� ����� �����
// **v_ncons - ������� ����������� ���������� � ������� �����
// ���������� ��� �������������� �� �������
double get_time_step_3d( const struct ParametersCommon* params, const struct Parameters2d* params2d, const double* x, const double* y, const double* z, double ***u, double ***v, double ***w );

void current_block_number_2d(struct ParametersCommon *params, struct Parameters2d *params2d, int step_number_x, int step_number_y, int *number_of_block);

#endif /* __UTILS_2D_H_ */