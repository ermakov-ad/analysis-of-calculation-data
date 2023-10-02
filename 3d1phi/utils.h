// utils.h
// ������� �������, ����������� ��� ���������������� ������ ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// ��� ����������� �� ����������� ������
// (c) ����� �����, 2018
// ������: 26 ������� 2018 �.

#ifndef __UTILS_H_
#define __UTILS_H_

#include "struct.h"







// ���������� ���������� ������ � ������������� ����
// paramsc - ��������� � ��������� ����������� ��������������� ������������, ����� ����� �������� ��������� ���������� � ������������� ���� (in/out)
// params1d - ��������� � ����������� ���������� ������, ����� ����� �������� ��������� ���������� � ������������� ���� (in/out)
void dimensionalization( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d ) ;

//void current_block_number(struct ParametersCommon *params, struct Parameters1d *params1d, struct Parameters2d *params2d, int step_number_x, int step_number_y, int step_number_z, int *number_of_block);

#endif // __UTILS_H_