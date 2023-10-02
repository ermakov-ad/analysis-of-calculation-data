// memory.h
// ������� ��� ������ � �������
// (c) ����� �����, 2013
// ������: 12 ������� 2013 �.

#ifndef __MEMORY_H_
#define __MEMORY_H_

#include <stdio.h>
#include <stdlib.h>

// ��������� ������ ��� ���������� ������ ��������� ���� double
// size - ������ ������� (in) 
// **a - ��������� �� ������ (out)
void get_memory_for_1D_double_array( int size, double **a );

void get_memory_for_1D_int_array( int size, int **a ) ;

// ��������� ������ ��� ��������� ������ ��������� ���� double
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// ***a - ��������� �� ������ (out)
void get_memory_for_2D_double_array( int size1, int size2, double ***a );

// ��������� ������ ��� ��������� ������ ��������� ���� int
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// ***a - ��������� �� ������ (out)
void get_memory_for_2D_int_array( int size1, int size2, int ***a );

// ��������� ������ ��� ���������� ������ ��������� ���� double
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// size3 - ������ ����������� ������� (in)
// ****a - ��������� �� ������ (out)
void get_memory_for_3D_double_array( int size1, int size2, int size3, double ****a );

// ��������� ������ ��� ���������� ������ ��������� ���� double
// size1 - ������ ����������� ������� (in)
// size2 - ������ ����������� ������� (in)
// size3 - ������ ����������� ������� (in)
// ****a - ��������� �� ������ (out)
void get_memory_for_3D_int_array( int size1, int size2, int size3, int ****a );

#endif /* __MEMORY_H_ */