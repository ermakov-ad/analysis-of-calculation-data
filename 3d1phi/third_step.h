#include "ex5.h"
#include "struct.h"

//������� ������� �������� ���� ����� �����������

void calc_third_step(struct ParametersCommon paramsc, struct Parameters2d params2d, double ***u_next_3D, double ***v_next_3D, double ***w_next_3D, double ***u_first_step_3D, double ***v_first_step_3D, double ***w_first_step_3D,
                     double ***p_next_3D, double hx, double hy, double hz, double dt );