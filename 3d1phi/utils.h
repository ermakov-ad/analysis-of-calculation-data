// utils.h
// Рабочие функции, специфичные для рассматриваемого класса систем уравнений типа Baer-Nunziato и Saurel-Abgrall
// вне зависимости от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 26 августа 2018 г.

#ifndef __UTILS_H_
#define __UTILS_H_

#include "struct.h"







// Приведение параметров задачи к безразмерному виду
// paramsc - структура с основными параметрами вычислительного эксперимента, часть полей проходит процедуру приведения к безразмерному виду (in/out)
// params1d - структура с параметрами одномерной задачи, часть полей проходит процедуру приведения к безразмерному виду (in/out)
void dimensionalization( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d ) ;

//void current_block_number(struct ParametersCommon *params, struct Parameters1d *params1d, struct Parameters2d *params2d, int step_number_x, int step_number_y, int step_number_z, int *number_of_block);

#endif // __UTILS_H_