﻿// struct.h
// Общие структуры данных и константы для двухфазных сжимаемых солверов
// вне завимисоти от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 28 августа 2018 г.

#ifndef __STRUCT_H_
#define __STRUCT_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

//#define NX 101
//#define NZ 81
//#define dlx0 10
//#define dlz0 4

using namespace std;

// неизменные константы не для вынесения в конфигурационный файл


#define M1D 4 // размерность вектора-решения полной системы 1d задачи
#define M2D 4 // размерность вектора-решения полной системы 2d задачи

#define MAX_MATRIX_SIZE 1000000000 // максимальная размерность матрицы для решения уравнения Пуассона
#define SMALL_NUMBER 1.e-8  // маленькое число
#define	MAX_STRING_SIZE 1000 // максимальное число символов в строковых переменных
#define INFINITY1 1.e10 // большое число для процедуры выбора шага по времени
#define M_REDUCTION 3 // размерность вектора-решения редуцированной системы для случая отсутствия градиента объемной доли дисперсной фазы
#define K_GENERAL_CASE 4 // размерность вектора для описания условия на контактном разрыве в дисперсной фазе в общем случае наличия дисперсной фазы по обе стороны                          // от первоначального разрыва
#define K_SPECIAL_CASE 3 // размерность вектора для описания условия на контактном разрыве в дисперсной фазе в специальном случае отсутствия дисперсной фазы по одну из сторон от первоначального разрыва
#define MAX_IC_BLOCKS_NUMBER 1 // максимальное количество блоков, на которые можно разбить расчетную область для задания начальных условий
#define MAX_SENSORS_NUM 100 // максимальное количество датчиков для записи величин в определенных точках
#define ADIABATIC_CURVE_PNUM 101 // количество точек на адиабате - множестве состояний, куда можно перевести среду, пустив по ней ударную волну или волну разрежения, для отладки метода Годунова
#define ADAIABATIC_CURVE_P_B 1.e-5 // нижняя граница графика адиабаты на оси ординат (давлений), для отладки метода Годунова
#define ADAIABATIC_CURVE_P_T 0.1 // верхняя граница графика адиабаты на оси ординат (давлений), для отладки метода Годунова
#define RELAX_CASES_NUM 10 // количество веток алгоритма в процедуре расчета релаксации давлений
#define GODUNOV_CASES_NUM 10 // количество веток алгоритма в методе Годунова
#define BOUN_TYPE 0 // 1 - периодические ГУ, 0 - стенка с прилипанием

// новые типы для работы с многомерными матрицами с использованием контейнера vector

//typedef vector<double> array1D; // одномерный динамический вещественный массив
//typedef vector<array1D> array2D; // двумерный динамический вещественный массив

// новые типы из перечислений


enum Program { ONED2PHC = 1, ONED3PHC = 2, TWOD2PHC = 3 }; // запускаемая программа




/* Вектор "консервативных" переменных для полной системы уравнений
   (0) объемная доля дисперсной фазы ( b1 )
   (1) b1 * плотность дисперсной фазы ( b1 * r1 )
   (2) b1 * r1 * скорость дисперсной фазы ( b1 * r1 * v1 )
   (3) b1 * r1 * полная удельная энергия дисперсной фазы ( b1 * r1 * E1 )
   (4) b2 * r2
   (5) b2 * r2 * v2
   (6) b2 * r2 * E2 */

// тип граничного условия
enum Boundary_condition { WALL = 1, // стенка
                          FREE = 2, // свободное втекание/истечение
                          INFLOW = 3, // втекание с заданными параметрами
                          COMPLEX_INFLOW = 4, // втекание с одними заданными параметрами до заданного момента времени и с другими заданными после него
                          WALL_DISP = 5  }; // стенка для дисперсной фазы

enum Boundary_direction { LEFT_BOUNDARY = 1, // левая граница
                          RIGHT_BOUNDARY = 2, // правая граница
                          UP_BOUNDARY = 3, // верхняя граница
                          DOWN_BOUNDARY = 4 }; // нижняя граница

enum Block_types { NONE = 0, INNER = 1, OUTER = 2 }; // внутренний или внешний блок. 0 - нет ячеек



// структура с параметрами вычислительного эксперимента, не зависящих от размерности задачи
struct ParametersCommon {
    
    // поля, которые считываются из конфигурационного файла
    int program_name; // название запускаемой программы


    // шаг по времени
    bool constant_time_step; // является ли шаг по времени постоянным или выбирается динамически
    double cfl; // число Куранта-Фридрихса-Леви (имеет смысл, если шаг выбирается динамически)
    double dt; // шаг по времени (имеет смысл, если расчет ведется с постоянным шагом)


    int approximation_order; // порядок аппроксимации метода - 1 или 2

    // начальные условия
    bool isContinue; // является ли запускаемый расчет продолжанием старого?

    double mu; // коэффициент вязкости газа

    
    // режим расчета
    bool is_debug; // требуется ли выполнять проверку корректности матричных операций и решения систем уравнений

    // вывод результатов
    bool is_output_on_time; // является ли критерием вывода результатов достижение заданного момента времени
    double stop_time; // момент времени окончания расчета (если is_output_on_time == true)
    int stop_steps_number; // требуемое количество шагов по времени (если is_output_on_time == false)
    int output_number; // количество файлов с промежуточными результатами
    bool is_tecplot_header; // добавить ли заголовок для визуализации в Tecplot?

    // датчики
    bool are_sensors; // требуется ли использовать датчики в заданных точках?
    int sensors_num; // фактическое количество датчиков

    // "малые" константы в коде
    double eps_general; // регулярное малое число для сравнения вещественных чисел, использования в математических функциях,
                        // контроля сходимости итерационных процессов (если иное не оговорено особо)
    double eps_ludcmp; // малое число в функции ludcmp (math_utilc.cc), оставлено как было в оригинальной работе
    double eps_decouple; // максимально допустимый перепад по объемной доле слева и справа, при котором учитываются неконсервативные члены,
                         // иначе - расщепление конвективной части уравнений для фаз
    double eps_thin_layer; // точность решения системы нелинейных алгебраических уравнений "тонкого слоя" для проверки
    double eps_contact; // малое значение для задании скорости в случае стационарного контактного разрыва
    double eps_disp_abs; // малое значение объемной доли дисперсной фазы, начиная с которой считается, что она отсутствует
    double eps_cut_out; // малое значение объемной доли дисперсной фазы, начиная с которой при выводе результатов параметры
                        // дисперсной фазы кладуся равными фоновым параметрам
    double eps_relax_compaction; // малое число для обнуления переменных или выражений в релаксации давлений с компактированием,
                                 // равных машинному нулю

    // константы обезразмеривания
    bool use_dimensions; // задача будет определяться размерными параметрами?
    
    double mass_scale; // характерный масштаб массы
    double time_scale; // характерный масштаб времени
    double length_scale; // характерный масштаб длины
    double temperature_scale; // характерный масштаб температуры


    double omega;
    double f0;        
    // поля, которые заполняются в программе
    char input_file_directory[MAX_STRING_SIZE]; // директория, в которой находится конфигурационный файл задачи
    char output_file_directory[MAX_STRING_SIZE]; // директория, куда будут записаны результаты решения задачи

};

// структура с параметрами вычислительного эксперимента в одномерной постановке
struct Parameters1d {
    
    // поля, которые считываются из конфигурационного файла
       
    // информация о лагранжевых скалярах
    int number_of_scalars; // число скаляров 
    double scalar_values[MAX_IC_BLOCKS_NUMBER][5]; // массив характеристик скаляров

    // расчетная сетка (равномерная)
    double left_boundary_x; // координата левой границы расчетной области
    double right_boundary_x; // координата правой границы расчетной области
    int cells_number; // число ячеек в расчетной области

    int ic_blocks_number; // количество блоков, на которые разбита область для задания начальных условий
    int cell_begin[MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек начала блока (нумерация ячеек начинается с 0)
    int cell_end[MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек конца блока (нумерация ячеек начинается с 0)
    double block_values[MAX_IC_BLOCKS_NUMBER][M1D]; //  векторы примитивных переменных в блоках

    // граничные условия
    Boundary_condition left_bc; // граничное условие на левой границе
    Boundary_condition right_bc; // граничное условие на правой границе

    double inflow_values[M1D+5]; // параметры втекающего потока
    double change_inflow_time; // момент времени переключения на другие параметры втекания
    double changed_inflow_values[M1D+5]; // другие параметры втекающего потока после заданного момента времени

    // точное решение, только для БН
    bool build_exact_solution; // требуется ли строить точное решение
    int cells_number_for_exact_solution; // число ячеек для построения точного решения,
                                         // только для начальных условий из двух блоков
    
    double sensors_location[MAX_SENSORS_NUM]; // координаты датчиков
    

};

// структура с параметрами вычислительного эксперимента для 2d задачи
struct Parameters2d {    
    
    // расчетная сетка (равномерная)
    double left_boundary_x; // абсцисса левой границы расчетной области
    double right_boundary_x; // абсцисса правой границы расчетной области
    double down_boundary_y; // ордината нижней границы расчетной области
    double up_boundary_y; // ордината верхней границы расчетной области
    double front_boundary_z; // аппликата передней границы расчетной области
    double rear_boundary_z; // аппликата задней границы расчетной области
    int cells_number_x; // число ячеек по x
    int cells_number_y; // число ячеек по y
    int cells_number_z; // число ячеек по z
       
    int x_blocks_number; // количество блоков, на которые разбита область для задания начальных условий по x
    int y_blocks_number; // количество блоков, на которые разбита область для задания начальных условий по y
    int z_blocks_number; // количество блоков, на которые разбита область для задания начальных условий по z

    int x_begin[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек начала блока по x (нумерация ячеек начинается с 0) 
    int x_end[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек конца блока по x (нумерация ячеек начинается с 0)

    int y_begin[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек начала блока по y (нумерация ячеек начинается с 0)
    int y_end[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек конца блока по y (нумерация ячеек начинается с 0)

    int z_begin[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек начала блока по z (нумерация ячеек начинается с 0)
    int z_end[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек конца блока по z (нумерация ячеек начинается с 0)

    int block_type[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив типов блоков

    double block_values_p[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // векторы примитивных переменных в блоках
    double block_values_u[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER];
    double block_values_v[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER];
    double block_values_w[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER];
    // граничные условия
    Boundary_condition left_bc;  // граничное условие на левой границе
    Boundary_condition right_bc; // граничное условие на правой границе
    Boundary_condition down_bc;  // граничное условие на нижней границе
    Boundary_condition up_bc;    // граничное условие на верхней границе
    Boundary_condition front_bc;  // граничное условие на нижней границе
    Boundary_condition rear_bc;    // граничное условие на верхней границе

    double inflow_values_left_bc[M2D];  // параметры втекающего потока через левую ганицу
    double inflow_values_right_bc[M2D]; // параметры втекающего потока через правую границу
    double inflow_values_down_bc[M2D];  // параметры втекающего потока через нижнюю границу
    double inflow_values_up_bc[M2D];    // параметры втекающего потока через верхнюю границу
        
    // координаты датчиков
    double sensors_location_x[MAX_SENSORS_NUM];
    double sensors_location_y[MAX_SENSORS_NUM];
    
};

// для расщепления по пространственным направлениям
enum Direction2d { X_DIRECTION, Y_DIRECTION };

// структура, определяющий текущий момент расчета - либо момент времени, либо количество шагов
struct TimeMoment {
    int steps_num; // текущее количество шагов
    double curr_t; // текущее время
};

// структура с данными для отладки кода и исследования используемых численных методов
struct DebugInfo {

    int current_cell; // номер текущей ячейки
    double current_cell_x; // координата центра текущей ячейки
    double current_cell_vncons[M1D]; // вектор неконсервативных переменных в текущей ячейке
    
    // для отладки функций расчета потока
    int neighbour_cell; // номер соседней ячейки
    double neighbour_cell_vncons[M1D]; // вектор неконсервативных переменных в соседней ячейке

    // для отладки процедуры релаксации давлений
    FILE *relaxation_out; // дескриптор файла для записи информации о работе алгоритма релаксации давлений
    int relaxation_cases[RELAX_CASES_NUM]; // счетчики попаданий в различные ветки алгоритма при расчете релаксации давлений

    // для статистики частоты работы различных частей алгоритма метода Годунова
    FILE *godunov_out; // дескриптор файла для записи информации о работе алгоритма метода Годунова
    int godunov_cases[GODUNOV_CASES_NUM];
};

#endif // __STRUCT_H_