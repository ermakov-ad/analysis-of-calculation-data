// параметр, задающий количество ячеек в пограничном слое
count_of_boundary_cells = 10;
// параметр, задающий количество ячеек основном объёме
count_of_inner_cells = 100;
// параметр, задающий толщину пограничного слоя
thickness_of_boundary_layer = 0.0632;
// размер рассчетной области
size_of_calculation_volume = 6.2832;
half_of_size = size_of_calculation_volume / 2.0;

//+
Point(1) = {-half_of_size + thickness_of_boundary_layer, -half_of_size + thickness_of_boundary_layer, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Point(2) = {half_of_size - thickness_of_boundary_layer, -half_of_size + thickness_of_boundary_layer, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Point(3) = {half_of_size - thickness_of_boundary_layer, half_of_size - thickness_of_boundary_layer, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Point(4) = {-half_of_size + thickness_of_boundary_layer, half_of_size - thickness_of_boundary_layer, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Point(5) = {-half_of_size, -half_of_size, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Point(6) = {half_of_size, -half_of_size, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Point(7) = {half_of_size, half_of_size, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Point(8) = {-half_of_size, half_of_size, -half_of_size + thickness_of_boundary_layer, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {1, 5};
//+
Line(11) = {2, 6};
//+
Line(12) = {3, 7};
//+
Line(13) = {4, 8};
//+
Curve Loop(1) = {13, -7, -12, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, -6, -11, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, -5, -9, 1};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, -8, -13, 4};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 1, 2, 3};
//+
Plane Surface(5) = {5};
//+
Transfinite Surface {5};
//+
Transfinite Surface {1};
//+
Transfinite Surface {4};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
// количество тетраэдров вдоль каждой оси в объёме слое
Transfinite Curve {1, 2, 3, 4, 8, 7, 6, 5} = count_of_inner_cells Using Progression 1;
//+
// количество тетраэдров в погран слое
Transfinite Curve {11, 9, 12, 13} = count_of_boundary_cells Using Progression 1;
//+
Extrude {0, 0, size_of_calculation_volume - 2.0 * thickness_of_boundary_layer} {
  Surface{1}; Curve{8}; Curve{5}; Curve{6}; Curve{7}; Curve{3}; Curve{4}; Curve{1}; Curve{2}; Curve{12}; Surface{4}; Curve{13}; Curve{9}; Curve{11}; Point{3}; Point{7}; Point{8}; Point{4}; Point{1}; Point{5}; Surface{3}; Surface{5}; Point{2}; Point{6}; Surface{2}; Layers {count_of_inner_cells}; Recombine;
}
//+
Extrude {0, 0, thickness_of_boundary_layer} {
  Curve{44}; Curve{40}; Curve{36}; Curve{16}; Curve{56}; Curve{52}; Curve{48}; Curve{18}; Curve{17}; Curve{82}; Curve{61}; Curve{15}; Point{14}; Point{18}; Point{22}; Point{26}; Point{20}; Point{24}; Point{9}; Point{10}; Surface{151}; Surface{35}; Surface{107}; Surface{81}; Surface{129}; Layers {count_of_boundary_cells}; Recombine;
}
//+
Extrude {0, 0, -thickness_of_boundary_layer} {
  Curve{5}; Curve{6}; Curve{7}; Curve{8}; Curve{1}; Curve{2}; Curve{3}; Curve{4}; Curve{9}; Curve{11}; Curve{12}; Curve{13}; Point{5}; Point{1}; Point{6}; Point{2}; Point{7}; Point{3}; Point{4}; Point{8}; Surface{3}; Surface{4}; Surface{2}; Surface{1}; Surface{5}; Layers {count_of_boundary_cells}; Recombine;
}
//+
Physical Surface("top", 468) = {309, 287, 265, 221, 243};
//+
Physical Surface("bottom", 469) = {467, 379, 401, 423, 445};
//+
Physical Surface("wall", 470) = {39, 26, 47, 43, 313, 317, 321, 325, 159, 155, 167, 163};
//+
Physical Volume("inner_volume", 471) = {4};
//+
Physical Volume("top_boundary_layer", 472) = {10, 7, 9, 8, 6};
//+
Physical Volume("bottom_boundary_layer", 473) = {15, 11, 13, 14, 12};
//+
Physical Volume("wall_boundary_layer", 474) = {2, 1, 5, 3};
