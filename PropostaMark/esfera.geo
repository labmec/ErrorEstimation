
//+
rad = DefineNumber[ 10, Name "Parameters/rad" ];
//+
ndens = DefineNumber[ 1, Name "Parameters/ndens" ];
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {rad, 0, 0, 1.0};
//+
Point(3) = {-rad, 0, 0, 1.0};
//+
Point(4) = {0, rad, 0, 1.0};
//+
Point(5) = {0, -rad, 0, 1.0};
//+
Point(6) = {0, 0, rad, 1.0};
//+
Point(7) = {0, 0, -rad, 1.0};
//+
Circle(1) = {2, 1, 4};
//+
Circle(2) = {4, 1, 3};
//+
Circle(3) = {3, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Circle(5) = {2, 1, 6};
//+
Circle(6) = {6, 1, 3};
//+
Circle(7) = {3, 1, 7};
//+
Circle(8) = {7, 1, 2};
//+
Circle(9) = {4, 1, 6};
//+
Circle(10) = {6, 1, 5};
//+
Circle(11) = {5, 1, 7};
//+
Circle(12) = {7, 1, 4};
//+
Curve Loop(1) = {2, 7, 12};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {2, -6, -9};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {3, -10, 6};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {3, 11, -7};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {4, -8, -11};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {4, 5, 10};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {1, 9, -5};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {1, -12, 8};
//+
Surface(8) = {8};
//+
Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
//+
Volume(1) = {1};
//+
Physical Volume("domain") = {1};
//+
Physical Surface("boundary") = {8, 7, 6, 5, 2, 1, 4, 3};
