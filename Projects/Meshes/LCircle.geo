//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {-1, 0, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {4, 1, 3};//+
Circle(3) = {5, 1, 4};
//+
Line(4) = {1, 5};
//+
Line(5) = {1, 2};
//+
Curve Loop(1) = {5, 1, -2, -3, -4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("dirichlet") = {1, 2, 3, 4, 5};
//+
Physical Surface("domain") = {1};
