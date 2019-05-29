// Gmsh project created on Fri May 24 11:11:30 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {-1, 0, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {5, 1};
//+
Line(3) = {4, 1};
//+
Line(4) = {3, 1};
//+
Circle(5) = {3, 1, 2};
//+
Circle(6) = {2, 1, 5};
//+
Circle(7) = {5, 1, 4};
//+
Circle(8) = {4, 1, 3};
//+
Curve Loop(1) = {4, -1, -5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, -2, -6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, -3, -7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, -4, -8};
//+
Plane Surface(4) = {4};
//+
Physical Surface("omega1") = {1};
//+
Physical Surface("omega2") = {2};
//+
Physical Surface("omega3") = {3};
//+
Physical Surface("omega4") = {4};
//+
Physical Curve("boundary1") = {5};
//+
Physical Curve("boundary2") = {6};
//+
Physical Curve("boundary3") = {7};
//+
Physical Curve("boundary4") = {8};
