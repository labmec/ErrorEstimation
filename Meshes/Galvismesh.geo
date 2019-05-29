// Gmsh project created on Mon May 27 10:09:47 2019
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
Point(6) = {1, 1, 0, 1.0};
//+
Point(7) = {-1, 1, 0, 1.0};
//+
Point(8) = {-1, -1, 0, 1.0};
//+
Point(9) = {1, -1, 0, 1.0};
//+
Line(1) = {2, 6};
//+
Line(2) = {6, 3};
//+
Line(3) = {3, 7};
//+
Line(4) = {7, 4};
//+
Line(5) = {4, 8};
//+
Line(6) = {8, 5};
//+
Line(7) = {5, 9};
//+
Line(8) = {9, 2};
//+
Line(9) = {2, 1};
//+
Line(10) = {1, 3};
//+
Line(11) = {1, 4};
//+
Line(12) = {1, 5};
//+
Curve Loop(1) = {1, 2, -10, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, -11, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 5, 6, -12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 12, 7, 8};
//+
Plane Surface(4) = {4};
//+
Physical Surface("omega1") = {3};
//+
Physical Surface("omega2") = {4};
//+
Physical Surface("omega3") = {2};
//+
Physical Surface("omega4") = {1};
//+
Physical Curve("boundary1") = {1};
//+
Physical Curve("boundary2") = {2};
//+
Physical Curve("boundary3") = {3};
//+
Physical Curve("boundary4") = {4};
//+
Physical Curve("boundary5") = {5};
//+
Physical Curve("boundary6") = {6};
//+
Physical Curve("boundary7") = {7};
//+
Physical Curve("boundary8") = {8};
