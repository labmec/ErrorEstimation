Point(1) = {-1, -1, 0, 2.0};
//+
Point(2) = {1, -1, 0, 2.0};
//+
Point(3) = {1, 1, 0, 2.0};
//+
Point(4) = {-1, 1, 0, 2.0};
//+
Point(5) = {0, 0, 0, 2.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {1, 5};
//+
Line(6) = {2, 5};
//+
Line(7) = {3, 5};
//+
Line(8) = {4, 5};
//+
Curve Loop(1) = {1, 6, -5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 7, -6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 8, -7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 5, -8};
//+
Plane Surface(4) = {4};
//+
Physical Surface("Omega1") = {2};
//+
Physical Surface("Omega2") = {3};
//+
Physical Surface("Omega3") = {4};
//+
Physical Surface("Omega4") = {1};
//+
Physical Curve("neumann1") = {2};
//+
Physical Curve("neumann2") = {3};
//+
Physical Curve("neumann3") = {4};
//+
Physical Curve("neumann4") = {1};
