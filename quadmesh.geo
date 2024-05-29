//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Surface("dom", 5) = {1};
//+
Physical Curve("bcD", 6) = {1};
//+
Physical Curve("bcR", 7) = {2};
//+
Physical Curve("bcT", 8) = {3};
//+
Physical Curve("bcL", 9) = {4};
//+
Transfinite Curve {1, 2, 3, 4} = 5 Using Progression 1;
//+
SetFactory("OpenCASCADE");
//+
Transfinite Surface {1};
