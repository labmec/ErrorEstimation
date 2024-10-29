// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {-1, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Point(5) = {-1, 1, 0, 1.0};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};
//+
Curve Loop(1) = {1,2,3,4,5};
//+
Plane Surface(1) = {1};
//+

Physical Surface("Domain", 1) = {1};
//+
Physical Curve("Dirichlet", 2) = {1,3,4,5};
//+
Physical Curve("Neumman", 3) = {2};
//+
Physical Point("Trailingedge", 4) = {2};
//+

Transfinite Surface {1} = {1,3,4,5};
//+
Transfinite Curve {1,2} = 5 Using Progression 1;
//+
Transfinite Curve {3,5} = 5 Using Progression 1;
//+
Transfinite Curve {4} = 9 Using Progression 1;
//+
Recombine Surface {1};
