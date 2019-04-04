// Gmsh project created on Thu Apr  4 17:30:58 2019
SetFactory("OpenCASCADE");

DefineConstant[ arcdiv = {5, Min 3, Max 50, Step 1, Name "Parameters/Arc Division"} ];
DefineConstant[ crossdiv = {3, Min 3, Max 50, Step 1, Name "Parameters/Cross Division"} ];
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, -1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {-1, 0, 0, 1.0};
//+
Circle(1) = {4, 1, 5};
//+
Circle(2) = {4, 1, 2};
//+
Circle(3) = {5, 1, 3};
//+
Line(4) = {3, 1};
//+
Line(5) = {1, 2};
//+
Line(6) = {1, 4};
//+
Line(7) = {1, 5};
//+
Curve Loop(1) = {7, -1, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 2, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 3, 4};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {2, 1, 3} = arcdiv Using Progression 1;
//+
Transfinite Curve {5, 6, 7, 4} = crossdiv Using Progression 1;
//+
Physical Surface("domain") = {1, 2, 3};
//+
Physical Curve("dirichlet") = {1, 2, 3, 4, 5};
