// Gmsh project created on Thu May  2 10:05:54 2024
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {-1, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, -1, 0, 1.0};
//+
Point(5) = {-1, -1, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};
//+
Point(7) = {0, -1, 0, 1.0};
//+
Point(8) = {-1, 0, 0, 1.0};
//+
Point(9) = {1, 0.00001, 0, 1.0};
Point(10) = {1, -0.00001, 0, 1.0};
//+
Line(1) = {3, 9};
//+
Line(2) = {10, 4};
//+
Line(3) = {4, 7};
//+
Line(4) = {7, 5};
//+
Line(5) = {5, 8};
//+
Line(6) = {8, 1};
//+
Line(7) = {1, 6};
//+
Line(8) = {6, 2};
//+
Line(9) = {1, 7};
//+
Line(10) = {3, 6};
//+
Line(11) = {9, 1};
//+
Line(12) = {10, 1};
//+
Line(13) = {2, 8};
//+
Curve Loop(1) = {1, 11, 7, -10};
//+
Plane Surface(1) = {1};

//+
Curve Loop(2) = {8, 13, 6, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, 9, 4, 5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 9, -3, -2};
//+
Plane Surface(4) = {4};
//+
Physical Curve("Boundaries", 14) = {10, 8, 4, 3, 1, 2, 13, 5, 11, 12};
//+
Physical Surface("Domain", 19) = {1, 2, 4, 3};
//+
Transfinite Curve {1, 10, 7, 8, 6, 11, 12, 2, 3, 9, 4, 5, 13} = 11 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3} Right;
//+
Transfinite Surface {4};

Recombine Surface {1,2,3,4};
