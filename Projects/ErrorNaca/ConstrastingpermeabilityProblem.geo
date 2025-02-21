// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};

//+
Point(5) = {-1, 0, 0, 1.0};
//+
Point(6) = {-1, 1, 0, 1.0};

//+
Point(7) = {-1, -1, 0, 1.0};
//+
Point(8) = {0, -1, 0, 1.0};

//+
Point(9) = {1, -1, 0, 1.0};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};

//+
Line(5) = {5, 1};
//+
Line(6) = {4, 6};
//+
Line(7) = {6, 5};

//+
Line(8) = {7, 8};
//+
Line(9) = {8, 1};
//+
Line(10) = {5, 7};

//+
Line(11) = {8, 9};
//+
Line(12) = {9, 2};


Curve Loop(1) = {1,2,3,4};
//+
Curve Loop(2) = {5,-4,6,7};
//+
Curve Loop(3) = {8,9,-5,10};
//+
Curve Loop(4) = {11,12,-1,-9};


Plane Surface(1) = {1};
//+
Plane Surface(2) = {2};
//+
Plane Surface(3) = {3};
//+
Plane Surface(4) = {4};


Physical Surface("Domain1", 1) = {1,3};
//+
Physical Surface("Domain2", 2) = {2,4};
//+
Physical Curve("Dirichlet1", 3) = {2,3,8,10};  
//+
Physical Curve("Dirichlet2", 4) = {6,7,11,12};
//+
Physical Curve("Neumman1", 5) = {};
//+
Physical Curve("Neumman2", 6) = {};
//+
Physical Point("Trailingedge", 7) = {1};

//+
Transfinite Curve {1,2,3,4,5,6,7,8,9,10,11,12} = 5 Using Progression 1;

Transfinite Surface {1} = {1,2,3,4};
//+
Transfinite Surface {2} = {5,1,4,6};
//+
Transfinite Surface {3} = {7,8,1,5};
//+
Transfinite Surface {4} = {8,9,2,1};
//+

Recombine Surface {1,2,3,4};