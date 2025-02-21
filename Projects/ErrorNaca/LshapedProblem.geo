// Gmsh project created on Mon Jan 20 10:30:05 2025
SetFactory("OpenCASCADE");

// Define points
Point(1) = {-1, -1, 0, 1.0};
Point(2) = {0, -1, 0, 1.0};
Point(3) = {0, 0, 0, 1.0};
Point(4) = {-1, 0, 0, 1.0};

Point(5) = {1, 0, 0, 1.0};
Point(6) = {1, 1, 0, 1.0};
Point(7) = {-1, 1, 0, 1.0};

// Define lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {3, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 4};

// Define curve loops and surfaces
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {-3, 5, 6, 7, 8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Physical groups
Physical Surface("Domain", 1) = {1, 2};
Physical Curve("Dirichlet", 2) = {1, 2, 5, 6, 7, 8, 4};
Physical Curve("Neumann", 3) = {};
Physical Point("Trailingedge", 4) = {3};
Physical Point("Fixedpoint", 5) = {1};

// Transfinite definitions
Transfinite Curve {1, 2, 3, 4} = 5 Using Progression 1;
Transfinite Curve {5, 6, 8} = 5 Using Progression 1;
Transfinite Curve {7} = 9 Using Progression 1;

Transfinite Surface {1} = {1, 2, 3, 4};
Transfinite Surface {2} = {4, 5, 6, 7};

// Recombine surfaces
Recombine Surface {1, 2};