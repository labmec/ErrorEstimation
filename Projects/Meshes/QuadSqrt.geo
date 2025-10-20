elem = 4; // number of divisions 
//-----------------------------------------------------------
// Two adjacent squares forming a 2x1 rectangle mesh
//-----------------------------------------------------------

SetFactory("OpenCASCADE");

// --- Define points ---
Point(1) = {-1, 0, 0};
Point(2) = {0, 0, 0};
Point(3) = {1, 0, 0};
Point(4) = {-1, 1, 0};
Point(5) = {0, 1, 0};
Point(6) = {1, 1, 0};
Point(7) = {-1, -1, 0};
Point(8) = {0, -1, 0};
Point(9) = {1, -1, 0};
Point(10) = {1,0,0};

// --- Define lines for left square ---
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {4, 1};

// --- Define lines for right square ---
Line(5) = {2, 3};
Line(6) = {3, 6};
Line(7) = {6, 5};



Line(8) = {7,8};
Line(9) = {8,2};
Line(10) = {1,7};
Line(11) = {8,9};
Line(12) = {9,10};
Line(13) = {10,2};

// --- Define line loops and surfaces ---
Line Loop(1) = {1, 2, 3, 4};  // left square above
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, -2}; // right square above
Plane Surface(2) = {2};

Line Loop(3) = {8, 9, -1, 10};  // left square below
Plane Surface(3) = {3};

Line Loop(4) = {11, 12, 13, -9}; // right square below
Plane Surface(4) = {4};

// --- Mesh control (transfinite for rectangular mesh) ---
Transfinite Line {1, 3, 5, 7, 8, 11, 13} = elem+1 Using Progression 1; // horizontal divisions
Transfinite Line {2, 4, 6, 10, 9, 12} = elem+1 Using Progression 1;  // vertical divisions

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Recombine Surface {1, 2, 3, 4}; // make quads instead of triangles

// --- Physical groups (optional) ---
Physical Surface("upper") = {1,2};
Physical Surface("lower") = {3,4};
Physical Line("neumannupper") = {5}; // common edge between squares
Physical Line("neumannlower") = {13}; // common edge between squares
Physical Line("contourupper") = {6,7,3,4};
Physical Line("contourlower") = {10,8,11,12};
Physical Point("point_singularity") = {2};

// --- Generate 2D mesh ---
Mesh 2;//+

