elem = 10;
//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Surface("domain") = {1};
//+
Physical Line("dirichlet") = {1, 2, 3, 4};
//+
Transfinite Surface {1} = {1, 2, 3, 4} Alternated;

Transfinite Line {1, 3} = elem+1 Using Progression 1;
Transfinite Line {2, 4} = elem+1 Using Progression 1;
//+
Recombine Surface {1};
