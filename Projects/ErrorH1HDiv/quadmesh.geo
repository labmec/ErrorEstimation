elem = 2;
//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Surface("dom") = {1};
//+
Physical Line("bcD") = {1};
Physical Line("bcR") = {2};
Physical Line("bcT") = {3};
Physical Line("bcL") = {4};
//+
Transfinite Surface {1} = {1, 2, 3, 4} Alternated;

Transfinite Line {1, 3} = elem+1 Using Progression 1;
Transfinite Line {2, 4} = elem+1 Using Progression 1;
//+
Recombine Surface {1};
//+
Show "*";
