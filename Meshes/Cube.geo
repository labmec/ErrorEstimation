// Gmsh project created on Sun Jun 24 10:44:31 2018
SetFactory("OpenCASCADE");

nelem = 1;
//+
Block(1) = {0, 0, 0, 1, 1, 1};

Physical Volume("domain") = {1};
Physical Surface("dirichlet") = {1,2,3,4,5,6};

Transfinite Line {2,4,6,8,9,10,11,12} = 1+1*nelem Using Progression 1;
Transfinite Line {1,3,5,7} = 1+nelem Using Progression 1;

Transfinite Surface {1,2,3,4,5,6};

Recombine Surface {1,2,3,4,5,6};

Transfinite Volume {1};
