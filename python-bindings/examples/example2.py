from neopz import *
from errorestimation import *

coord0 = TPZManVectorReal(3, 0.);
coord1 = TPZManVectorReal(3, 1.);
coord1[2] = 0.;

bcIDs = TPZManVectorInt(4);
bcIDs[0] = 1;
bcIDs[1] = 2;
bcIDs[2] = 3;
bcIDs[3] = 4;

x_nel = 2;
y_nel = 3;

gmesh = Create2DGridMesh(x_nel, y_nel, coord0, coord1, bcIDs)
PrintGMeshToVTK(gmesh, "GMeshBeforeRefinement.vtk")

UniformRefinement(2, gmesh)
PrintGMeshToVTK(gmesh, "GMeshAfterRefinement.vtk")


