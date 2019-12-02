from neopz import *
from errorestimation import *

InitializePZLog()

coord0 = TPZManVectorReal(3, 0.)
coord1 = TPZManVectorReal(3, 1.)
coord1[2] = 0.

bcIDs = TPZManVectorInt(4)
bcIDs[0] = -1
bcIDs[1] = -1
bcIDs[2] = -1
bcIDs[3] = -1

x_nel = 2
y_nel = 3

gmesh = Create2DGridMesh(x_nel, y_nel, coord0, coord1, bcIDs)
PrintGMeshToVTK(gmesh, "GMeshBeforeRefinement.vtk")

#UniformRefinement(2, gmesh)

cfg = ProblemConfig()
cfg.gmesh = gmesh
cfg.Porder = 1
cfg.Dimension = 2
cfg.Prefine = 0
cfg.Makepressurecontinuous = 1
cfg.TensorNonConst = 0
cfg.Hdivmais = 1
cfg.Materialids = {1}
cfg.Bcmaterialids = {-1}
cfg.Problemname = "EPython"

multiphysicsCMesh = CreateHDivMesh(cfg)
multiphysicsCMesh.InitializeBlock()

HybridizeCompMesh(multiphysicsCMesh)
SolveHybridProblem(multiphysicsCMesh, cfg)
#EstimateErrorsWithH1Reconstruction(multiphysicsCMesh, cfg)

