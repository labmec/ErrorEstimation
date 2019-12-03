import os

from neopz import *
from errorestimation import *

InitializePZLog()
Initialize2DUniformRefPatterns()

coord0 = TPZManVectorReal(3, 0.)
coord1 = TPZManVectorReal(3, 1.)
coord1[2] = 0.

bcIDs = TPZManVectorInt(4)
bcIDs[0] = -1
bcIDs[1] = -1
bcIDs[2] = -1
bcIDs[3] = -1

x_nel = 2
y_nel = 2

gmesh = Create2DGridMesh(x_nel, y_nel, coord0, coord1, bcIDs)

cfg = ProblemConfig()
cfg.NDivisions = 2

UniformRefinement(cfg.NDivisions, gmesh)

PrintGMeshToVTK(gmesh, "GMeshAfterRefinement.vtk")
PrintGMeshToTXT(gmesh, "GMeshAfterRefinementPython.txt")

cfg.gmesh = gmesh
cfg.Porder = 1
cfg.Hdivmais = 0
cfg.Dimension = 2
cfg.Prefine = False
cfg.Makepressurecontinuous = True
cfg.AdaptivityStep = 2
cfg.TensorNonConst = False

cfg.Materialids = {1}
cfg.Bcmaterialids = {-1}
cfg.Problemname = "EPython"
cfg.DirName = "TestePythonSinSin"
cfg.Exact.Exact = TLaplaceExample.ExactSol.ESinSin

multiphysicsCMesh = CreateHybridMultiphysicsMesh(cfg)

SolveHybridProblem(multiphysicsCMesh, cfg)

try:
    os.mkdir(cfg.DirName)
except OSError:
    print ("Creation of the directory %s failed. The directory already exists!" % cfg.DirName)
else:
    print ("Successfully created the directory %s" % cfg.DirName)

resultsFileH1 = cfg.DirName + "/H1ReconstructionErrors.vtk"
EstimateErrorWithH1Reconstruction(multiphysicsCMesh, cfg, resultsFileH1)

resultsFileHdiv = cfg.DirName + "/HdivReconstructionErrors.vtk"
EstimateErrorWithHdivReconstruction(multiphysicsCMesh, cfg, resultsFileHdiv)





