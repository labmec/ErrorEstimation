# #%% Importing the libraries
from TPZMeshModeling import TPZMeshModeling
from numpy import cos, pi 
import gmsh

def main():
    #%% Initial gmsh settings 
    TPZMeshModeling.Begin() # initializing gmsh

    TPZMeshModeling.TurnOnLabels('curve', 'surface') # turning on the labels for points and curves

    #%% -------- INPUT DATA --------
    # Any changes you wann do, do it here

    # file name
    file_name: str = "naca" 

    # z-coordinate
    depth: float = 0. 

    # mesh size (later the mesh size will be set, but gmsh requires a value to create the points and lines)
    lc: float = 1e-1 

    # number of elements / edge (that's how we set the mesh size)
    nEl: int = 1 
    # number of nodes / edge
    nNodes: int = nEl + 1 
    # mesh dimension
    mesh_dim: int = 2 
    # -------- END OF INPUT DATA --------

    #%% Creating the naca profile
    point_coord: list[float] = [ 
    [ 0, 0 , depth], # p1 
    [ 0.40000000000000036, 0.322772, depth], # p2, p10
    [ 1.6000000000000005, 0.544174, depth], # p3, p9
    [ 3.5999999999999996, 0.592421, depth], # p4, p8
    [ 6.4000000000000012, 0.420066, depth], # p5, p7
    [ 10, 0, depth] # p6
    ]

    # the botton part of the profile
    point_coord += [ [x[0], -x[1], depth] for x in point_coord[::-1] if x[1] != 0 ]

    # Creating points
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = TPZMeshModeling.CreatePoints(point_coord, lc)


    # Creating lines on the surface of the profile
    lines: list[int] = [
        [p1, p2], # l1
        [p2, p3], # l2
        [p3, p4], # l3
        [p4, p5], # l4
        [p5, p6], # l5
        [p6, p7], # l6
        [p7, p8], # l7
        [p8, p9], # l8
        [p9, p10], # l9
        [p10, p1] # l10
    ]

    l1, l2, l3, l4, l5, l6, l7, l8, l9, l10 = TPZMeshModeling.CreateLines(lines)

    #%% Creating the boundary (these lines you use to apply BCs)
    boundary_coord: list[float] = [
        [-2, 0, depth], # p11
        [cos(pi/4) * -2, cos(pi/4) * -2, depth], # p12
        [point_coord[2][0], -2, depth], # p13
        [point_coord[3][0], -2, depth], # p14
        [point_coord[4][0], -2, depth], # p15
        [point_coord[5][0], -2, depth], # p16
        [14, -2, depth], # p17
        [14, 0, depth], # p18
    ]

    boundary_coord += [ [x[0], -x[1], depth] for x in boundary_coord[::-1] if x[1] != 0 ]

    # Create points
    p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24 = TPZMeshModeling.CreatePoints(boundary_coord, lc)

    lines_boundary: list[int] = [
        [p11, p12], # l11
        [p12, p13], # l12
        [p13, p14], # l13
        [p14, p15], # l14
        [p15, p16], # l15
        [p16, p17], # l16
        [p17, p18], # l17
        [p18, p19], # l18
        [p19, p20], # l19
        [p20, p21], # l20
        [p21, p22], # l21
        [p22, p23], # l22
        [p23, p24], # l23
        [p24, p11] # l24
    ]

    l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24 = TPZMeshModeling.CreateLines(lines_boundary)

    #%% Creating inner lines, from the profile to the boundary
    inner_lines: list[int] = [
        [p1, p11], # l25  Discontinuty - first try
        [p10, p12], # l26
        [p9, p13], # l27
        [p8, p14], # l28
        [p7, p15], # l29
        [p6, p16], # l30
        [p6, p18], # l31 Discontinuty - second try
        [p6, p20], # l32
        [p5, p21], # l33
        [p4, p22], # l34
        [p3, p23], # l35
        [p2, p24] # l36
    ]

    l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = TPZMeshModeling.CreateLines(inner_lines)

    #%% Creating the surfaces
    curve_loops: list[int] = [
        [l11, l26, l10, l25], # cl1
        [l12, l27, l9, l26], # cl2
        [l13, l28, l8, l27], # cl3
        [l14, l29, l7, l28], # cl4
        [l15, l30, l6, l29], # cl5
        [l16, l17, l31, l30], # cl6
        [l31, l18, l19, l32], # cl7
        [l32, l20, l33, l5], # cl8
        [l33, l21, l34, l4], # cl9
        [l34, l22, l35, l3], # cl10
        [l35, l23, l36, l2], # cl11
        [l36, l24, l25, l1] # cl12
    ]

    # curve loops
    curves: list[int] = TPZMeshModeling.CreateCurveLoops(curve_loops)
    curves_: list[list[int]] = [[c] for c in curves]

    # and surfaces
    surfaces: list[int] = TPZMeshModeling.CreatePlanes(curves_)

    # synchronize the model with gmsh
    TPZMeshModeling.Synchronize() 

    #%% -------- PHYSICAL GROUPS --------
    # No physical group has been created yet for the boundary conditions
    # to do so, use the following command:
    dim = 0 # dimension of the entity
    tag = [p1] # tag of the entity
    ID = 4 
    gp = [
#        [(dim, tag), ID, "name"]
        [(dim, tag), ID, "FixPoint"]
    ]
    TPZMeshModeling.CreatePhysicalGroup(gp)
    dim = 0 # dimension of the entity
    tag = [p6] # tag of the entity
    ID = 6 
    gp = [
#        [(dim, tag), ID, "name"]
        [(dim, tag), ID, "Trailingedge"]
    ]
    TPZMeshModeling.CreatePhysicalGroup(gp)

    dim = 1 # dimension of the entity
    tag = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10] # tag of the entity
    ID = 2
    gp = [
#        [(dim, tag), ID, "name"]
        [(dim, tag), ID, "Profile"]
    ]
    TPZMeshModeling.CreatePhysicalGroup(gp)

    dim = 1 # dimension of the entity
    tag = [l31] # tag of the entity
    ID = 5
    gp = [
#        [(dim, tag), ID, "name"]
        [(dim, tag), ID, "Cut"]
    ]
    TPZMeshModeling.CreatePhysicalGroup(gp)
 
    dim = 1 # dimension of the entity
    tag = [l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24] # tag of the entity
    ID = 3
    gp = [
#        [(dim, tag), ID, "name"]
        [(dim, tag), ID, "OuterBoundary"]
    ]
    TPZMeshModeling.CreatePhysicalGroup(gp)

    # dim: dimension of the entity
    # tag: tag of the entity
    # ID: ID of the physical group
    # name: name of the physical group

    # Example for creation of domain

    gp: list[tuple[int, list[int]], int, str] = [
        [(2, surfaces), 1, "Domain"],
        # keep adding the physical groups...
    ]

    TPZMeshModeling.CreatePhysicalGroup(gp)

    #%% Setting the mesh

    # domain lines
    domain_lines: list[int] = [i + 1 for i in range(36)] 

    # creating the transfinte lines and surfaces
    TPZMeshModeling.TransfiniteCurve(domain_lines, nNodes)
    TPZMeshModeling.TransfiniteSurface(surfaces, "Left") # it really doesn't matter whether you use "Left" or "Right"

    # recombining elements to have quadrilateral ones
    TPZMeshModeling.RecombineSurface(surfaces)

    # generating mesh
    TPZMeshModeling.CreateMesh(mesh_dim) 

    # showing model on gmsh
    # TPZMeshModeling.ShowModel()

    # writing .msh file
    TPZMeshModeling.WriteMeshFiles(file_name, ".msh")

    path = "/home/cordeiro/projects/ErrorEstimation/build/Projects/ErrorNaca/"
    TPZMeshModeling.MoveFiles(file_name, None, path, ".msh")

    # closing gmsh
    TPZMeshModeling.End()

#%% Defining main function
if __name__ == "__main__":
    main()




# #%% Importing the libraries
# from TPZMeshModeling import TPZMeshModeling
# from numpy import cos, pi 
# import gmsh

# def main():
#     #%% Initial gmsh settings 
#     TPZMeshModeling.Begin() # initializing gmsh

#     TPZMeshModeling.TurnOnLabels('curve', 'surface') # turning on the labels for points and curves

#     #%% -------- INPUT DATA --------
#     # Any changes you wann do, do it here

#     # file name
#     file_name: str = "naca" 

#     # z-coordinate
#     depth: float = 0. 

#     # mesh size (later the mesh size will be set, but gmsh requires a value to create the points and lines)
#     lc: float = 1e-1 

#     # number of elements / edge (that's how we set the mesh size)
#     nEl: int = 1 
#     # number of nodes / edge
#     nNodes: int = nEl + 1 
#     # mesh dimension
#     mesh_dim: int = 2 
#     # -------- END OF INPUT DATA --------

#     #%% Creating the naca profile
#     point_coord: list[float] = [ 
#     [ 0, 0 , depth], # p1 
#     [ 0.40000000000000036, 0.322772, depth], # p2, p11
#     [ 1.6000000000000005, 0.544174, depth], # p3, p10
#     [ 3.5999999999999996, 0.592421, depth], # p4, p9
#     [ 6.4000000000000012, 0.420066, depth], # p5, p8
#     [ 10, 0, depth], # p6
#     [ 10, 0, depth], # p7
#     ]

#     # the botton part of the profile
#     point_coord += [ [x[0], -x[1], depth] for x in point_coord[::-1] if x[1] != 0 ]

#     # Creating points
#     p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 = TPZMeshModeling.CreatePoints(point_coord, lc)


#     # Creating lines on the surface of the profile
#     lines: list[int] = [
#         [p1, p2], # l1
#         [p2, p3], # l2
#         [p3, p4], # l3
#         [p4, p5], # l4
#         [p5, p6], # l5
#         [p7, p8], # l6
#         [p8, p9], # l7
#         [p9, p10], # l8
#         [p10, p11], # l9
#         [p11, p1] # l10
#     ]

#     l1, l2, l3, l4, l5, l6, l7, l8, l9, l10 = TPZMeshModeling.CreateLines(lines)

#     #%% Creating the boundary (these lines you use to apply BCs)
#     boundary_coord: list[float] = [
#         [-2, 0, depth], # p12
#         [cos(pi/4) * -2, cos(pi/4) * -2, depth], # p13, p26
#         [point_coord[2][0], -2, depth], # p14, p25
#         [point_coord[3][0], -2, depth], # p15, p24
#         [point_coord[4][0], -2, depth], # p16, p23
#         [point_coord[5][0], -2, depth], # p17, p22
#         [14, -2, depth], # p18, p21
#         [14, 0, depth], # p19
#         [14, 0, depth], # p20
#     ]

#     boundary_coord += [ [x[0], -x[1], depth] for x in boundary_coord[::-1] if x[1] != 0 ]

#     # Create points
#     p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26 = TPZMeshModeling.CreatePoints(boundary_coord, lc)

#     lines_boundary: list[int] = [
#         [p12, p13], # l11
#         [p13, p14], # l12
#         [p14, p15], # l13
#         [p15, p16], # l14
#         [p16, p17], # l15
#         [p17, p18], # l16
#         [p18, p20], # l17
#         [p19, p21], # l18
#         [p21, p22], # l19
#         [p22, p23], # l20
#         [p23, p24], # l21
#         [p24, p25], # l22
#         [p25, p26], # l23
#         [p26, p12], # l24
#     ]

#     l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24 = TPZMeshModeling.CreateLines(lines_boundary)

#     #%% Creating inner lines, from the profile to the boundary
#     inner_lines: list[int] = [
#         [p1, p12], # 125 
#         [p11, p13], # 126 
#         [p10, p14], # l27 
#         [p9, p15], # l28 
#         [p8, p16], # l29 
#         [p17, p7], # l30 
#         [p6, p22], # l31 (old l32)
#         [p23, p5], # l32 (old l33)
#         [p4, p24], # l33 (old l34)
#         [p3, p25], # l34 (old l35)
#         [p2, p26] # l35 (old l36)
#     ]

#     l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35 = TPZMeshModeling.CreateLines(inner_lines)

#     #%% Creating the cut lines: cut_line + and cut_line -
#     cut_lines: list[int] = [
#         [p6, p19], # l36  Cut line + (old 131)
#         [p20, p7], # l37 Cut line - (new line)
#     ]

#     l36, l37 = TPZMeshModeling.CreateLines(cut_lines)

#     #%% Creating the surfaces
#     curve_loops: list[int] = [
#         [l11, l26, l10, l25], # cl1
#         [l12, l27, l9, l26], # cl2
#         [l13, l28, l8, l27], # cl3
#         [l14, l29, l7, l28], # cl4
#         [l15, l30, l6, l29], # cl5
#         [l16, l17, l37, l30], # cl6
#         [l36, l18, l19, l31], # cl7
#         [l31, l20, l32, l5], # cl8
#         [l32, l21, l33, l4], # cl9
#         [l33, l22, l34, l3], # cl10
#         [l34, l23, l35, l2], # cl11
#         [l35, l24, l25, l1] # cl12
#     ]

#     # curve loops
#     curves: list[int] = TPZMeshModeling.CreateCurveLoops(curve_loops)
#     curves_: list[list[int]] = [[c] for c in curves]

#     # and surfaces
#     surfaces: list[int] = TPZMeshModeling.CreatePlanes(curves_)

#     # synchronize the model with gmsh
#     TPZMeshModeling.Synchronize() 

#     #%% -------- PHYSICAL GROUPS --------
#     # No physical group has been created yet for the boundary conditions
#     # to do so, use the following command:

#     dim = 0 # dimension of the entity
#     tag = [p1] # tag of the entity
#     ID = 1 
#     gp = [
# #        [(dim, tag), ID, "name"]
#         [(dim, tag), ID, "FixPoint"]
#     ]
#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     dim = 0 # dimension of the entity
#     tag = [p6] # tag of the entity
#     ID = 2 
#     gp = [
# #        [(dim, tag), ID, "name"]
#         [(dim, tag), ID, "Trailingedge_plus"]
#     ]
#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     dim = 0 # dimension of the entity
#     tag = [p7] # tag of the entity
#     ID = 3 
#     gp = [
# #        [(dim, tag), ID, "name"]
#         [(dim, tag), ID, "Trailingedge_minus"]
#     ]
#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     dim = 1 # dimension of the entity
#     tag = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10] # tag of the entity
#     ID = 4
#     gp = [
# #        [(dim, tag), ID, "name"]
#         [(dim, tag), ID, "Profile"]
#     ]
#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     dim = 1 # dimension of the entity
#     tag = [l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24] # tag of the entity
#     ID = 5
#     gp = [
# #        [(dim, tag), ID, "name"]
#         [(dim, tag), ID, "OuterBoundary"]
#     ]
#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     dim = 1 # dimension of the entity
#     tag = [l36] # tag of the entity
#     ID = 6
#     gp = [
# #        [(dim, tag), ID, "name"]
#         [(dim, tag), ID, "Cutplus"]
#     ]
#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     dim = 1 # dimension of the entity
#     tag = [l37] # tag of the entity
#     ID = 7
#     gp = [
# #        [(dim, tag), ID, "name"]
#         [(dim, tag), ID, "Cutminus"]
#     ]
#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     dim = 2 # dimension of the entity
#     tag = surfaces # tag of the entity
#     ID = 8
#     gp: list[tuple[int, list[int]], int, str] = [
#         [(dim, tag), ID, "Domain"],
#     ]

#     TPZMeshModeling.CreatePhysicalGroup(gp)

#     #%% Setting the mesh

#     # domain lines
#     domain_lines: list[int] = [i + 1 for i in range(37)] 

#     # creating the transfinte lines and surfaces
#     TPZMeshModeling.TransfiniteCurve(domain_lines, nNodes)
#     TPZMeshModeling.TransfiniteSurface(surfaces, "Left") # it really doesn't matter whether you use "Left" or "Right"

#     # recombining elements to have quadrilateral ones
#     TPZMeshModeling.RecombineSurface(surfaces)

#     # generating mesh
#     TPZMeshModeling.CreateMesh(mesh_dim) 

#     # showing model on gmsh
#     # TPZMeshModeling.ShowModel()

#     # writing .msh file
#     TPZMeshModeling.WriteMeshFiles(file_name, ".msh")

#     path = "/home/cordeiro/projects/ErrorEstimation/build/Projects/ErrorNaca/"
#     TPZMeshModeling.MoveFiles(file_name, None, path, ".msh")

#     # closing gmsh
#     TPZMeshModeling.End()

# #%% Defining main function
# if __name__ == "__main__":
#     main()
