from dataclasses import dataclass
from typing import ClassVar
import json
import gmsh
import sys
import os

@dataclass
class TPZMeshModeling:
    kernel: ClassVar[str] = 'occ'

    @staticmethod
    def PrintJson(JsonData: dict, FileName: str)->None:
        """ 
        Reads the information from the 'JsonData' dictionary and writes them in a Json file
        named 'FileName'.
        """
        json_object = json.dumps(JsonData, indent=4)

        with open(FileName+".json", "w") as outfile:
            outfile.write(json_object)

    @staticmethod
    def WriteMeshFiles(FileName: str, *extensions: str)->None:
        """
        Creates files with given extensions from a gmsh model.
        Currently are available, .geo, .msh, and .vtk 
        """
        for extension in extensions:
            gmsh.write(FileName + extension)

    @staticmethod
    def MoveFiles(fileName: str, JsonNewPath: str = None, MeshNewPath:str = None, *extensions: str)->None:
        """
        Move the files with extensions provided to a given directory
        Currently are available: .json, .geo, .msh, and .vtk 
        """
        for extension in extensions:
            if extension == ".json":
                oldJson = fileName + extension
                os.rename(oldJson, JsonNewPath+oldJson)

            else:
                old_file = fileName + extension
                os.rename(old_file, MeshNewPath+old_file)

    @staticmethod
    def CreatePoints(PointCoordinates: list[int], lc: float)->list[int]:
        """
        Return a list with the tags of the points created from 'PointCoordinates' with mesh size lc. 
        """
        points = []
        for coord in PointCoordinates:
            x,y,z = coord
            if TPZMeshModeling.kernel == 'occ': 
                p = gmsh.model.occ.addPoint(x,y,z,lc)

            elif TPZMeshModeling.kernel == 'built':
                p = gmsh.model.geo.addPoint(x,y,z,lc)
            
            points.append(p)

        return points

    @staticmethod
    def CreateLines(LineIndexes: list[int])->list[int]:
        """
        Returns a list with the tags of the lines created from the indexes of points in 'LineIndexes'
        """
        lines = []
        for index in LineIndexes:
            init, end = index
            if TPZMeshModeling.kernel == 'occ':
                l = gmsh.model.occ.addLine(init, end)

            elif TPZMeshModeling.kernel == 'built':
                l = gmsh.model.geo.addLine(init, end)

            lines.append(l)

        return lines

    @staticmethod
    def CreateCurveLoops(CurveLoopIndexes: list[int])->list[int]:
        """
        Returns a list of curve loops created from the line indexes in 'CurveLoopIndexes'
        """
        curves = []
        for  index in (CurveLoopIndexes):
            if TPZMeshModeling.kernel == 'occ':
                c = gmsh.model.occ.addCurveLoop(index)

            elif TPZMeshModeling.kernel == 'built':
                c = gmsh.model.geo.addCurveLoop(index)

            curves.append(c)

        return curves

    @staticmethod
    def CreatePlanes(PlaneIndexes: list[int])->list[int]:
        """
        Returns a list of integers as the created plane tags. 
        """
        planes = []
        for index in (PlaneIndexes):
            if TPZMeshModeling.kernel == 'occ': 
                plane = gmsh.model.occ.addPlaneSurface(index) 

            elif TPZMeshModeling.kernel == 'built':
                plane = gmsh.model.geo.addPlaneSurface(index) 

            planes.append(plane)

        return planes

    @staticmethod
    def CreateSurfaceLoop(SurfaceLoopIndexes: list[int])->list[int]:
        """
        Returns a list of integers as the created surface loop tags
        """
        surface_loop = []
        for index in (SurfaceLoopIndexes):
            if TPZMeshModeling.kernel == 'occ':
                s = gmsh.model.occ.addSurfaceLoop(index)

            elif TPZMeshModeling.kernel == 'built':
                s = gmsh.model.geo.addSurfaceLoop(index)

            surface_loop.append(s)

        return surface_loop

    @staticmethod
    def CreateVolumes(VolumesIndexes: list)->None:
        """
        Returns a list of integers as the created volume tags
        """

        volume = []
        for index in (VolumesIndexes):
            if TPZMeshModeling.kernel == 'occ':
                v = gmsh.model.occ.addVolume([index])

            elif TPZMeshModeling.kernel == 'built':
                v = gmsh.model.geo.addVolume([index])

            volume.append(v)

        return volume

    @staticmethod
    def CreatePhysicalGroup(GroupData: list)->None:
        """
        Creates the physical groups from a list contaning the group information. 
        Proide a tuple with 'dimTags', an integer ID, and the group name.
        """

        for dimTag, id, name in GroupData:
            dimension, tag = dimTag
            gmsh.model.addPhysicalGroup(dimension, tag, tag=id, name=name)

    @staticmethod
    def CreateCircleArcs(arc_points: list[int])->list[int]:
        arcs = []
        for coord in arc_points:
            start, center, end = coord

            if TPZMeshModeling.kernel == 'occ':
                arc = gmsh.model.occ.addCircleArc(start, center, end)

            elif TPZMeshModeling.kernel == 'built':
                arc = gmsh.model.occ.addCircleArc(start, center, end)

            arcs.append(arc)

        return arcs

    @staticmethod
    def CreateCircles(Xcenter: float, Ycenter: float, Zcenter: float, Radius: float) -> int:
        """
        Returns the tag of the created surface circle, using the OPENCASCADE kernel
        """

        circle = gmsh.model.occ.addCircle(Xcenter,Ycenter,Zcenter,Radius)
        curveLoopCircle = gmsh.model.occ.addCurveLoop([circle])
        surfaceCircle = gmsh.model.occ.addPlaneSurface([curveLoopCircle])
        
        return surfaceCircle

    @staticmethod
    def CreateRectangles(coordinates:tuple[float], sideX:float, sideY:float)->list[int]:
        """
        Returns a list of integers as the rectangle surface tags, using the OPENCASCADE kernel
        """
        square_list = []
        for coord in coordinates:
            x, y, z = coord
            square = gmsh.model.occ.addRectangle(x, y, z, sideX, sideY)
            square_list.append(square)

        return square_list

    @staticmethod
    def MakeHoles(object: int, holesList: list, holeDim: int)->list[int]:
        """
        Makes holes in a surface domain. Given a 'domain' tag, the 'holeList' tags, and the 'meshDim', it uses the 
        gmsh module cut to calculate the boolean difference the object domain and the object to be cut from it. Returns
        the new surface tags.
        """
        
        holesTuple = [(holeDim, hole) for hole in holesList]
        holes = gmsh.model.occ.cut([(holeDim,object)], holesTuple)

        return holes

    @staticmethod
    def TurnOnRepresentation(*variables: str):
        """
        Turn on the selected entity CAD representation. 
            - points
            - curves
            - surfaces
            - volumes
        """
        for var in variables:
            var = var.capitalize()
            gmsh.option.setNumber("Geometry." + var, 1)

    def TurnOnLabels(*variables: str):
        """
        Turn on the selected entities' labels. 
            - points
            - curves
            - surfaces
            - volumes
        """
        for var in variables:
            var = var.capitalize()  + "Numbers"
            gmsh.option.setNumber("Geometry." + var, 1)

    @staticmethod
    def TurnOnNormals(size: int=50)->None:
        """
        Display the normal vectors with 'size'
        """
        gmsh.option.setNumber("Geometry.Normals", size)

    @staticmethod
    def TurnOnTangents(size: int=50)->None:
        """
        Display the tnagent vectors with 'size'
        """
        gmsh.option.setNumber("Geometry.Tangents", size)

    @staticmethod
    def Synchronize():
        """
        Synchronizes the gmsh CAD with the gmsh kernel
        """
        if TPZMeshModeling.kernel == 'occ':
            gmsh.model.occ.synchronize()
        
        elif TPZMeshModeling.kernel == 'built':
            gmsh.model.geo.synchronize()

    @staticmethod
    def ShowModel(meshDim: int =-1)->None:
        """
        Show the model
        """
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run() 

    @staticmethod
    def SetMeshSize(fieldID, surfacesList, meshSize, bigNumber = 1e12):
        """
        Defines the mesh size in a surface.
        """
        fieldOperator = gmsh.model.mesh.field 

        fieldOperator.add("Constant", fieldID)
        fieldOperator.set_number(fieldID, "IncludeBoundary",1)
        fieldOperator.set_numbers(fieldID, "SurfacesList", surfacesList)
        fieldOperator.set_number(fieldID, "VIn", meshSize)
        fieldOperator.set_number(fieldID, "VOut", bigNumber)

        fieldOperator.setAsBackgroundMesh(fieldID)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    @staticmethod
    def Begin():
        """
        Initializes gmsh
        """
        gmsh.initialize()

    @staticmethod
    def End():
        """
        Finalizes gmsh
        """
        gmsh.finalize()

    @staticmethod
    def SetGmshKernel(kernel: int)->None:
        """
        Sets the gmsh kernel used during the modeling. Default: OPENCASCADE
        0 -> built-in
        1 -> OPENCASCADE
        """
        option = {0: 'built', 1: 'occ'}
        TPZMeshModeling.kernel = option[kernel]

    @staticmethod
    def CreateMesh(meshDim: int):
        """
        Creates the model mesh with dimension meshDim
        """
        gmsh.model.mesh.generate(meshDim)

    @staticmethod
    def TransfiniteCurve(curves: list[int], nPoints: int)->None:
        """
        Transfinite the curves in the list with nPoints
        """
        for curve in curves:
            gmsh.model.mesh.setTransfiniteCurve(curve, nPoints)

    @staticmethod
    def TransfiniteSurface(surfaces: list[int], transfiniteType: str)->None:
        """
        Transfinite the surfaces in the list with the transfiniteType
        """
        for surface in surfaces:
            gmsh.model.mesh.setTransfiniteSurface(surface, transfiniteType)

    @staticmethod
    def RecombineSurface(surfaces: list[int])->None:
        """
        Recombine the surfaces in the list
        """
        for surface in surfaces:
            gmsh.model.mesh.setRecombine(2, surface)