# triangles or squares
# stores nodes (numbers), probably jacobian so on
# processes element
# converts it to the form to assemble system of
from enum import Enum
import numpy as np
import finel
import local
import triangulation as tr

class TypeOfElement(Enum):
    TRIANGLE = 3
    SQUARE = 4

class Element:
    type: TypeOfElement
    nodes: np.array(int)
    vertices: np.array(object)
    ksieta_coords: np.array(object)
    jacobian: float
    approx: int

    def __init__(self, type, nodes, vertices, approx=1):
        if (type!=3 and type!=4):
            raise('wrong type of element')
        self.type = type
        self.nodes = nodes
        self.vertices = []
        for node in self.nodes:
            self.vertices.append(vertices[node])
        self.approx = approx
        if (type==3):
            self.calculate_jacobian(local.LTriangle)
        else:
            self.calculate_jacobian(local.LSquare)

    def calculate_jacobian(self, localFigure):
        self.jacobian = localFigure.jacobi(self.vertices)

    def calculate_stress(self, D, strain):
        return D*strain

    def calculate_strain(self, B, u_e):
        return B*u_e
