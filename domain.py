# use delaunay for triangulation
# smth else for squares
# stores array of elements of nodes
# maps coordinates to numbers of nodes
import math
from cmath import cos, sin


class Domain:
    def __init__(self, vertices):
        self.vertices = vertices
        self.bounds = []
        for i in range(len(self.vertices) - 1):
            self.bounds.append([self.vertices[i], self.vertices[i+1]])
        self.bounds.append([self.vertices[0], self.vertices[len(self.vertices) - 1]])

    def get_bounds(self):
        return self.bounds

    @staticmethod
    def check_on_edge(start, end, coords):
        tol = 0.001
        return abs((coords[1] - start[1]) * (end[0] - start[0]) - (end[1] - start[1]) * (coords[0] - start[0])) < tol

    def find_nodes_at_bound(self, bound, vertices):
        nodes = []
        for i in range(len(vertices)):
            if self.check_on_edge(bound[0], bound[1], vertices[i]):
                nodes.append(i)
        nodes = sorted(nodes, key= lambda i: (vertices[i][0]-bound[0][0])**2 + (vertices[i][1]-bound[0][1])**2)
        return nodes
    def find_elems_at_bound(self, nodes_at_bound, elements):
        belems = []
        for i in range(len(nodes_at_bound)-1):
            elem = [el for el in elements if (nodes_at_bound[i] in el.nodes and nodes_at_bound[i+1] in el.nodes)]
            assert(len(elem)==1)
            belems.append(elem[0])
        return belems

class CircularDomain:
    # vertices:
    def __init__(self, vertices, density=10):
        self.bounds = []
        for i in range(len(vertices) - 1):
            self.bounds.append([vertices[i], vertices[i+1]])
        self.bounds.append([vertices[0], vertices[len(vertices) - 1]])

        self.vertices = []
        self.vertices.append(vertices[0])
        self.vertices.append(vertices[1])
        a = math.pi/(2*density)
        self.r1 = vertices[1][1]-vertices[2][1]
        self.c = [vertices[1][0], vertices[1][1]-self.r1]
        for i in range(1, density):
            self.vertices.append([self.r1*math.sin(a*i)+self.c[0], self.r1*math.cos(a*i)+self.c[1]])
        self.vertices.append(vertices[2])
        self.vertices.append(vertices[3])
        self.r2 = vertices[0][1]-vertices[3][1]
        for i in range(1, density):
            self.vertices.append([self.r2 * math.cos(a * i) + self.c[0], self.r2 * math.sin(a * i) + self.c[1]])

    def get_bounds(self):
        return self.bounds

    def find_elems_at_bound(self, nodes_at_bound, elements):
        # print(nodes_at_bound)
        belems = []
        for i in range(len(nodes_at_bound)-1):
            elem = [el for el in elements if (nodes_at_bound[i] in el.nodes and nodes_at_bound[i+1] in el.nodes)]
            # print(nodes_at_bound[i], nodes_at_bound[i+1])
            assert(len(elem)==1)
            belems.append(elem[0])
        return belems

    def find_nodes_at_bound(self, bound, vertices):
        nodes = []
        for i in range(len(vertices)):
            if self.check_on_edge(bound[0], bound[1], vertices[i]):
                nodes.append(i)
        nodes = sorted(nodes, key= lambda i: (vertices[i][0]-bound[0][0])**2 + (vertices[i][1]-bound[0][1])**2)
        return nodes

    def check_on_edge(self, start, end, coords):
        tol_lin = 0.1
        tol_circ = 0.05
        if (start==self.bounds[0][0] and end==self.bounds[0][1]) or (start==self.bounds[2][0] and end==self.bounds[2][1]):
            return abs((coords[1] - start[1]) * (end[0] - start[0]) - (end[1] - start[1]) * (coords[0] - start[0])) < tol_lin
        elif (start==self.bounds[0][1] and end==self.bounds[2][0]):
            return abs((coords[0]-self.c[0])**2 + (coords[1]-self.c[1])**2 - self.r1**2)/self.r1**2 < tol_circ
        else:
            return abs((coords[0] - self.c[0]) ** 2 + (coords[1] - self.c[1]) ** 2 - self.r2 ** 2)/self.r2**2 < tol_circ

import element

class Mesh:
    # points stores coordinates and numbers as indices
    def __init__(self, type, domain, approx_degree, density_params):
        self.type = type
        self.points = self.generate_mesh(domain, approx_degree, density_params)

    # def generate_mesh(self, domain, approx_degree, density_params):
