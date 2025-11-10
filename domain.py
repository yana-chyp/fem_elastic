# use delaunay for triangulation
# smth else for squares
# stores array of elements of nodes
# maps coordinates to numbers of nodes

class Domain:
    def __init__(self, vertices ):
        # store bounds as well?
        self.vertices = vertices

    def get_bounds(self):
        bounds = []
        for i in range(len(self.vertices) - 1):
            bounds.append([self.vertices[i], self.vertices[i+1]])
        bounds.append([self.vertices[0], self.vertices[len(self.vertices) - 1]])
        return bounds

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

import element

class Mesh:
    # points stores coordinates and numbers as indices
    def __init__(self, type, domain, approx_degree, density_params):
        self.type = type
        self.points = self.generate_mesh(domain, approx_degree, density_params)

    # def generate_mesh(self, domain, approx_degree, density_params):
