# transforms form of bounds of needed
# applies boundary conditions
# Dirichlet, Neumann, Robin

from enum import Enum

import numpy as np
import scipy
import sympy as sp

import domain
import local

x = sp.symbols('x')
y =  sp.symbols('y')
ksi = sp.symbols('ksi')
eta = sp.symbols('eta')

class Type(Enum):
    DIRICHLET = 1
    NEUMANN = 2
    ROBIN = 3

class BoundaryConds:
    def applyDirichlet(self, matrix, vector, nodes, vertices, g):
        for node in nodes:
            coords = vertices[node]
            matrix[2*node] = [0 for el in matrix[2*node]]
            matrix[2*node][2*node] = 1
            vector[2*node] = g(coords[0], coords[1])[0]

            matrix[2*node+1] = [0 for el in matrix[2*node+1]]
            matrix[2*node+1][2*node+1] = 1
            vector[2*node+1] = g(coords[0], coords[1])[1]
        return matrix, vector

    def applyNeumann(self, NT, vector, p, nodes_at_bound, elems_at_bound):
        # system = local.LTriangle.system(element.vertices)
        # NT (ksi, eta) ~> NT (x, y) because we integrate over bound
         # # NT*p = [ 6x2 ] * [ 2 ] = [ 6 ]
        # NTp = np.matvec(NT, p(x, y))

        t = sp.symbols('t')



        for i in range(len(elems_at_bound)):
            elem = elems_at_bound[i]
            ios = np.where(elem.nodes==nodes_at_bound[i])[0][0]
            start = local.local_triangle[ios]
            ioe = np.where(elem.nodes==nodes_at_bound[i+1])[0][0]
            end = local.local_triangle[ioe]

            system = local.LTriangle.system(elem.vertices)
            pksi = sp.sympify(p(x, y)[0]).subs('x', system[0]).subs('y', system[1])
            peta = sp.sympify(p(x, y)[1]).subs('x', system[0]).subs('y', system[1])
            NTp = np.matvec(NT, [pksi, peta])

            ksi_t = start[0] + (end[0]-start[0])*t
            eta_t = start[1] + (end[1]-start[1])*t
            NTp = [el.subs('ksi', ksi_t).subs('eta', eta_t) for el in NTp]

            integral = [elem.jacobian * scipy.integrate.quad(sp.lambdify(t, el), 0, 1)[0] for el in NTp]
            ion = nodes_at_bound[i+1]
            vector[2*ion] += integral[2*ioe]
            vector[2*ion+1] += integral[2*ioe+1]
        return vector

    def applyRobin(self, matrix, vector, nodes, beta, sigma, g):
        # beta * Nu + sigma * u = g

        return matrix, vector