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
            # print(elem.nodes)
            ios = np.where(elem.nodes==nodes_at_bound[i])[0][0]
            start = local.local_triangle[ios]
            ioe = np.where(elem.nodes==nodes_at_bound[i+1])[0][0]
            end = local.local_triangle[ioe]

            system = local.LTriangle.system(elem.vertices)
            # print(system)
            pksi = sp.sympify(p(x, y)[0])
            peta = sp.sympify(p(x, y)[1])
            # print("pksi = ", pksi, ", peta =  ", peta)
            pksi = pksi.subs('x', system[0]).subs('y', system[1])
            peta = peta.subs('x', system[0]).subs('y', system[1])
            # print("pksi = ", pksi, ", peta =  ", peta)
            NTp = np.matvec(NT, [pksi, peta])
            # print("NTp = ", NTp)

            ksi_t = start[0] + (end[0]-start[0])*t
            eta_t = start[1] + (end[1]-start[1])*t
            NTp = [el.subs('ksi', ksi_t).subs('eta', eta_t) for el in NTp]
            # print("NTp = ", NTp)
            integral = [elem.jacobian * scipy.integrate.quad(sp.lambdify(t, el), 0, 1)[0] for el in NTp]
            ion = nodes_at_bound[i+1]
            vector[2*ion] += integral[2*ioe]
            vector[2*ion+1] += integral[2*ioe+1]
        return vector

    def applyRobin(self, NT, matrix, vector, m_e, nodes_at_bound, elems_at_bound, alpha, u_0):
        # Nu + alpha * u = alpha * u_0(x)

        vector = self.applyNeumann(NT, vector,  lambda x, y: [alpha*u_0(x, y)[0], alpha*u_0(x, y)[1]], nodes_at_bound, elems_at_bound)

        t = sp.symbols('t')
        for i in range(len(elems_at_bound)):
            elem = elems_at_bound[i]
            # print(elem.nodes)
            ios = np.where(elem.nodes == nodes_at_bound[i])[0][0]
            start = local.local_triangle[ios]
            ioe = np.where(elem.nodes == nodes_at_bound[i + 1])[0][0]
            end = local.local_triangle[ioe]
            ksi_t = start[0] + (end[0] - start[0]) * t
            eta_t = start[1] + (end[1] - start[1]) * t
            m_e = [[el.subs('ksi', ksi_t).subs('eta', eta_t) for el in row] for row in m_e]
            integral = [[scipy.integrate.quad(sp.lambdify(t, el), 0, 1)[0] for el in row] for row in m_e]
            print(elem.jacobian)
            ion = nodes_at_bound[i]
            # jon = nodes_at_bound[i]
            for j in range(len(elems_at_bound[i].nodes)):
                jon = elems_at_bound[i].nodes[j]
                matrix[2 * ion][2 * jon] += alpha * integral[2 * i][2 * j]
                matrix[2 * ion + 1][2 * jon + 1] += alpha * integral[2 * i + 1][2 * j + 1]

        return matrix, vector