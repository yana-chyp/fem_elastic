from abc import ABC, abstractmethod

import numpy as np
import sympy as sp
import scipy.integrate as scin
from sympy.core import numbers

local_triangle = [[0,0], [1,0], [0,1]]
local_det = 1

ksi = sp.symbols('ksi')
eta = sp.symbols('eta')
x = sp.symbols('x')
y = sp.symbols('y')

class Local(ABC):
    @abstractmethod
    def system(self, element):
        pass
    def jacobi_sys(self, system):
        pass
    def jacobi(self, figure):
        # for squares?
        pass
    def base(self, approx_degree=1):
        pass
    @abstractmethod
    def matrix(self):
        pass
    def integrate(self, f):
        pass


class LTriangle(Local):
    @staticmethod
    def system(triangle):
        x = triangle[0][0] + (triangle[1][0] - triangle[0][0]) * ksi + (triangle[2][0] - triangle[0][0]) * eta
        y = triangle[0][1] + (triangle[1][1] - triangle[0][1]) * ksi + (triangle[2][1] - triangle[0][1]) * eta
        return (x, y)

    @staticmethod
    def ksieta_from_xy(triangle, x, y):
        a = [[triangle[1][0] - triangle[0][0], triangle[2][0] - triangle[0][0]],
                  [triangle[1][1] - triangle[0][1], triangle[2][1] - triangle[0][1]]]
        b = [x - triangle[0][0],
            y - triangle[0][1]]
        return np.linalg.solve(a, b)

    def jacobi_sys(self, system):
        x = system[0]
        y = system[1]
        jacobi = [[sp.diff(x, 'ksi'), sp.diff(y, 'ksi')], [sp.diff(x, 'eta'), sp.diff(y, 'eta')]]
        return jacobi[0][0] * jacobi[1][1] - jacobi[0][1] * jacobi[1][0]

    @staticmethod
    def jacobi(triangle):
        det = (triangle[1][0] - triangle[0][0]) * (triangle[2][1] - triangle[0][1])
        det -= (triangle[2][0] - triangle[0][0]) * (triangle[1][1] - triangle[0][1])
        return det

    @staticmethod
    def base(approx_degree=1):
        if approx_degree == 1:
            base = [1 - ksi - eta, ksi, eta]
        elif approx_degree == 2:
            base = [(1 - ksi - eta) * (1 - 2 * ksi - 2 * eta),
                    ksi * (2 * ksi - 1),
                    eta * (2 * eta - 1),
                    4 * ksi * (1 - ksi - eta),
                    4 * ksi * eta,
                    4 * eta * (1 - ksi - eta)]
        elif approx_degree == 3:
            # todo
            base = []
        else:
            raise ('degree not impemented')
        return base

    def matrix(self):
        return [[]]

    @staticmethod
    def integrate(f):
        if not callable(f):
            const_val = float(f) if isinstance(f, numbers.Number) else f
            f = lambda ksi, eta: const_val
        # f = f(ksi, eta)
        return scin.dblquad(lambda ksi, eta: f(ksi, eta), 0, 1, lambda ksi: 0, lambda ksi: 1-ksi)[0]

# class LSquare(Local):
#     # def __init__(self):
#     def system(self, polyfour):
#
#     def jacobi_sys(self, system):
#
#     # def jacobi(self, figure):
#
#     def base(self, approx_degree=1):


def f(ksi, eta):
    return ksi+eta

# ltr = LTriangle()
# triangle = [[0, 0], [2, 1], [-1, 2]]
# system = ltr.system(triangle)
# print('system:')
# print(system)
#
# print('jacobi from system: ', local_jacobi_sys(system))
# print('jacobi from triangle: ', local_jacobi_tri(triangle))

# f_lamb = sp.lambdify([ksi, eta], f(ksi, eta))
# res = scin.dblquad(lambda ksi, eta: f_lamb(ksi, eta), 0, 1, lambda ksi: 0, lambda ksi: 1-ksi)[0]
# print(res)