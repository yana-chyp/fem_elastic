import sympy as sp
import numpy as np

import local

x = sp.symbols('x')
y =  sp.symbols('y')
ksi = sp.symbols('ksi')
eta = sp.symbols('eta')

class ElasticDeform:
    def __init__(self, base):
        # for degree = 1
        E = 1
        nu = 0
        # plane stress
        self.D = [[E/(1-nu**2) * element for element in row] for row in [[ 1, nu,        0],
                                                                    [nu,  1,        0],
                                                                    [ 0,  0, (1-nu)/2]]]
        L = [[self.d_dksi,         0],
             [        0, self.d_deta],
             [self.d_deta, self.d_dksi]]
        self.N = [[base[0],       0, base[1],       0, base[2],       0],
             [      0, base[0],       0, base[1],       0, base[2]]]
        # B = L*N
        self.B = [[self.d_dksi(base[0]),                  0, self.d_dksi(base[1]),                  0, self.d_dksi(base[2]),                  0],
             [                 0, self.d_deta(base[0]),                  0, self.d_deta(base[1]),                  0, self.d_deta(base[2])],
             [self.d_deta(base[0]), self.d_dksi(base[0]), self.d_deta(base[1]), self.d_dksi(base[1]), self.d_deta(base[2]), self.d_dksi(base[2])]]
        # BT = np.transpose(self.B)
        # BT * D = [ 6x3 ]
        # BT * D * B = [ 6x6 ]
        # K_e = integral(BT*D*B)_over(V_e) = [ 6x6 ]

        # sigma = D * epsilon

    def calculate_stress(self, D, strain):
        # sigma = D * strain = [ 3x3 ] * [ 3 ] = [ 3 ]
        return D*strain

    def calculate_strain(self, B, u_e):
        # epsilon = B * u_e = [ 3x6] * [ 6 ] = [ 3 ]
        return B*u_e

    def d_dksi(self, func):
        return sp.diff(func, ksi)
    def d_deta(self, func):
        return sp.diff(func, eta)


    def stiffness_matrix(self):
        BT = np.transpose(self.B)
        # BT*D
        # n = len(BT)
        # m = len(self.D[0])
        # k = len(BT[0])
        prod = np.matmul(np.matmul(BT, self.D), self.B)
        k_e = [[local.LTriangle.integrate(el) for el in row] for row in prod]

        return k_e

    def load_vector(self, b, element):
        # r_e = integral(NT*b)_over(V_e) + integral(NT*p)_over(S_e)
        # b = [b_x, b_y]; p = [p_x, p_y]
        system = local.LTriangle.system(element.vertices)
        bksi = sp.sympify(b(x, y)[0]).subs('x', system[0]).subs('y', system[1])
        beta = sp.sympify(b(x, y)[1]).subs('x', system[0]).subs('y', system[1])
        # NT*b = [ 6x2 ] * [ 2 ] = [ 6 ]
        NT = np.transpose(self.N)
        NTb = np.matvec(NT, [bksi, beta])
        r_e = [local.LTriangle.integrate(sp.lambdify([ksi, eta], el)) for el in NTb]

        return r_e

# base = local.LTriangle.base()
# eldef = ElasticDeform(base)
# k_e = eldef.stiffness_matrix()
#
# for row in k_e:
#     print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )
