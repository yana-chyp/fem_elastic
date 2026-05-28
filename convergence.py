import math

import numpy as np

import domain
from padapter import PAdapter
from solver import apply_bounds_lame, assemble_matvec, create_elements, init_approx_val_grads_lame, init_exact_val_grads_lame, initialize_base, plot_mesh, plot_solution, solve
import triangulation as tr
import system as sys
import bound_cond as bc
from norms import h1_error_elements

domain = domain.CircularDomain([[0, 1], [0, 2], [2, 0], [1, 0]], 4)

alpha=1000
young=2500
degree = 1

def b(x, y):
    return [0, 0]
def f1(x, y):
    return [0, 0]
def f2(x, y):
    return [0, 0]
def g1(x, y):
    return [0, 0]
def g2(x, y):
    return [x, y]

tr = tr.Triangulator(domain)
min_angle = 30

max_area_2 = 0.05
mesh = tr.triangulate(len(domain.vertices), min_angle, max_area_2)
vertices = mesh['vertices']
triangles = mesh['triangles']
plot_mesh(vertices, triangles)
print('vertices: ', len(vertices))
print('triangles: ', len(triangles))
elements, vertices = create_elements(vertices, triangles, degree)
ed, base, k_e, m_e, base_integrals = initialize_base(young)
s = sys.System(ed)
matrix, vector = assemble_matvec(s, k_e, elements, len(vertices), b)
bound_types = [bc.Type.ROBIN, bc.Type.NEUMANN, bc.Type.ROBIN, bc.Type.NEUMANN]
matrix, vector = apply_bounds_lame(domain, bound_types, elements, vertices, [lambda x, y: [0, 0], g1, lambda x, y: [0, 0], g2],
                            base, m_e, base_integrals, matrix, vector, alpha)
u = solve(matrix, vector)
u_exact, u_exact_grad = init_exact_val_grads_lame(ed, elements, vertices)
u_approx, u_approx_grad = init_approx_val_grads_lame(u, elements)
plot_solution(u, elements, vertices, u_exact)

errors1 = h1_error_elements(elements, u_exact, u_approx, u_exact_grad, u_approx_grad)
sum005 = sum(errors1)



max_area = 0.1
mesh = tr.triangulate(len(domain.vertices), min_angle, max_area)
vertices = mesh['vertices']
triangles = mesh['triangles']
plot_mesh(vertices, triangles)
print('vertices: ', len(vertices))
print('triangles: ', len(triangles))
elements, vertices = create_elements(vertices, triangles, degree)
ed, base, k_e, m_e, base_integrals = initialize_base(young)
s = sys.System(ed)
matrix, vector = assemble_matvec(s, k_e, elements, len(vertices), b)
bound_types = [bc.Type.ROBIN, bc.Type.NEUMANN, bc.Type.ROBIN, bc.Type.NEUMANN]
matrix, vector = apply_bounds_lame(domain, bound_types, elements, vertices, [lambda x, y: [0, 0], g1, lambda x, y: [0, 0], g2],
                            base, m_e, base_integrals, matrix, vector, alpha)
u = solve(matrix, vector)
u_exact, u_exact_grad = init_exact_val_grads_lame(ed, elements, vertices)
u_approx, u_approx_grad = init_approx_val_grads_lame(u, elements)
plot_solution(u, elements, vertices, u_exact)

errors2 = h1_error_elements(elements, u_exact, u_approx, u_exact_grad, u_approx_grad)
sum02 = sum(errors2)

print(f'error for {max_area} is {sum02}\nerror for {max_area_2} is {sum005}')
p = (math.log(sum02) - math.log(sum005))/math.log(2)
print(f'convergence rate: {p}')


max_area2 = 0.2
mesh = tr.triangulate(len(domain.vertices), min_angle, max_area2)
vertices = mesh['vertices']
triangles = mesh['triangles']
plot_mesh(vertices, triangles)
print('vertices: ', len(vertices))
print('triangles: ', len(triangles))
elements, vertices = create_elements(vertices, triangles, degree)
ed, base, k_e, m_e, base_integrals = initialize_base(young)
s = sys.System(ed)
matrix, vector = assemble_matvec(s, k_e, elements, len(vertices), b)
bound_types = [bc.Type.ROBIN, bc.Type.NEUMANN, bc.Type.ROBIN, bc.Type.NEUMANN]
matrix, vector = apply_bounds_lame(domain, bound_types, elements, vertices, [lambda x, y: [0, 0], g1, lambda x, y: [0, 0], g2],
                            base, m_e, base_integrals, matrix, vector, alpha)
u = solve(matrix, vector)
u_exact, u_exact_grad = init_exact_val_grads_lame(ed, elements, vertices)
u_approx, u_approx_grad = init_approx_val_grads_lame(u, elements)
plot_solution(u, elements, vertices, u_exact)

errors3 = h1_error_elements(elements, u_exact, u_approx, u_exact_grad, u_approx_grad)
sum08 = sum(errors3)
print(f'error for {max_area2} is {sum08}')

p = math.log2((sum02 - sum08)/(sum005 - sum02))
print(f'convergence rate by eitkin: {p}')
