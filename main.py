import numpy as np
from matplotlib import pyplot as plt

import triangulation as tr
import element
import elasticdeform as elde
import local
import system as s
import domain
import bound_cond as bc

def b(x, y):
    return [0, 0]

# def p(x, y):
#     return [0, 0]

domain = domain.Domain([[-2, 0], [-2, 1], [2, 1], [2, 0]])
tr = tr.Triangulator(domain)
min_angle = 30
max_area = 0.5
mesh = tr.triangulate(min_angle, max_area)
tr.plot(mesh)
print(mesh)
vertices = mesh['vertices']
triangles = mesh['triangles']

elements = []
for triangle in triangles:
    elements.append(element.Element(3, triangle, vertices))

base = local.LTriangle.base()
ed = elde.ElasticDeform(base)
NT = np.transpose(ed.N)
k_e = ed.stiffness_matrix()
# print(k_e)

matrix = s.System.assemble_matrix(elements, k_e, len(vertices))
# for row in matrix:
#     print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )

# r_e = ed.load_vector(b, p, elements[0])
# print('load vector r_e: ')
# print(r_e)

vector = s.System.assemble_vector(elements, len(vertices), base, b)
print('global load vector: ')
for el in vector:
    print(el)

bounds = domain.get_bounds()
bound_cond_types = [bc.Type.DIRICHLET, bc.Type.NEUMANN, bc.Type.DIRICHLET, bc.Type.NEUMANN]
bound_nodes = []
bound_elems = []
for bound in bounds:
    # print('bound: ', bound)
    nodes_at_bound = domain.find_nodes_at_bound(bound, vertices)
    # print('nodes at bound: ', nodes_at_bound)
    bound_nodes.append(nodes_at_bound)
    bound_elems.append(domain.find_elems_at_bound(nodes_at_bound, elements))
    # print('elems at bound')
    # for el in bound_elems[-1]:
    #     print(el.nodes)

def f1(x, y):
    return [0, 0]

def f2(x, y):
    return [0, 0]

def g1(x, y):
    return [0, 10]

def g2(x, y):
    return [0, -10]

bconds = bc.BoundaryConds()

bconds.applyNeumann(NT, vector, g1, bound_nodes[1], bound_elems[1])
bconds.applyNeumann(NT, vector, g2, bound_nodes[3], bound_elems[3])
bconds.applyDirichlet(matrix, vector, bound_nodes[0], vertices, f1)
bconds.applyDirichlet(matrix, vector, bound_nodes[2], vertices, f2)


for row in matrix:
    print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )
for el in vector:
    print(el)

u = np.linalg.solve(matrix, vector)
print('u = ')
print('[' + ', '.join([f"{u[i]:.4f}" for i in range(len(u)) if i%2==0]) + ']' )
print('[' + ', '.join([f"{u[i]:.4f}" for i in range(len(u)) if i%2==1]) + ']' )




scale = 1.0  # deformation scale factor (tune as needed)
# Extract displacements
u_x = u[0::2]
u_y = u[1::2]

# Compute deformed coordinates
nodes_def = np.copy(vertices)
nodes_def[:, 0] += scale * u_x
nodes_def[:, 1] += scale * u_y

# --- Visualization ---
plt.figure(figsize=(6, 6))

# Plot undeformed mesh (gray)
for elem in triangles:
    elem_nodes = np.append(elem, elem[0])  # close the polygon
    plt.plot(vertices[elem_nodes, 0], vertices[elem_nodes, 1], 'k--', linewidth=0.5)

# Plot deformed mesh (colored)
for elem in triangles:
    elem_nodes = np.append(elem, elem[0])
    plt.plot(nodes_def[elem_nodes, 0], nodes_def[elem_nodes, 1], 'r-', linewidth=1.5)

plt.title("Mesh Deformation Visualization")
plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.legend(["undeformed", "deformed"])
plt.show()