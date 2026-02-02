# assembles system (matrix) of elements
# independently of their kind
import elasticdeform as elde
import local
from elasticdeform import ElasticDeform
from local import LTriangle



class System:
    elde: ElasticDeform
    def __init__(self, elde):
        self.elde = elde
    def assemble_matrix(self, elements, k_e, nn):
        # nn - number of nodes
        matrix = [[0 for j in range(2*nn)] for i in range(2*nn)]
        for element in elements:
            #supposed to be equal to len of stiffness matrix - no
            m = len(element.nodes)
            for i in range(m):
                # ion - index of node
                ion = element.nodes[i]
                for j in range(m):
                    # jon index of another node
                    jon = element.nodes[j]
                    matrix[2*ion][2*jon] += element.jacobian * k_e[2*i][2*j]
                    matrix[2*ion+1][2*jon+1] += element.jacobian * k_e[2*i+1][2*j+1]
        return matrix

    def assemble_vector(self, elements, nn, b):
        # elde = ElasticDeform(base)
        vector = [0 for i in range(2*nn)]
        for element in elements:
            m = len(element.nodes)
            r_e = self.elde.load_vector(b, element)
            for i in range(m):
                ion = element.nodes[i]
                vector[2*ion] += element.jacobian * r_e[2*i]
                vector[2*ion+1] += element.jacobian * r_e[2*i+1]
        return vector


# ed = elde.ElasticDeform(base=LTriangle.base())
# k_e = ed.stiffness_matrix()
# print(k_e)
# is [ 6x6 ] for 3-node triangle

