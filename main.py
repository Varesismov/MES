import numpy as np
from math import sqrt

class UniversalElement:
    def __init__(self, num_nodes, node_ids, coordinate_list, bc_list, nodes_number, initial_temp, steptime):
        self.num_nodes = num_nodes
        self.coordinate_list = coordinate_list
        self.coordinate_list_help = []
        self.node_ids = node_ids
        self.node_ids_help = []
        self.bc_list = bc_list
        self.bc_list_help = []
        self.nodes = nodes_number
        self.coordinates_x = []
        self.coordinates_y = []
        self.ksi_values = []
        self.eta_values = []
        self.weights = []
        self.alfa = 300
        self.step_time = steptime

        #DO LICZENIA MACIERZY H
        self.dNdKsi_matrix = []
        self.dNdEta_matrix = []

        self.dXdKsi_help = []
        self.dYdKsi_help = []
        self.dXdEta_help = []
        self.dYdEta_help = []
        self.matrix_jakobian_list = [[] for _ in range(16)]
        self.det = []
        self.det_bc = []

        self.matrix_x_pc = [[] for _ in range(16)]
        self.matrix_y_pc = [[] for _ in range(16)]
        self.matrix_H_pc = [[] for _ in range(16)]
        self.matrix_H_final = []


        #DO LICZENIA MACIERZY HBC
        self.bc_matrix = []
        self.ksi_values_bc = []
        self.eta_values_bc = []
        self.matrix_N_bc = []
        self.matrix_Hbc = [[] for _ in range(self.num_nodes*4)]

        self.dlugosc_scian = []
        self.matrix_Hbc_sum = []
        self.all_hbc_matrix = []
        self.matrix_Hbc_final = []

        #DO LICZENIA WEKTORA P
        self.initial_temp = initial_temp
        self.matrix_P_pre = []
        self.matrix_P_sum = []
        self.matrix_P = []

        #AGREGACJA

        self.lista_2d = [[0] * self.nodes for _ in range(self.nodes)]
        self.lista_2d_hbc = [[0] * self.nodes for _ in range(self.nodes)]
        self.lista_2d_p = [[0] for _ in range(self.nodes)]
        self.lista_2d_c = [[0] * self.nodes for _ in range(self.nodes)]



        #Macierz C
        self.ksi_values_c = []
        self.eta_values_c = []
        self.matrix_N_c = []
        self.matrix_C_pc = [[[0] * 4 for _ in range(4)] for _ in range(16)]
        self.matrix_C_final = [[0] * 4 for _ in range(4)]

        #konwertowanie na 2D

        self.matrix_H_final_2d = [[]]
        self.matrix_Hbc_final_2d = [[]]
        self.matrix_P_2d = [[0] for _ in range(4)]


        #LINIOWE

        self.wektor_temperatur = [[0] for _ in range(data.NodesNumber)]
        self.dosumowany_P = []




    def fill_wektor_temperatur(self):
        for i in range(data.NodesNumber):
            self.wektor_temperatur[i][0] = self.initial_temp
           # print(self.wektor_temperatur[i]) drukowanie wektora temperatur
    def upload_coordinates(self):
        self.coordinate_list_help = []
        for o in range(data.NodesNumber):
            self.coordinate_list_help.append(self.coordinate_list[o])
            #print(self.coordinate_list_help[o])

    def upload_node_ids(self): #upload nodów oraz bc do poszczegolnych elementow
        self.bc_matrix = []
        self.node_ids_help = []
        for i in range(4):
            self.node_ids_help.append(self.node_ids)
        print(self.node_ids_help[0].ID)

        for k in range(4):
            found = False

            for j in range(len(self.bc_list)):
                if self.node_ids_help[0].ID[k] == self.bc_list[j]:
                    found = True
                    break

            if found:
                self.bc_matrix.append(1)
            else:
                self.bc_matrix.append(0)
        print("BC - warunki brzegowe:")
        print(self.bc_matrix)
    #
    # def nodes_print(self):
    #     print(self.node_ids_help[0].ID)

    def set_coordinates(self):
        # for k, node_id in enumerate(self.node_ids_help):
        for k in range(4):
            for i in range(data.NodesNumber):
                if self.node_ids_help[0].ID[k] == i+1:
                    self.coordinates_x.append(self.coordinate_list_help[i][0])
                    self.coordinates_y.append(self.coordinate_list_help[i][1])
                # print(self.coordinates_x)
                # print(self.coordinates_y)

    def initialize_ksi_eta_weights(self):   #inicjalizacja dla macierzy H
        self.ksi_values = []
        self.eta_values = []
        self.weights = []

        if self.num_nodes == 2:
            self.ksi_values = np.array([-1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3)])
            self.eta_values = np.array([-1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)])
            self.weights = np.array([1, 1])
        elif self.num_nodes == 3:
            self.ksi_values = np.array([-sqrt(3 / 5), 0, sqrt(3 / 5), -sqrt(3 / 5), 0, sqrt(3 / 5),
                                       -sqrt(3 / 5), 0, sqrt(3 / 5)])
            self.eta_values = np.array([-sqrt(3 / 5), -sqrt(3 / 5), -sqrt(3 / 5), 0, 0, 0,
                                       sqrt(3 / 5), sqrt(3 / 5), sqrt(3 / 5)])
            self.weights = np.array([5/9, 8/9, 5/9])
        elif self.num_nodes == 4:
            self.ksi_values = np.array([-0.861136, -0.339981, 0.339981, 0.861136, -0.861136, -0.339981, 0.339981, 0.861136,
                                       -0.861136, -0.339981, 0.339981, 0.861136, -0.861136, -0.339981, 0.339981, 0.861136])
            self.eta_values = np.array([-0.861136, -0.861136, -0.861136, -0.861136, -0.339981, -0.339981, -0.339981, -0.339981,
                                       0.339981, 0.339981, 0.339981, 0.339981, 0.861136, 0.861136, 0.861136, 0.861136])
            self.weights = np.array([(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36])
        else:
            raise ValueError("Nieprawidłowa liczba węzłów. Ilość węzłów musi wynosić 2, 3 lub 4")

    def initialize_ksi_eta_BC(self): #inicjalizacja dla macierzy HBC
        self.ksi_values_bc = []
        self.eta_values_bc = []
        self.weights = []

        if self.num_nodes == 2:
            self.ksi_values_bc = np.array([-1 / sqrt(3), 1 / sqrt(3), 1, 1, 1 / sqrt(3), -1 / sqrt(3), -1, -1])
            self.eta_values_bc = np.array([-1, -1, -1 / sqrt(3), 1 / sqrt(3), 1, 1, 1 / sqrt(3), -1 / sqrt(3)])
            self.weights = np.array([1, 1])
        elif self.num_nodes == 3:
            self.ksi_values_bc = np.array(
                [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5),
                 1, 1, 1,
                 np.sqrt(3 / 5), 0, -np.sqrt(3 / 5),
                 -1, -1, -1])
            self.eta_values_bc = np.array(
                [-1, -1, -1,
                 -np.sqrt(3 / 5), 0, np.sqrt(3 / 5),
                 1, 1, 1,
                 np.sqrt(3 / 5), 0, -np.sqrt(3 / 5)])
            self.weights = np.array([5 / 9, 8 / 9, 5 / 9])
        elif self.num_nodes == 4:
            self.ksi_values_bc = np.array(
                [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526,
                 1, 1, 1, 1,
                 0.8611363115940526, 0.3399810435848563, -0.3399810435848563, -0.8611363115940526,
                 -1, -1, -1, -1])
            self.eta_values_bc = np.array([
                -1, -1, -1, -1,
                -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526,
                1, 1, 1, 1,
                0.8611363115940526, 0.3399810435848563, -0.3399810435848563, -0.8611363115940526])
            self.weights = np.array([(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36])

    def initalize_ksi_eta_c(self):      #inicjalizacja dla macierzy C
        self.ksi_values_c = []
        self.eta_values_c = []
        self.weights = []

        if self.num_nodes == 2:
            self.ksi_values_c = np.array([-1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3)])
            self.eta_values_c = np.array([-1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)])
            self.weights = np.array([1, 1])
        elif self.num_nodes == 3:
            self.ksi_values_c = np.array([-sqrt(3 / 5), 0, sqrt(3 / 5), -sqrt(3 / 5), 0, sqrt(3 / 5),
                                       -sqrt(3 / 5), 0, sqrt(3 / 5)])
            self.eta_values_c = np.array([-sqrt(3 / 5), -sqrt(3 / 5), -sqrt(3 / 5), 0, 0, 0,
                                        sqrt(3 / 5), sqrt(3 / 5), sqrt(3 / 5)])
            self.weights = np.array([5 / 9, 8 / 9, 5 / 9])
        elif self.num_nodes == 4:
            self.ksi_values_c = np.array(
                [-0.861136, -0.339981, 0.339981, 0.861136, -0.861136, -0.339981, 0.339981, 0.861136,
                 -0.861136, -0.339981, 0.339981, 0.861136, -0.861136, -0.339981, 0.339981, 0.861136])
            self.eta_values_c = np.array(
                [-0.861136, -0.861136, -0.861136, -0.861136, -0.339981, -0.339981, -0.339981, -0.339981,
                 0.339981, 0.339981, 0.339981, 0.339981, 0.861136, 0.861136, 0.861136, 0.861136])
            self.weights = np.array(
                [(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36])
        else:
            raise ValueError("Nieprawidłowa liczba węzłów. Ilość węzłów musi wynosić 2, 3 lub 4")

    def calcualte_N_bc(self):        #kalkulowanie tabelki z "N" dla macierzy HBc
        self.matrix_N_bc = []
        for i in range(self.num_nodes*4):

            ksi_bc = []

            ksi_bc.append(1/4 * (1 - self.ksi_values_bc[i]) * (1 - self.eta_values_bc[i]))
            ksi_bc.append(1 / 4 * (1 + self.ksi_values_bc[i]) * (1 - self.eta_values_bc[i]))
            ksi_bc.append(1 / 4 * (1 + self.ksi_values_bc[i]) * (1 + self.eta_values_bc[i]))
            ksi_bc.append(1 / 4 * (1 - self.ksi_values_bc[i]) * (1 + self.eta_values_bc[i]))
            self.matrix_N_bc.append(ksi_bc)
          # print(self.matrix_N_bc[i])

    def calculate_N_c(self):    #kalkulowanie tabelki z "N" dla macierzy C
        self.matrix_N_c = []
        for i in range(self.num_nodes * self.num_nodes):
            ksi_n = []

            ksi_n.append(1 / 4 * (1 - self.ksi_values_c[i]) * (1 - self.eta_values_c[i]))
            ksi_n.append(1 / 4 * (1 + self.ksi_values_c[i]) * (1 - self.eta_values_c[i]))
            ksi_n.append(1 / 4 * (1 + self.ksi_values_c[i]) * (1 + self.eta_values_c[i]))
            ksi_n.append(1 / 4 * (1 - self.ksi_values_c[i]) * (1 + self.eta_values_c[i]))

            self.matrix_N_c.append(ksi_n)
            # print(self.matrix_N_c[i])
        #print(self.matrix_N_c)
    def calcualte_distance_det(self):
        self.dlugosc_scian = []
        x = self.coordinates_x
        y = self.coordinates_y

        for j in range(4):
            if j < 3:
                dx = x[j] - x[j+1]
                dy = y[j] - y[j+1]
            else:
                dx = x[j] - x[0]
                dy = y[j] - y[0]
            self.dlugosc_scian.append(sqrt(dx**2 + dy**2))
            self.det_bc.append(self.dlugosc_scian[j] / 2)


        print("X'y dla tego elementu:")
        print(x)
        print("Y'y dla tego elementu:")
        print(y)
        print("dlugosci L:")
        print(self.dlugosc_scian)
        print("det[J]:")
        print(self.det_bc)

    def calculate_matrix_HBC(self):
        self.matrix_Hbc = []
        for i in range(self.num_nodes * 4):
            matrix_helper = []
            for j in range(4):
                matrix_helper.extend(self.matrix_N_bc[i][k] * self.matrix_N_bc[i][j] for k in range(4))

            self.matrix_Hbc.append(matrix_helper)

            # print(f"Macierz HBC cząstkowe-{i + 1}")
            # print(self.matrix_Hbc[-1][:4])
            # print(self.matrix_Hbc[-1][4:8])
            # print(self.matrix_Hbc[-1][8:12])
            # print(self.matrix_Hbc[-1][12:16])
            # print("\n")
        # print("\nwszystkie macierze razem wzięte:")
        # for j in range(self.num_nodes * 4):
        #     print(self.matrix_Hbc[j])

    def calculate_matrix_P(self):  #calkowanie wektora P
        self.matrix_P_pre = []
        for i in range(self.num_nodes * 4):
            matrix_helper = []
            for j in range(4):
                matrix_helper.extend(self.matrix_N_bc[i][k] * 1200 for k in range(4))

            self.matrix_P_pre.append(matrix_helper)

            # print(f"Wektor cząstkowy P-{i + 1}")
            # print(self.matrix_P_pre[-1][:4])
            # print(self.matrix_P_pre[-1][4:8])
            # print(self.matrix_P_pre[-1][8:12])
            # print(self.matrix_P_pre[-1][12:16])
            #print("\n")
        # print("\nwszystkie macierze razem wzięte:")
        # for j in range(self.num_nodes * 4):
        #     print(self.matrix_Hbc[j])
        print("\n")
    def example_loop_hbc(self): #hbc do agregacji
        self.matrix_Hbc_sum = []

        for i in range(0, self.num_nodes * 4, self.num_nodes):
            matrix_helper = []
            for k in range(16):
                element_sum = sum(self.weights[j] * self.matrix_Hbc[i + j][k] for j in range(self.num_nodes))
                matrix_helper.append(self.alfa * element_sum)
            self.matrix_Hbc_sum.append(matrix_helper)
        print("\n")
        for l in range(4):
            for i in range(len(self.matrix_Hbc_sum[l])):
                self.matrix_Hbc_sum[l][i] *= self.det_bc[l]
            #print(self.det_bc[l])
            # print(f"\nMacierz HBC{l + 1}:")
            # for i in range(0, len(self.matrix_Hbc_sum[l]), 4):
            #     print(self.matrix_Hbc_sum[l][i:i + 4])

    def example_loop_wektor_P(self): #wektor P agregacja
        self.matrix_P_sum = []
        for i in range(0, self.num_nodes * 4, self.num_nodes):
            matrix_helper = []
            for k in range(16):
                element_sum = sum(self.weights[j] * self.matrix_P_pre[i + j][k] for j in range(self.num_nodes))
                matrix_helper.append(self.alfa * element_sum)
            self.matrix_P_sum.append(matrix_helper)
        for l in range(4):
            for i in range(len(self.matrix_P_sum[l])):
                self.matrix_P_sum[l][i] *= self.det_bc[l]

            # print(f"\nWektor P{l + 1}:")
            # print(self.matrix_P_sum[l][:4])
    def calcualte_final_matrix_hbc(self):
        self.matrix_Hbc_final = []
        value = 0.0

        for j in range(len(self.matrix_Hbc_sum[0])):
            for k in range(4):
                if k < 3:
                    if self.bc_matrix[k] and self.bc_matrix[k+1] == 1:
                        value += self.matrix_Hbc_sum[k][j]
                else:
                    if self.bc_matrix[0] and self.bc_matrix[k] == 1:
                        value += self.matrix_Hbc_sum[k][j]

            self.matrix_Hbc_final.append(value)
            value = 0.0

        # print(f"\nFinalna Macierz HBC:")
        # for i in range(0, len(self.matrix_Hbc_final), 4):
        #     print(self.matrix_Hbc_final[i:i + 4])
        # print("\n")

    def calculate_final_wektor_P(self):  #wektor P
        self.matrix_P = []
        value = 0.0

        for j in range(len(self.matrix_P_sum[0])):
            for k in range(4):
                if k < 3:
                    if self.bc_matrix[k] and self.bc_matrix[k + 1] == 1:
                        value += self.matrix_P_sum[k][j]
                else:
                    if self.bc_matrix[0] and self.bc_matrix[k] == 1:
                        value += self.matrix_P_sum[k][j]

            self.matrix_P.append(value)
            value = 0.0
        # print("\nFinalna Wektor P:")
        # print(self.matrix_P[:4])
        # print("\n")
    def calculate_dNdKsi_matrix(self):  #kalkulowanie tabelki z "N" dla macierzy H
        self.dNdKsi_matrix = []

        for i in range(self.num_nodes * self.num_nodes):
            dNdKsi_row = []

            dNdKsi_row.append(-1 / 4 * (1 - self.eta_values[i])) #dNdKsi[0]
            dNdKsi_row.append(1 / 4 * (1 - self.eta_values[i])) #(...)
            dNdKsi_row.append(1 / 4 * (1 + self.eta_values[i]))
            dNdKsi_row.append(-1 / 4 * (1 + self.eta_values[i]))

            self.dNdKsi_matrix.append(dNdKsi_row)

    def calculate_dNdEta_matrix(self):      #kalkulowanie tabelki z "N" dla macierzy H
        self.dNdEta_matrix = []

        for i in range(self.num_nodes * self.num_nodes):
            dNdEta_row = []

            dNdEta_row.append(-1 / 4 * (1 - self.ksi_values[i]))
            dNdEta_row.append(-1 / 4 * (1 + self.ksi_values[i]))
            dNdEta_row.append(1 / 4 * (1 + self.ksi_values[i]))
            dNdEta_row.append(1 / 4 * (1 - self.ksi_values[i]))

            self.dNdEta_matrix.append(dNdEta_row)

    def calculate_Jakob(self):  #Jacob
        for pc in range(self.num_nodes * self.num_nodes):
            self.matrix_jakobian_list[pc] = []

            self.dXdKsi_help = [0.0] * self.num_nodes * self.num_nodes
            self.dYdKsi_help = [0.0] * self.num_nodes * self.num_nodes
            self.dXdEta_help = [0.0] * self.num_nodes * self.num_nodes
            self.dYdEta_help = [0.0] * self.num_nodes * self.num_nodes

            for j in range(4):
                self.dXdKsi_help[pc] += (
                    self.dNdKsi_matrix[pc][j] * self.coordinates_x[j])
                self.dYdKsi_help[pc] += (
                    self.dNdKsi_matrix[pc][j] * self.coordinates_y[j])
                self.dXdEta_help[pc] += (
                    self.dNdEta_matrix[pc][j] * self.coordinates_x[j])
                self.dYdEta_help[pc] += (
                    self.dNdEta_matrix[pc][j] * self.coordinates_y[j])


            self.matrix_jakobian_list[pc].append(self.dXdKsi_help[pc])
            self.matrix_jakobian_list[pc].append(self.dYdKsi_help[pc])
            self.matrix_jakobian_list[pc].append(self.dXdEta_help[pc])
            self.matrix_jakobian_list[pc].append(self.dYdEta_help[pc])
            # self.det = 0.0
            # self.det = self.matrix_jakobian_list[pc][0] * self.matrix_jakobian_list[pc][3] - \
            #            self.matrix_jakobian_list[pc][2] * self.matrix_jakobian_list[pc][1]

            # print(f"\nMacierz Jakobianu pc{pc + 1}:")
            # for i in range(0, len(self.matrix_jakobian_list[pc]), 2):
            #     print(self.matrix_jakobian_list[pc][i:i + 2])


    def flip_matrix_jakob(self): #jacob odwrocenie
        for pc in range(self.num_nodes *self.num_nodes):
            tmp = []
            tmp = self.matrix_jakobian_list[pc][0]
            self.matrix_jakobian_list[pc][0] = self.matrix_jakobian_list[pc][3]
            self.matrix_jakobian_list[pc][3] = tmp
            self.matrix_jakobian_list[pc][1] *= -1
            self.matrix_jakobian_list[pc][2] *= -1

            # print(f"\nMacierz Jakobianu po odwroceniu dla pc{pc + 1}:")
            # for i in range(0, len(self.matrix_jakobian_list[pc]), 2):
            #     print(self.matrix_jakobian_list[pc][i:i + 2])

    def calculate_jakob_matrix_det(self):  #jacob wyznacznik
        self.det = []
        for pc in range(self.num_nodes * self.num_nodes):
            self.det.append(self.matrix_jakobian_list[pc][0] * self.matrix_jakobian_list[pc][3] - self.matrix_jakobian_list[pc][2] * self.matrix_jakobian_list[pc][1])
            #print(self.det)
            var = 1.0/self.det[pc]
            for j in range(4):
                self.matrix_jakobian_list[pc][j] *= var
            # print(f"\nMacierz Jakobianu po wymnozeniu dla pc{pc + 1}:")
            # for i in range(0, len(self.matrix_jakobian_list[pc]), 2):
            #     print(self.matrix_jakobian_list[pc][i:i + 2])


    def matrix_multiplication_dx(self):
            self.matrix_help_x = []
            for pc in range(self.num_nodes * self.num_nodes):
                matrix_help_x_row = []
                for j in range(4):
                    result = self.matrix_jakobian_list[pc][0] * self.dNdKsi_matrix[pc][j] + self.matrix_jakobian_list[pc][1] * \
                             self.dNdEta_matrix[pc][j]
                    matrix_help_x_row.append(result)
                self.matrix_help_x.extend(matrix_help_x_row)
            # print("\nMacierz dNdX:")
            # for i in range(0, len(self.matrix_help_x), 4):
            #     print(self.matrix_help_x[i:i + 4])

    def matrix_multiplication_dy(self):

        self.matrix_help_y = []
        for pc in range(self.num_nodes * self.num_nodes):
            matrix_help_y_row = []
            for j in range(4):
                result = self.matrix_jakobian_list[pc][2] * self.dNdKsi_matrix[pc][j] + self.matrix_jakobian_list[pc][3] * \
                        self.dNdEta_matrix[pc][j]
                matrix_help_y_row.append(result)
            self.matrix_help_y.extend(matrix_help_y_row)
        # print("\nMacierz dNdY:")
        # for i in range(0, len(self.matrix_help_y), 4):
        #     print(self.matrix_help_y[i:i + 4])

    def matrix_multiplication_matrixT_X(self):  #mnozienie przez transponowaną macierz dla Y
        for pc in range(self.num_nodes * self.num_nodes):
            self.matrix_x_pc[pc] = []
            for i in range(pc * 4, (pc + 1) * 4):
                matrix_helper = []
                for j in range(pc * 4, (pc + 1) * 4):
                    matrix_helper.append(self.matrix_help_x[i] * self.matrix_help_x[j])
                self.matrix_x_pc[pc].extend(matrix_helper)
            # print(f"\n[X] dNdX * dNdX[T] dla pc{pc + 1}:")
            # for i in range(0, len(self.matrix_x_pc[pc]), 4):
            #     print(self.matrix_x_pc[pc][i:i + 4])

    def matrix_multiplication_matrixT_Y(self): #mnozienie przez transponowaną macierz dla X
        for pc in range(self.num_nodes * self.num_nodes):
            self.matrix_y_pc[pc] = []
            for i in range(pc * 4, (pc + 1) * 4):
                matrix_helper = []
                for j in range(pc * 4, (pc + 1) * 4):
                    matrix_helper.append(self.matrix_help_y[i] * self.matrix_help_y[j])
                self.matrix_y_pc[pc].extend(matrix_helper)
            # print(f"\n[Y] dNdY * dNdY[T] dla pc{pc + 1}:")
            # for i in range(0, len(self.matrix_y_pc[pc]), 4):
            #     print(self.matrix_y_pc[pc][i:i + 4])

    def calculate_matrix_H(self): #całkowanie macierzy H
        self.matrix_H_pc = [[] for _ in range(self.num_nodes * self.num_nodes)]  # Inicjalizacja listy dla pc
        for pc in range(self.num_nodes * self.num_nodes):
            for i in range(16):
                matrix_helper = []
                matrix_helper.append(
                    25 * (self.matrix_x_pc[pc][i] + self.matrix_y_pc[pc][i]) * self.det[pc])
                self.matrix_H_pc[pc].extend(matrix_helper)
            # print(f"\nMacierz H dla pc{pc + 1}:")
            # for i in range(0, len(self.matrix_H_pc[pc]), 4):
            #     print(self.matrix_H_pc[pc][i:i + 4])

    def calculate_matrix_H_nodes(self): #całkowanie maciery H
        self.matrix_H_final = []

        for i in range(16):
            matrix_helper = []
            value = 0.0

            for k in range(self.num_nodes):
                for l in range(self.num_nodes):
                    value += self.matrix_H_pc[k * self.num_nodes + l][i] * self.weights[k] * self.weights[l]

            matrix_helper.append(value)
            self.matrix_H_final.extend(matrix_helper)

        # print("Macierz H finalna: ")
        # for i in range(0, len(self.matrix_H_final), 4):
        #     print(self.matrix_H_final[i:i + 4])
        # print("\n")

    def calculate_matrix_C(self): #całkowanie macierzy C
        #self.matrix_C_pc
        x = 0
        for k in range(self.num_nodes):
            for l in range(self.num_nodes):
                for i in range(4):
                    for j in range(4):
                        self.matrix_C_pc[k * self.num_nodes + l][i][j] = (700 * 7800 * self.det[k * self.num_nodes + l] *
                            self.matrix_N_c[k * self.num_nodes + l][i] *
                            self.matrix_N_c[k * self.num_nodes + l][j]) * self.weights[k] * self.weights[l]

                        self.matrix_C_final[i][j] += self.matrix_C_pc[k * self.num_nodes + l][i][j]



        #drukowanie:
        # for matrix_2d in self.matrix_C_pc:
        #     print(f"macierz C dla pc{x + 1}", x)
        #     for row in matrix_2d:
        #         print(row)
        #     print()
        #     x += 1

        # self.matrix_C_final
        print("Macierz C finalna:")
        for row in self.matrix_C_final:
            print(row)


    def agregacja(self):
        x = 0.0
        y = 0.0
        counter = 0
        for i in range(4):
            for j in range(4):
                x = self.node_ids_help[0].ID[i]
                y = self.node_ids_help[0].ID[j]
                sum = 0.0
                sum2 = 0.0
                sum = self.matrix_H_final[counter]
                self.lista_2d[x-1][y-1] += sum
                sum2 = self.matrix_Hbc_final[counter]
                self.lista_2d_hbc[x-1][y-1] += sum2
                counter += 1


        return self.lista_2d, self.lista_2d_hbc

    def agregacja_p(self):
        x = 0.0
        for i in range(4):
            x = self.node_ids_help[0].ID[i]
            sum = self.matrix_P[i]
            self.lista_2d_p[x-1][0] += sum

        return self.lista_2d_p

    def agregacja_c(self):
        x = 0.0
        y = 0.0
        for i in range(4):
            for j in range(4):
                x = self.node_ids_help[0].ID[i]
                y = self.node_ids_help[0].ID[j]
                sum = 0.0
                sum = self.matrix_C_final[i][j]

                self.lista_2d_c[x-1][y-1] += sum / self.step_time

        return self.lista_2d_c

    def zmiana_temp(self, temperatura):
        self.wektor_temperatur = temperatura
    def dosumowanie_do_P(self):
        self.dosumowany_P = np.dot(self.lista_2d_c ,self.wektor_temperatur)
        # print("dosumowany_P dziwne")
        # print(self.dosumowany_P)
        # for i in range(self.num_nodes):
        #     print(self.wektor_temperatur[i])
        self.dosumowany_P += self.lista_2d_p


    def conv_H(self):
        self.matrix_H_final_2d = [self.matrix_H_final[i:i + 4] for i in range(0, len(self.matrix_H_final), 4)]
        print("Macierz H")
        #print(self.matrix_H_final_2d)
        for row in self.matrix_H_final_2d:
            print(row)

    def conv_Hbc(self):
        self.matrix_Hbc_final_2d = [self.matrix_Hbc_final[i:i + 4] for i in range(0, len(self.matrix_Hbc_final), 4)]
        print("Macierz Hbc")
        #print(self.matrix_Hbc_final_2d)
        for row in self.matrix_Hbc_final_2d:
            print(row)

    def conv_P(self):
        print("Wektor P:")
        for i in range(4):
            self.matrix_P_2d[i] = self.matrix_P[i:i+1]
            #print(self.matrix_P[i:i+1])
        print(self.matrix_P_2d)

    def method(self):

        self.fill_wektor_temperatur()
        self.upload_coordinates()
        self.upload_node_ids()
        self.set_coordinates()
        self.initialize_ksi_eta_weights()

        #metody dla macierzy H
        self.calculate_dNdKsi_matrix()
        self.calculate_dNdEta_matrix()

        self.calculate_Jakob()
        self.flip_matrix_jakob()
        self.calculate_jakob_matrix_det()

        self.matrix_multiplication_dx()
        self.matrix_multiplication_dy()
        self.matrix_multiplication_matrixT_X()
        self.matrix_multiplication_matrixT_Y()

        self.calculate_matrix_H()
        self.calculate_matrix_H_nodes()

        #metody od hbc
        self.initialize_ksi_eta_BC()
        self.calcualte_N_bc()
        self.calcualte_distance_det()
        self.calculate_matrix_HBC()
        self.example_loop_hbc()
        self.calcualte_final_matrix_hbc()

        #metody dla wektora P
        self.calculate_matrix_P()
        self.example_loop_wektor_P()
        self.calculate_final_wektor_P()


        #macierz C

        self.initalize_ksi_eta_c()
        self.calculate_N_c()
        self.calculate_matrix_C()



        #konwertowanie na 2D poniewaz zrobilem wszystko 1D
        #konwertuje macierz H, Hbc, oraz wektor P

        self.conv_H()
        self.conv_Hbc()
        self.conv_P()


        #agregacja

        self.agregacja()
        self.agregacja_p()
        self.agregacja_c()

        #liniowe
        self.dosumowanie_do_P()


class GlobalData:
    def __init__(self):
        initial_values = {
            'SimulationTime': 0,
            'SimulationStepTime': 0,
            'Conductivity': 0,
            'Alfa': 0,
            'Tot': 0,
            'InitialTemp': 0,
            'Density': 0,
            'SpecificHeat': 0,
            'NodesNumber': 0,
            'ElementsNumber': 0
        }

        for attr, value in initial_values.items():
            setattr(self, attr, value)

class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Element:
    def __init__(self, node_ids):
        self.ID = node_ids

class Grid:
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.elements = elements


class Matrix20X20:
    def __init__(self):
        self.matrix = []

def load_data(filename):
    global_data = GlobalData()

    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            key = " ".join(parts[:-1])
            value = parts[-1]
            if key == 'SimulationTime':
                global_data.SimulationTime = int(value)
            elif key == 'SimulationStepTime':
                global_data.SimulationStepTime = int(value)
            elif key == 'Conductivity':
                global_data.Conductivity = int(value)
            elif key == 'Alfa':
                global_data.Alfa = int(value)
            elif key == 'Tot':
                global_data.Tot = int(value)
            elif key == 'InitialTemp':
                global_data.InitialTemp = int(value)
            elif key == 'Density':
                global_data.Density = int(value)
            elif key == 'SpecificHeat':
                global_data.SpecificHeat = int(value)
            elif key == 'Nodes number':
                global_data.NodesNumber = int(value)
            elif key == 'Elements number':
                global_data.ElementsNumber = int(value)
    return global_data

def load_nodes(filename):
    nodes = []

    with open(filename, 'r') as file:
        inside_node_section = False

        for line in file:
            if line.strip() == '*Node':
                inside_node_section = True
                continue

            if inside_node_section:
                parts = line.split(',')
                if len(parts) == 3:
                    x = float(parts[1].strip())
                    y = float(parts[2].strip())
                    node = Node(x, y)
                    nodes.append(node)
                else:
                    inside_node_section = False  # Koniec sekcji węzłów

    return nodes

def load_elements(filename):
    elements = []
    is_element_section = False

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('*Element'):
                is_element_section = True
                continue
            elif is_element_section and line.startswith('*BC'):
                break
            elif is_element_section and line:
                parts = line.split(',')
                node_ids = [int(part.strip()) for part in parts[1:]]
                element = Element(node_ids)
                elements.append(element)
            elif not line:
                is_element_section = False  # Koniec sekcji elementów

    return elements

def load_bc(filename):
    bc_values = []

    with open(filename, 'r') as file:
        is_bc_section = False
        for line in file:
            line = line.strip()
            if line.startswith('*BC'):
                is_bc_section = True
                continue
            elif is_bc_section and line:
                bc_parts = [int(part.strip()) for part in line.split(',')]
                bc_values.extend(bc_parts)
            elif not line:
                is_bc_section = False  # Koniec sekcji BC

    return bc_values

#zmiana testowego
filename = "Test2_4_4_MixGrid.txt"

data = load_data(filename)
nodes = load_nodes(filename)
elements = load_elements(filename)
bc_values = load_bc(filename)

nodes_number = []
InitialTemp = []
StepTime = []


if data is not None:
    print("SimulationTime:", data.SimulationTime)
    print("SimulationStepTime:", data.SimulationStepTime)
    print("Conductivity:", data.Conductivity)
    print("Alfa:", data.Alfa)
    print("Tot:", data.Tot)
    print("InitialTemp:", data.InitialTemp)
    print("Density:", data.Density)
    print("SpecificHeat:", data.SpecificHeat)
    print("Nodes number:", data.NodesNumber)
    print("Elements number:", data.ElementsNumber)
    InitialTemp = data.InitialTemp
    StepTime = data.SimulationStepTime
    NodesNumber = data.NodesNumber

if nodes is not None:
    for i, node in enumerate(nodes, start=1):
        print(f"Node {i}; x{node.x}, y{node.y}")

if elements is not None:
    for i, element in enumerate(elements, start=1):
        print(f"Element {i}: ID={element.ID}")

if nodes and elements:
    grid = Grid(nodes, elements)
    print("Liczba węzłów:", len(grid.nodes))
    print("Liczba elementów:", len(grid.elements))
    nodes_number = len(grid.nodes)


if bc_values:
    print("BC Values:", bc_values)

bc_list = []

var3 = len(bc_values)
for l in range(var3):
    bc_list.append(bc_values[l])
# for l in range(bc_count):
#     bc_list.append(bc_values[l])


coordinate_list = []
var2 = len(nodes)
for j in range(var2):
    coordinate_list.append([nodes[j].x, nodes[j].y])
    #print(coordinate_list[j])
    # print(f"Koordynaty x{coordinate_list[j].x}")

element_list = []

var = len(elements)
for i in range(var):
    element_list.append(elements[i])

    #print(nodes[0].x)
    # for k in range(4):
    #     print(element_list[i].ID[k])


#print(InitialTemp)
macierz_wynikowa_h = [[0] * nodes_number for _ in range(nodes_number)]
macierz_wynikowa_hbc = [[0] * nodes_number for _ in range(nodes_number)]
macierz_wynikowa_SUMA = [[0] * nodes_number for _ in range(nodes_number)]
macierz_wynikowa_p = [[0] for _ in range(nodes_number)]
macierz_wynikowa_c = [[0] for _ in range(nodes_number)]
macierz_wynikowa_p_help = [[0] for _ in range(nodes_number)]

temperatura_pomoc = []
wektor_p_pomoc = []


lista_elementow_cala = []
for i in range(len(elements)):
    element = UniversalElement(4, element_list[i], coordinate_list, bc_list, nodes_number, InitialTemp, StepTime)
    lista_elementow_cala.append(element)
    # print(f"\nmacierz {i + 1}")
    lista_elementow_cala[i].method()
    temperatura_pomoc = element.lista_2d_c


    macierz_wynikowa_h += np.array(element.lista_2d)
    macierz_wynikowa_hbc += np.array(element.lista_2d_hbc)
    macierz_wynikowa_p += np.array(element.lista_2d_p)
    macierz_wynikowa_c += np.array(element.lista_2d_c)

        # print("macierz H w elemencie")
        # for wiersz in element.lista_2d:
        #     print(wiersz)
        # print("macierz Hbc w elemencie")
        # for wiersz in element.lista_2d_hbc:
        #     print(wiersz)
    macierz_wynikowa_SUMA += np.array(element.lista_2d)
    macierz_wynikowa_SUMA += np.array(element.lista_2d_hbc)
    macierz_wynikowa_SUMA += np.array(element.lista_2d_c)
    macierz_wynikowa_p_help += np.array(element.dosumowany_P)

max_line_width = nodes_number

#DRUKOWANIE=== AGREGACJI

# print("Agregacja dla macierzy H:")
# print(np.array_str(macierz_wynikowa_h, max_line_width=float('inf')))
#
# print("Agregacja dla macierzy Hbc:")
# print(np.array_str(macierz_wynikowa_hbc, max_line_width=float('inf')))

print("Agregacja dla macierzy C/timestep:")
print(np.array_str(macierz_wynikowa_c, max_line_width=float('inf')))

print("aggregacja wektor P:")
print(macierz_wynikowa_p)

print("Aggregacja macierz zsumowana Hbc oraz H oraz C")
print(np.array_str(macierz_wynikowa_SUMA, max_line_width=float('inf')))
#
# print("zsumowana ")
# print(macierz_wynikowa_p_help)


#temp = np.dot(t, macierz_wynikowa_c)
# k = np.linalg.solve(macierz_wynikowa_SUMA, temp)
# print(k)

#testy
# print("rezultat p:")
# print(macierz_wynikowa_p)
#
# print("rezultat c")
# print(macierz_wynikowa_c)


t = np.linalg.solve(macierz_wynikowa_SUMA, macierz_wynikowa_p_help)
print("Tempratura po czasie 50 sekund")
print("MIN:", min(t), "MAX:", max(t))


t0 = np.full(NodesNumber, InitialTemp)
#print(t0)

k = 2

for i in range(data.SimulationStepTime, data.SimulationTime, data.SimulationStepTime):
    x = np.dot(macierz_wynikowa_c, t)
    x += macierz_wynikowa_p
    t = np.linalg.solve(macierz_wynikowa_SUMA, x)
    sek = 50
    print("Temperatura po czasie:",sek*k)
    print("MIN:", min(t), "MAX:", max(t))
    k += 1




# print("=======================do zadania:")
# lista_help = [1, 2, 3, 4]
# bc_list_help = [1, 1, 1, 1]
# obietk_help = UniversalElement(2, lista_help, coordinate_list, bc_list_help, 4, 100)
# obietk_help.method()
# print("\ndNdKsi:")
# for row in obiekt.dNdKsi_matrix:
#     print(row)
#
# print("\ndNdEta:")
# for row in obiekt.dNdEta_matrix:
#     print(row)
#
# print("\n N1..N4xN1...N4:")
# for row in obiekt.matrix_N_bc:
#     print(row)

