from levels import levels_e_ex
import numpy as np

class particle:
    def __init__(self):
        self.mass = 0
        self.num_elex_levels = 1
        self.num_vibr_levels = []
        self.we = []
        self.wexe = []
        self.weye = []
        self.ev_i = []

    def add_ev_i(self, arr):
        self.ev_i.append(arr)


CO = particle()
CO.num_elex_levels = 3                     # number of electronical levels
CO.num_vibr_levels = [68, 34, 18]          # number of vibrational levels
CO.mass = 4.651236272601599e-26            # molecular mass, kg
CO.we = [216981.358, 174341, 151824]       # DB, м-1
CO.wexe = [1328.831, 1436, 1940]
CO.weye = [1.0511, -4.5, 76.6]

e_i = levels_e_ex(CO, 0)
CO.add_ev_i(e_i-e_i[0])

e_i = levels_e_ex(CO, 1)
CO.add_ev_i(e_i-e_i[0])

e_i = levels_e_ex(CO, 2)
CO.add_ev_i(e_i-e_i[0])

# print(CO.ev_i)


# m1 = CO
# m2 = CO
#
# m_red = m1.mass * m2.mass / (m1.mass + m2.mass)
# print(m_red)

