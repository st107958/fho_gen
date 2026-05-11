from levels import levels_e_ex
import numpy as np

class particle:
    def __init__(self):
        self.mass = 0
        self.diameter = 0
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
CO.diameter = 3.65e-10                     # m
CO.we = [216981.358, 174341, 151824]       # DB, 1/m
CO.wexe = [1328.831, 1436, 1940]
CO.weye = [1.0511, -4.5, 76.6]

e_i = levels_e_ex(CO, 0)
CO.add_ev_i(e_i-e_i[0])

e_i = levels_e_ex(CO, 1)
CO.add_ev_i(e_i-e_i[0])

e_i = levels_e_ex(CO, 2)
CO.add_ev_i(e_i-e_i[0])

N2 = particle()
N2.num_elex_levels = 1                     # number of electronical levels
N2.num_vibr_levels = [47]                  # number of vibrational levels
N2.mass = 4.651236272601599e-26            # molecular mass, kg
N2.diameter = 3.4039e-10                   # m
N2.we = [235960]                           # DB, 1/m
N2.wexe = [1445.6]
N2.weye = [-0.226]

e_i = levels_e_ex(N2, 0)
N2.add_ev_i(e_i-e_i[0])


O2 = particle()
O2.num_elex_levels = 1                     # number of electronic levels
O2.num_vibr_levels = [42]                  # number of vibrational levels (обычно ~40–45)
O2.mass = 5.313525e-26                     # molecular mass, kg
O2.diameter = 3.46e-10                     # m

# Vibrational constants for O2 (ground electronic state)
O2.we = [158000]                           # 1/m  (≈ 1580 cm^-1)
O2.wexe = [1198]                           # 1/m  (≈ 11.98 cm^-1)
O2.weye = [-3.0]                           # 1/m  (малый третий член, часто отрицательный)

e_i_O2 = levels_e_ex(O2, 0)
O2.add_ev_i(e_i_O2 - e_i_O2[0])
