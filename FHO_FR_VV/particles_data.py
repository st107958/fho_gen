from FHO_FR_VV.levels import levels_e_ex
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
N2.we = [235857]                           # DB, 1/m
N2.wexe = [1432]
N2.weye = [-0.226]

e_i = levels_e_ex(N2, 0)
N2.add_ev_i(e_i-e_i[0])

# print(N2.mass)
# print(e_i)
# print(N2.ev_i)

# print(CO.ev_i)


# m1 = CO
# m2 = CO
#
# m_red = m1.mass * m2.mass / (m1.mass + m2.mass)
# print(m_red)

# print(N2 == N2)
#
# m1 = N2
# m2 = N2
#
# print(m1 == m2)


class VSSData:
    def __init__(self, dref, Tref, omega):
        self.dref = dref  # Reference diameter, m
        self.Tref = Tref  # Reference temperature, K
        self.omega = omega  # Viscosity temperature exponent

    def __repr__(self):
        return f"VSS_data(dref={self.dref:.3e} m, Tref={self.Tref} K, omega={self.omega})"


# База данных параметров VSS
VSS_DATABASE = {
    "N2-N2": VSSData(dref=4.04e-10, Tref=273, omega=0.686),
    "N2-O2": VSSData(dref=3.604e-10, Tref=273, omega=0.703),
    "N2-NO": VSSData(dref=4.391e-10, Tref=273, omega=0.756),
    "N2-N": VSSData(dref=4.088e-10, Tref=273, omega=0.762),
    "N2-O": VSSData(dref=3.222e-10, Tref=273, omega=0.702),
    "O2-N2": VSSData(dref=3.604e-10, Tref=273, omega=0.703),
    "O2-O2": VSSData(dref=3.896e-10, Tref=273, omega=0.7),
    "O2-NO": VSSData(dref=4.054e-10, Tref=273, omega=0.718),
    "O2-N": VSSData(dref=3.721e-10, Tref=273, omega=0.757),
    "O2-O": VSSData(dref=3.734e-10, Tref=273, omega=0.76),
    "NO-N2": VSSData(dref=4.391e-10, Tref=273, omega=0.756),
    "NO-O2": VSSData(dref=4.054e-10, Tref=273, omega=0.718),
    "NO-NO": VSSData(dref=4.218e-10, Tref=273, omega=0.737),
    "NO-N": VSSData(dref=4.028e-10, Tref=273, omega=0.788),
    "NO-O": VSSData(dref=3.693e-10, Tref=273, omega=0.752)
}


# Функция для получения параметров
def get_vss_data(collision_pair):
    """Get VSS parameters for collision pair"""
    key = f"{collision_pair.particle1.name}-{collision_pair.particle2.name}"
    return VSS_DATABASE.get(key)