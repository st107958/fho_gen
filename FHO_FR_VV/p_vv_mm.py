from scipy.integrate import nquad, trapezoid
from scipy.integrate import simpson
import numpy as np
np.set_printoptions(threshold=np.inf)
import cupy as cp

from constants import *
from p_vv_mm_ij import p_vv, g, gamma
from particles_data import *


def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapez'):
    e_in_J = h * c * 100
    E = E * e_in_J   # 1/cm --> J

    # print('E', E)

    if m1 == m2 and i1 == f2 and i2 == f1:
        raise ValueError("ksi = 0, resonance process")

    if method == 'trapez':
        maxdiv = 5  # макс. кол-во делений по координате: 18 (больше - много памяти)

        # пределы интегрирования
        eps1 = np.linspace(0, 1, maxdiv)
        eps2 = np.linspace(0, 1, maxdiv)
        y = np.linspace(0, 1, maxdiv)
        v1 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        phi1 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        v2 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        phi2 = np.linspace(-np.pi/2, np.pi/2, maxdiv)

        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2, indexing='ij')

        mask = (EPS1 + EPS2) <= (1-np.power(Y, 2))/(2*(1-np.power(Y, 2)/2))
        # mask = (EPS1 + EPS2) <= -1
        #############################################

        EPS1_filtered = EPS1[mask]
        EPS2_filtered = EPS2[mask]
        Y_filtered = Y[mask]
        V1_filtered = V1[mask]
        PHI1_filtered = PHI1[mask]
        V2_filtered = V2[mask]
        PHI2_filtered = PHI2[mask]

        F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1_filtered, EPS2_filtered, Y_filtered,
                 V1_filtered, PHI1_filtered, V2_filtered, PHI2_filtered)

        # print(F.shape)

        F_ = np.zeros_like(EPS1)  # Исходная форма (maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv)
        F_[mask] = F  # допустимые точки
        F_ = np.nan_to_num(F_)

        # print(F_.shape)

        result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
            F_, eps1, axis=6), eps2, axis=5), y, axis=4), v1, axis=3),
            phi1, axis=2), v2, axis=1), phi2, axis=0)

        # print(result)

        #############################################################

        # F = np.nan_to_num(p_vv(m1, m2, i1, f1, i2, f2, E, EPS1, EPS2, Y, V1, PHI1, V2, PHI2) * mask, nan=0.0)
        #
        # result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
        #     F, eps1, axis=6), eps2, axis=5), y, axis=4), v1, axis=3),
        #     phi1, axis=2), v2, axis=1), phi2, axis=0)
        #
        # print(result)

        result = result / pi ** 4

    return result

print(p_vv_int(CO, CO, 2, 0, 1, 3, 10000, 'trapez'))



print(p_vv_int(N2, N2, 1, 0, 2, 3, 10000, 'trapez'))
#
# print(p_vv_int(N2, N2, 4, 2, 0, 3, 10000, 'trapez'))
#
# print(p_vv_int(N2, N2, 5, 2, 0, 3, 10000, 'trapez'))
# print(p_vv_int(N2, N2, 5, 1, 0, 4, 10000, 'trapez'))
# print(p_vv_int(N2, N2, 6, 1, 0, 5, 10000, 'trapez'))



