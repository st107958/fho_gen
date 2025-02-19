from scipy.integrate import nquad, trapezoid
import numpy as np
import cupy as cp

from constants import *
from p_vv_mm_ij import p_vv
from particles_data import *


def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapez'):
    # max_angle = pi  the upper limit of theta_v1 and phi1 angles integration
    result = 0

    if method == 'trapez':
        maxdiv = 5  # кол-во делений по координате: 18 (больше - много памяти)

        # пределы интегрирования
        eps1 = np.linspace(0, 1, maxdiv)
        eps2 = np.linspace(0, 1, maxdiv)
        y = np.linspace(0, 1, maxdiv)
        v1 = np.linspace(0, pi, maxdiv)
        phi1 = np.linspace(0, pi, maxdiv)
        v2 = np.linspace(0, pi, maxdiv)
        phi2 = np.linspace(0, pi, maxdiv)
        # v1 = np.linspace(-pi/2, pi/2, maxdiv)
        # phi1 = np.linspace(-pi/2, pi/2, maxdiv)
        # v2 = np.linspace(-pi/2, pi/2, maxdiv)
        # phi2 = np.linspace(-pi/2, pi/2, maxdiv)

        # сеткa
        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2,
                                                        indexing='ij')

        mask = ((EPS1 + EPS2) <= 1)

        # вычисление значений функций на сетке
        F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1, EPS2, Y, V1, PHI1, V2, PHI2) * mask

        result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
            F, phi2, axis=6), v2, axis=5), phi1, axis=4), v1, axis=3),
            y, axis=2), eps2, axis=1), eps1, axis=0)

    return result

print(p_vv_int(CO, CO, 1, 0, 0, 1, 1000, 'trapez'))

