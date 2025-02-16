from scipy.integrate import nquad, trapz
import numpy as np

from constants import *
from p_vv_mm_ij import p_vv

def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapezoidal'):


    max_angle = pi  #the upper limit of theta_v1 and phi1 angles integration

    if method == 'quad':
        error_v = 1e-2 # integration error
        limits = [(0, 1), (0, 1), (0, 1), (0, pi), (0, pi), (0, pi), (0, pi)]
        result, error = nquad(lambda eps1, eps2, y, v1, phi1, v2, phi2:
                              p_vv(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2), limits)

    if method == 'trapezoidal':
        # Параметры интегрирования
        eps1 = np.linspace(0, 1, 100)  # 100 точек для eps
        eps2 = np.linspace(0, 1, 100)
        y = np.linspace(0, 1, 100)
        v1 = np.linspace(0, pi, 100)
        phi1 = np.linspace(0, pi, 100)
        v2 = np.linspace(0, pi, 100)
        phi2 = np.linspace(0, pi, 100)

        # Создаём сетку
        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2, indexing='ij')

        # Вычисляем значения функции на сетке
        F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1, EPS2, Y, V1, PHI1, V2, PHI2)

        # Многомерное интегрирование методом трапеций
        result = trapz(trapz(trapz(trapz(trapz(trapz(trapz(
            F, phi2, axis=6), v2, axis=5), phi1, axis=4), v1, axis=3),
            y, axis=2), eps2, axis=1), eps1, axis=0)


