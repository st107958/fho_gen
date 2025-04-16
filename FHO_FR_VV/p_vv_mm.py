from scipy.integrate import nquad, trapezoid
from scipy.integrate import simpson
import numpy as np
import cupy as cp

from constants import *
from p_vv_mm_ij import p_vv, g, gamma
from particles_data import *


def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapez'):
    # max_angle = pi  the upper limit of theta_v1 and phi1 angles integration

    e_in_J = h * c * 100
    E = E * e_in_J

    if m1 == m2 and i1 == f2 and i2 == f1:
        result = 0

    if method == 'simps':
        maxdiv = 12  # кол-во делений по координате: 18 (больше - много памяти)

        # пределы интегрирования
        eps1 = np.linspace(0, 1, maxdiv)
        eps2 = np.linspace(0, 1, maxdiv)
        y = np.linspace(0, 1, maxdiv)
        v1 = np.linspace(-pi / 2, pi / 2, maxdiv)
        phi1 = np.linspace(-pi / 2, pi / 2, maxdiv)
        v2 = np.linspace(-pi / 2, pi / 2, maxdiv)
        phi2 = np.linspace(-pi / 2, pi / 2, maxdiv)

        # сеткa
        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2, indexing='ij')  # indexing='ij'

        # вычисление значений функций на сетке
        F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1, EPS2, Y, V1, PHI1, V2, PHI2)

        integral = simpson(y=F, x=eps1, axis=0)
        integral = simpson(y=integral, x=eps2, axis=0)
        integral = simpson(y=integral, x=y, axis=0)
        integral = simpson(y=integral, x=v1, axis=0)
        integral = simpson(y=integral, x=phi1, axis=0)
        integral = simpson(y=integral, x=v2, axis=0)
        integral = simpson(y=integral, x=phi2, axis=0)
        result = integral


    if method == 'trapez':
        maxdiv = 5  # кол-во делений по координате: 18 (больше - много памяти)
        result = 0

        # пределы интегрирования
        eps1 = np.linspace(0, 1, maxdiv)
        eps2 = np.linspace(0, 1, maxdiv)
        y = np.linspace(0, 1, maxdiv)
        v1 = np.linspace(-pi/2, pi/2, maxdiv)
        phi1 = np.linspace(-pi/2, pi/2, maxdiv)
        v2 = np.linspace(-pi/2, pi/2, maxdiv)
        phi2 = np.linspace(-pi/2, pi/2, maxdiv)

        # сеткa
        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2) # indexing='ij'

        mask = ((EPS1 + EPS2) <= (1-np.power(Y, 2))/(2*(1-np.power(Y, 2)/2)))
        # mask = EPS1 + EPS2 <= 1

        # bad_points = np.where((1 - EPS1 - EPS2) * (1 - Y) < 0)
        # if len(bad_points[0]) > 0:
        #     print("Ошибка в точках:", bad_points)
        #     print("Значения EPS1, EPS2, Y:", EPS1[bad_points], EPS2[bad_points], Y[bad_points])

        # вычисление значений функций на сетке
        F = np.nan_to_num(p_vv(m1, m2, i1, f1, i2, f2, E, EPS1, EPS2, Y, V1, PHI1, V2, PHI2) * mask, nan=0.0)

        result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
            F, eps1, axis=6), eps2, axis=5), y, axis=4), v1, axis=3),
            phi1, axis=2), v2, axis=1), phi2, axis=0)

    return result

# obratnoe_m_in_J = 1.98e-23

print(p_vv_int(CO, CO, 2, 0, 1, 3, 10000, 'trapez'))



# print(p_vv_int(N2, N2, 1, 0, 2, 3, 10000*e_in_J, 'trapez'))
#
# print(p_vv_int(N2, N2, 4, 2, 0, 3, 10000*e_in_J, 'trapez'))
#
# print(p_vv_int(N2, N2, 5, 2, 0, 3, 10000*e_in_J, 'trapez'))
# print(p_vv_int(N2, N2, 5, 1, 0, 4, 10000*e_in_J, 'trapez'))
# print(p_vv_int(N2, N2, 6, 1, 0, 5, 10000*e_in_J, 'trapez'))


# print(p_vv_int(CO, CO, 2, 0, 1, 3, 100000, 'simps'))
# print(p_vv_int(N2, N2, 1, 0, 2, 3, 100000, 'simps'))
#
# print(p_vv_int(N2, N2, 4, 2, 0, 3, 100000, 'simps'))
#
# print(p_vv_int(N2, N2, 5, 2, 0, 3, 100000, 'simps'))
# print(p_vv_int(N2, N2, 5, 1, 0, 4, 100000, 'simps'))
# print(p_vv_int(N2, N2, 6, 1, 0, 5, 100000, 'simps'))


