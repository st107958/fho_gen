from scipy.integrate import trapezoid
import numpy as np


from constants import h, c
from p_vv_mm_ij import p_vv


def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapez'):
    e_in_J = h * c * 100
    E = E * e_in_J   # 1/cm --> J


    if method == 'trapez':
        maxdiv = 9  # макс. кол-во делений: 18

        # пределы интегрирования
        eps1 = np.linspace(0, 1, maxdiv)
        eps2 = np.linspace(0, 1, maxdiv)
        y = np.linspace(0, 1, maxdiv)
        v1 = np.linspace(0, np.pi, maxdiv)
        phi1 = np.linspace(0, np.pi, maxdiv)
        v2 = np.linspace(0, np.pi, maxdiv)
        phi2 = np.linspace(0, np.pi, maxdiv)

        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2, indexing='ij')

        mask = (EPS1 + EPS2) <= 1/2
        # eps1, eps2, y уже сетки
        # mask = (EPS1 >= 0) & (EPS2 >= 0) & (EPS2 <= 1 - EPS1)  # треугольник
        # # добавляем ограничение (5) из статьи
        # with np.errstate(divide='ignore', invalid='ignore'):
        #     limit = 0.5 * (1 - Y**2) / (1 - Y**2/2)
        #     # в точках y=1 знаменатель 1 - 1/2 = 0.5 → limit=0, ок
        #     # для y близких к sqrt(2) будет особенность, но y ∈ [0,1]
        # mask = mask & (EPS1 + EPS2 <= limit)

        EPS1_filtered = EPS1[mask]
        EPS2_filtered = EPS2[mask]
        Y_filtered = Y[mask]
        V1_filtered = V1[mask]
        PHI1_filtered = PHI1[mask]
        V2_filtered = V2[mask]
        PHI2_filtered = PHI2[mask]

        F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1_filtered, EPS2_filtered, Y_filtered,
                 V1_filtered, PHI1_filtered, V2_filtered, PHI2_filtered)

        F_ = np.zeros_like(EPS1)  # Исходная форма (maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv)
        F_[mask] = F  # допустимые точки
        F_ = np.nan_to_num(F_)

        result = trapezoid(
            trapezoid(
                trapezoid(
                    trapezoid(
                        trapezoid(
                            trapezoid(
                                trapezoid(F_, eps1, axis=6),
                                eps2, axis=5),
                            y, axis=4),
                        v1, axis=3),
                    phi1, axis=2),
                v2, axis=1),
            phi2, axis=0)

        result = result / (np.pi ** 4)

    return result




