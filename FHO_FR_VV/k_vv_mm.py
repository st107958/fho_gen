import numpy as np
from scipy.integrate import quad
from scipy.special import factorial
from scipy.integrate import quadrature
from scipy.integrate import romberg
import quadpy

from p_vv_mm import *
from constants import *

def integrand(x_array, *args):
    m1, m2, i1, f1, i2, f2, T_inv_cm = args

    # print(x_array.shape)


    results = np.zeros_like(x_array[0])
    for idx, x in enumerate(x_array[0]):
        E = x * T_inv_cm
        results[idx] = np.nan_to_num(p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez'), nan=0.0)
    return np.power(x_array, 3) * results * np.exp(-x_array)

# def integrand(E_array, m1, m2, i1, f1, i2, f2, T_inv_cm):
#     # m1, m2, i1, f1, i2, f2, T_inv_cm = args
#
#     results = np.zeros_like(E_array[0])
#     for idx, E in enumerate(E_array[0]):
#         results[idx] = np.nan_to_num(p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez'), nan=0.0)
#     return np.power(E_array/T_inv_cm, 3) * results * np.exp(-E_array/T_inv_cm)

# def E_coll(E, *args):
#     m1, m2, i1, f1, i2, f2, T_inv_cm = args
#
#     el_lvl = 1 - 1
#
#     e1_1 = m1.ev_i[el_lvl][i1]  # initial state, J
#     e1_2 = m1.ev_i[el_lvl][f1]  # final state, J
#     e2_1 = m2.ev_i[el_lvl][i2]  # initial state, J
#     e2_2 = m2.ev_i[el_lvl][f2]  # final state, J
#
#     delta_E = e1_1 - e2_2 - e1_2 + e2_1
#
#     return E + delta_E / 2 + (delta_E * delta_E) / (16*E)


# def integrand(E, *args):
#     m1, m2, i1, f1, i2, f2, T_inv_cm = args
#
#     E_ = E_coll(E, *args)
#
#     p_val = p_vv_int(m1, m2, i1, f1, i2, f2, E_, 'trapez')
#     p_val = np.nan_to_num(p_val, nan=0.0)
#     return np.power((E_/T_inv_cm), 3) * p_val * np.exp(-(E_/T_inv_cm))

def k_vv_mm(m1, m2, i1, f1, i2, f2, T):
    T_K = T
    T_inv_cm = T_K * k / h / c / 100  # 1/cm

    m_red = (m1.mass * m2.mass) / (m1.mass + m2.mass)  # приведенная масса, kg
    r = (m1.diameter + m2.diameter) / 2  # collision diameter

    mean_u = np.sqrt(8 * k * T_K / (np.pi * m_red))  # m/s

    if m1 == m2 and i1 == f2 and i2 == f1:

        s = np.absolute(i2 - f2)
        el_lvl = 1 - 1  # electronic level

        e1_1 = m1.ev_i[el_lvl][i1]  # initial state, J
        e1_2 = m1.ev_i[el_lvl][f1]  # final state, J

        omega = np.absolute(e1_1 - e1_2) / (s * h_red)  # J/(J*s)=1/s

        ns1 = np.power((factorial(max(i1, f1)) / factorial(min(i1, f1))), (1 / s))
        ns2 = np.power((factorial(max(i2, f2)) / factorial(min(i2, f2))), (1 / s))

        z = 3 * pi * np.power(r, 2) * mean_u

        f = ((np.power((1 + (1 / np.power(2, s-1))), 4) * factorial(s+3))
             / (np.power(2, s+8) * np.power(s+1, 2) * factorial(3)))

        f_2 = np.power(alpha / omega, 2) * k * T_K / (2 * m_red)

        result = (z * f * np.power(ns1 * ns2 * f_2, s) / np.power(factorial(s), 2)
                  / np.power((1 + ((2 * ns1 * ns2 * f * f_2) / (s + 1))), s+4))

        return result

    else:
        # args = (m1, m2, i1, f1, i2, f2, T_inv_cm)
        # result, error = quad(integrand, 0, np.inf, args=args, epsabs=1e-7, epsrel=1e-7, limit=100)
        # print(error)


        # E_max = find_adaptive_upper_limit(integrand, m1, m2, i1, f1, i2, f2, T_inv_cm)
        #
        # scheme = quadpy.c1.gauss_kronrod(100)
        # result = scheme.integrate(
        #     lambda E: integrand(E, m1, m2, i1, f1, i2, f2, T_inv_cm),
        #     [0, E_max]  # Пределы интегрирования (верхний предел вместо бесконечности)
        # )

        scheme = quadpy.e1r.gauss_laguerre(300, alpha=1.5)
        result = scheme.integrate(lambda E: integrand(E, m1, m2, i1, f1, i2, f2, T_inv_cm))

        # scheme = quadpy.c1.gauss_kronrod(100)
        # result = scheme.integrate(
        #     lambda E: integrand(E, m1, m2, i1, f1, i2, f2, T_inv_cm),
        #     [0, 120000000]  # Пределы интегрирования (верхний предел вместо бесконечности)
        # )

        # result1 = quad(integrand, 0, 1e10, args=args)
        # result2 = quad(integrand, 1e6, np.inf, args=args)[0]
        # result = result1 + result2

        # result = quadrature(integrand, 0, 1e5, args=args)[0]

        # maxdiv = 2560
        # x = np.linspace(0, np.inf, maxdiv)
        # F = integrand(x, m1, m2, i1, f1, i2, f2, T_inv_cm)
        # result = trapezoid(F, x, axis=0)


        # k_vv = np.pi * (r ** 2) * mean_u * result * 1e6 #sm^3/s
        k_vv = np.pi * (r ** 2) * mean_u * result * 1e6  # sm^3/s

        return k_vv




# print(1e6 * k_vv_mm(N2, N2, 41, 40, 40, 41, 3000))  # cm^3 / s
#
# print(1e6 * k_vv_mm(N2, N2, 1, 3, 5, 3, 3000))  # cm^3 / s
#
# print(1e6 * k_vv_mm(N2, N2, 1, 3, 3, 1, 3000))  # cm^3 / s
# print(1e6 * k_vv_mm(N2, N2, 1, 4, 3, 0, 3000))  # cm^3 / s

