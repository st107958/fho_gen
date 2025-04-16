import numpy as np
from scipy.integrate import quad
from scipy.special import factorial

from p_vv_mm import *
from constants import *


# def integrand_vv_mm(m1, m2, i1, f1, i2, f2, E, T_inv_cm):
#     # symm_E = E + delta_E/2 + np.power(delta_E, 2)/16/E
#     return (np.power((E / T_inv_cm), 3) * p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez')
#             * (np.exp(- (E / T_inv_cm))))


# def integrand(x, m1, m2, i1, f1, i2, f2, T_inv_cm):
#     E = x * T_inv_cm
#     return np.power(x, 3) * p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez') * np.exp(-x)

def integrand(x, *args):
    m1, m2, i1, f1, i2, f2, T_inv_cm = args
    E = x * T_inv_cm
    p_val = p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez')
    p_val = np.nan_to_num(p_val, nan=0.0)
    return x**3 * p_val * np.exp(-x)

def k_vv_mm(m1, m2, i1, f1, i2, f2, T):
    T_K = T
    T_inv_cm = T * k / h / c / 100  # 1/cm

    m_red = (m1.mass * m2.mass) / (m1.mass + m2.mass)  # приведенная масса, kg
    r = (m1.diameter + m2.diameter) / 2  # collision diameter

    mean_u = np.sqrt(8 * k * T_K / (pi * m_red))  # m/s

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
        args = (m1, m2, i1, f1, i2, f2, T_inv_cm)
        integral, error = quad(integrand, 0, np.inf, args=args)

        # integral, error = quad(integrand, 0, np.inf, args=(m1, m2, i1, f1, i2, f2, T_inv_cm))

        k_vv = np.pi * (r ** 2) * mean_u * integral

        return k_vv




print(1e6 * k_vv_mm(N2, N2, 41, 40, 40, 41, 3000))  # cm^3 / s

print(1e6 * k_vv_mm(N2, N2, 1, 3, 5, 3, 3000))  # cm^3 / s

print(1e6 * k_vv_mm(N2, N2, 0, 3, 3, 1, 3000))  # cm^3 / s
print(1e6 * k_vv_mm(N2, N2, 0, 3, 3, 0, 3000))  # cm^3 / s

