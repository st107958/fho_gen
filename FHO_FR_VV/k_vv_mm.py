
import numpy as np

from scipy.integrate import quad

from scipy.special import factorial, gamma
from FHO_FR_VV.p_vv_mm import *
from FHO_FR_VV.constants import *

# def get_crossection_VSS(velocity: float, collision_reduced_mass: float, VSS_data: VSSData):
#     gref = np.sqrt(2 * k * VSS_data.Tref / collision_reduced_mass) # Reference velocity, m/s
#     return np.pi * VSS_data.dref**2 * (velocity / gref) ** (1 - 2 * VSS_data.omega) / gamma(2.5 - VSS_data.omega)

# def integrand(x_array, *args):
#     m1, m2, i1, f1, i2, f2, T_inv_cm = args
#
#     # print(x_array.shape)
#
#     x_array = np.atleast_1d(x_array)
#     results = np.zeros_like(x_array)
#     for idx, x in enumerate(x_array):
#         E = x * T_inv_cm
#         results[idx] = p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez')
#     return results * np.exp(-x_array) * np.power(x_array, 3)


def integrand(x, *args):
    m1, m2, i1, f1, i2, f2, T_inv_cm = args

    E = x * T_inv_cm
    pvv_result = p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez')

    return pvv_result * np.exp(-x) * np.power(x, 3)


# def integrand1(x_array, *args):
#     m1, m2, i1, f1, i2, f2, T_inv_cm = args
#
#     m = (m1.mass * m2.mass) / (m1.mass + m2.mass)
#     x_array = np.atleast_1d(x_array)
#     probability = np.zeros_like(x_array)
#     crossection = np.zeros_like(x_array)
#
#     for idx1, x1 in enumerate(x_array):
#         E = np.power(x1, 2) * m / 2
#         probability[idx1] = p_vv_int(m1, m2, i1, f1, i2, f2, E, 'trapez')
#
#     for idx2, x2 in enumerate(x_array):
#         crossection[idx2] = get_crossection_VSS(x2, m)
#
#     return crossection * probability * np.exp(-np.power(x_array, 2) * m / 2 * k * T_inv_cm) * x_array


def k_vv_mm(m1, m2, i1, f1, i2, f2, T):
    T_K = T
    T_inv_cm = T_K * k / h / c / 100  # 1/cm

    m_red = (m1.mass * m2.mass) / (m1.mass + m2.mass)  # приведенная масса, kg
    r = (m1.diameter + m2.diameter) / 2  # collision diameter
    print(r)

    mean_u = np.sqrt(8 * k * T_K / (np.pi * m_red))  # m/s


    if m1 == m2 and i1 == f2 and i2 == f1:

        s = np.absolute(i2 - f2)
        el_lvl = 1 - 1  # electronic level

        e1_1 = m1.ev_i[el_lvl][i1]  # initial state, J
        e1_2 = m1.ev_i[el_lvl][f1]  # final state, J

        omega = np.absolute(e1_1 - e1_2) / (s * h_red)  # J/(J*s)=1/s

        ns1 = np.power((factorial(max(i1, f1)) / factorial(min(i1, f1))), (1 / s))
        ns2 = np.power((factorial(max(i2, f2)) / factorial(min(i2, f2))), (1 / s))

        z = 3 * pi * np.power(r, 2) * mean_u # m^3/s
        # z = 3 * np.power(r, 2) * mean_u  # m^3/s

        f = ((np.power((1 + (1 / np.power(2, s-1))), 4) * factorial(s+3))
             / (np.power(2, s+8) * np.power(s+1, 2) * factorial(3)))

        f_2 = np.power(alpha / omega, 2) * k * T_K / (2 * m_red)

        result = (z * f * np.power(ns1 * ns2 * f_2, s) / np.power(factorial(s), 2)
                  / np.power((1 + ((2 * ns1 * ns2 * f * f_2) / (s + 1))), s+4))

        return result*1e6 # sm^3/s

    else:
        args = (m1, m2, i1, f1, i2, f2, T_inv_cm)
        args1 = (m1, m2, i1, f1, i2, f2, T_K)
        result, error = quad(integrand, 0, np.inf, args=args, epsabs=1e-15, limit=1000)
        #result1, error1 = quad(integrand1, 0, np.inf, args=args1, epsabs=1e-15, limit=1000)
        # print(error)

        # result1 = quad(integrand, 0, 1e10, args=args)
        # result2 = quad(integrand, 1e6, np.inf, args=args)[0]
        # result = result1 + result2

        k_vv = np.pi * (r ** 2) * mean_u * result * 1e6  # sm^3/s
        #k_vv = np.pi * (r ** 2) * mean_u * result * 1e6   # sm^3/s

        return k_vv




# print(1e6 * k_vv_mm(N2, N2, 41, 40, 40, 41, 3000))  # cm^3 / s
#
# print(1e6 * k_vv_mm(N2, N2, 1, 3, 5, 3, 3000))  # cm^3 / s
#
# print(1e6 * k_vv_mm(N2, N2, 1, 3, 3, 1, 3000))  # cm^3 / s
# print(1e6 * k_vv_mm(N2, N2, 1, 4, 3, 0, 3000))  # cm^3 / s

