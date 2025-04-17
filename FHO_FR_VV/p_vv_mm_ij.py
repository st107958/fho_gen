from constants import *

import numpy as np
from scipy.special import factorial


# def gamma(eps1, eps2, y, v1, phi1, v2, phi2):
#     y_squared = y*y
#     threshold = (1 - y_squared) / (2 * (1 - y_squared / 2))
#     sum_eps = eps1 + eps2
#
#     condition = ((sum_eps < threshold) & (1 - y > 0))
#     gamma_values = np.where(condition, np.maximum(0, -0.5 * np.sin(2 * v1) * np.cos(phi1) * np.sqrt(eps1)
#                                                   - 0.5 * np.sin(2 * v2) * np.cos(phi2) * np.sqrt(eps2)
#                                                   + np.sqrt((1 - eps1 - eps2) * (1 - y))), 0.0)
#     return gamma_values

# def gamma(eps1, eps2, y, v1, phi1, v2, phi2):
#     y_squared = y*y
#     return np.where(eps1 + eps2 < (1 - y_squared) / (2 * (1 - y_squared / 2)),
#                     np.maximum(0, -0.5 * np.sin(2 * v1) * np.cos(phi1) * np.sqrt(eps1)
#                                - 0.5 * np.sin(2 * v2) * np.cos(phi2) * np.sqrt(eps2)
#                                + np.sqrt((1 - eps1 - eps2) * (1 - y))), 0)
#


def gamma(eps1, eps2, y, v1, phi1, v2, phi2):

    return np.maximum(0, -0.5 * np.sin(2 * v1) * np.cos(phi1) * np.sqrt(eps1)
                      - 0.5 * np.sin(2 * v2) * np.cos(phi2) * np.sqrt(eps2)
                      + np.sqrt((1 - eps1 - eps2) * (1 - y)))

# безразмерная

def g(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u):
    return (np.power((np.cos(v1) * np.cos(phi1) * np.cos(v2) * np.cos(phi2)
            * (gamma(eps1, eps2, y, v1, phi1, v2, phi2) * alpha * u * 0.5)), 2)
            * (1 / (omega1 * omega2)) * np.power((ksi / np.sinh(ksi)), 2))

# безразмерная


# def g(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u):
#     return np.where(ksi != 0, np.power((np.cos(v1) * np.cos(phi1) * np.cos(v2) * np.cos(phi2)
#                                          * (gamma(eps1, eps2, y, v1, phi1, v2, phi2) * alpha * u * 0.5)), 2)
#                     * (1 / (omega1 * omega2)) * np.power((ksi / np.sinh(ksi)), 2),
#                     (np.power((np.cos(v1) * np.cos(phi1) * np.cos(v2) * np.cos(phi2)
#                                * (gamma(eps1, eps2, y, v1, phi1, v2, phi2) * alpha * u * 0.5)), 2)
#                      * (1 / (omega1 * omega2)), 2))

def p_vv(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2):
    # s = np.absolute(i1 - f1)
    s = np.absolute(i2 - f2)

    m_r = (m1.mass*m2.mass)/(m1.mass+m2.mass)  # приведенная масса, kg
    el_lvl = 1 - 1  # electronic level

    e1_1 = m1.ev_i[el_lvl][i1]  # initial state, J
    e1_2 = m1.ev_i[el_lvl][f1]  # final state, J
    e2_1 = m2.ev_i[el_lvl][i2]  # initial state, J
    e2_2 = m2.ev_i[el_lvl][f2]  # final state, J

    omega1 = np.absolute(e1_1 - e1_2) / (s*h_red)       # J/(J*s)=1/s
    omega2 = np.absolute(e2_1 - e2_2) / (s*h_red)       # J/(J*s)=1/s

    u = np.sqrt(2*E/m_r)        # sqrt(J/kg)=sqrt(m^2/s^2)=m/s

    ksi = (pi*(omega1 - omega2)) / (alpha*u)  # (1/s) / (1/m) / (m/s) = 1

    ns1 = np.power((factorial(max(i1, f1)) / factorial(min(i1, f1))), (1 / s))
    ns2 = np.power((factorial(max(i2, f2)) / factorial(min(i2, f2))), (1 / s))

    return (np.power(ns1*ns2*g(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u), s)
            / (np.power(factorial(s), 2))
            * np.exp(-(2*ns1*g(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u) / (s+1))
                     - (np.power(ns1*g(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u)/(s+1), 2)
                        / (s+2))))

# print (g(0.2, 0.2, 0.2, pi/4, pi/4, pi/4, pi/4, -0.2053596293420575,
#          425390684400060.25 , 650224277299828.6, 1400))







