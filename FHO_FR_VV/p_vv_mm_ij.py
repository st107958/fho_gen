from constants import *

import numpy as np
from scipy.special import factorial


def p_vv(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2):
    s = np.absolute(i1 - f1)
    m_r = m1.mass*m2.mass/(m1.mass+m2.mass)
    el_lvl = 1 - 1  # electronic level

    e1_1 = m1.ev_i[el_lvl][i1+1]  # initial state, J
    e1_2 = m1.ev_i[el_lvl][f1+1]  # final state, J
    e2_1 = m2.ev_i[el_lvl][i2+1]  # initial state, J
    e2_2 = m2.ev_i[el_lvl][f2+1]  # final state, J

    omega1 = np.absolute(e1_1 - e1_2) / (s*h_red)
    omega2 = np.absolute(e2_1 - e2_2) / (s*h_red)

    def u(E, m):
        return np.sqrt(2*E/m)

    ksi = pi*(omega1 - omega2) / (alpha*u(E, m_r))

    def gamma(eps1, eps2, y):
        return np.where(eps1 + eps2 < 1, np.maximum(0, -0.5*np.sin(2*v1)*np.cos(phi1)*np.sqrt(eps1)
                                                    - 0.5*np.sin(2*v1)*np.cos(phi1)*np.sqrt(eps1)
                                                    + np.sqrt((1-eps1-eps2)*(1-y))), 0)

    def g(E, e1, e2, m, v1, phi1, v2, phi2):
        return (((np.cos(v1)*np.cos(phi1)*np.cos(v2)*np.cos(phi2)*(gamma(e1, e2, y)*alpha*u(E, m)*0.5))**2)
                / (omega1*omega2)) / ((ksi*np.sinh(ksi))**2)

    ns1 = np.power(factorial(max(i1, f1)) / factorial(min(i1, f1)), 1 / s)
    ns2 = np.power(factorial(max(i2, f2)) / factorial(min(i2, f2)), 1 / s)

    return (np.power(ns1*ns2*g(E, eps1, eps2, m_r, v1, phi1, v2, phi2), s) / (np.power(factorial(s), 2))
            * np.exp(-(2*ns1*g(E, eps1, eps2, m_r, v1, phi1, v2, phi2) / (s+1))) -
            (np.power(ns1*g(E, eps1, eps2, m_r, v1, phi1, v2, phi2), 2) / (((s+1)**2) * (s+2))))


# print(p_vv(N2, N2, ))