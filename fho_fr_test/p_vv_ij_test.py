import cupy as cp
from scipy.special import factorial  # можно заменить на cp.gamma(x+1) если нужно на GPU

from fho_fr_test.constants import *
from fho_fr_test.particles_data import *


def gamma(eps1, eps2, y, v1, phi1, v2, phi2):
    return cp.maximum(
        0,
        -0.5 * cp.sin(2 * v1) * cp.cos(phi1) * cp.sqrt(eps1)
        - 0.5 * cp.sin(2 * v2) * cp.cos(phi2) * cp.sqrt(eps2)
        + cp.sqrt((1 - eps1 - eps2) * (1 - y))
    )


def g(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u):
    cos_product = (
        cp.cos(v1) * cp.cos(phi1) * cp.cos(v2) * cp.cos(phi2)
    )
    gamma_val = gamma(eps1, eps2, y, v1, phi1, v2, phi2)
    return (
        (cos_product * gamma_val * alpha * u * 0.5) ** 2
        * (1 / (omega1 * omega2))
        * (ksi / cp.sinh(ksi)) ** 2
    )


def p_vv(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2):
    s_ = abs(i1 - f1)
    s = abs(i2 - f2)

    m_r = (m1.mass * m2.mass) / (m1.mass + m2.mass)

    el_lvl = 0
    e1_1 = m1.ev_i[el_lvl][i1]
    e1_2 = m1.ev_i[el_lvl][f1]
    e2_1 = m2.ev_i[el_lvl][i2]
    e2_2 = m2.ev_i[el_lvl][f2]

    omega1 = abs(e1_1 - e1_2) / (s * h_red)
    omega2 = abs(e2_1 - e2_2) / (s * h_red)
    # print(omega1, omega2)

    u = cp.sqrt(2 * E / m_r)
    ksi = (pi * (omega1 - omega2)) / (alpha * u)

    # print(f"i1: {i1}, f1: {f1}, s: {s}")
    # print(f"e1_1: {e1_1}, e1_2: {e1_2}, h_red: {h_red}")

    # Можно заменить на cp.power(cp.exp(cp.loggamma(...)), ...) для GPU, но scipy.factorial точнее
    ns1 = (factorial(max(i1, f1)) / factorial(min(i1, f1))) ** (1 / s)
    ns2 = (factorial(max(i2, f2)) / factorial(min(i2, f2))) ** (1 / s)

    g_val = g(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u)

    numerator = (ns1 * ns2 * g_val) ** s
    denominator = factorial(s) ** 2
    exp_part = cp.exp(
        -2 * ns1 * g_val / (s + 1)
        - (ns1 * g_val / (s + 1)) ** 2 / (s + 2)
    )

    return numerator / denominator * exp_part
