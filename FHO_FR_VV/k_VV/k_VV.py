"""
Численное вычисление k_VV методом матричных элементов (matrix method).

Содержит полную цепочку: p_vv → p_vv_int → k_vv_mm.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad, trapezoid
from scipy.special import factorial

from constants import alpha, h_red, pi, k, h, c


# ─── p_vv_mm_ij.py ────────────────────────────────────────

def gamma(eps1, eps2, y, v1, phi1, v2, phi2):
    return np.maximum(0, -0.5 * np.sin(2 * v1) * np.cos(phi1) * np.sqrt(eps1)
                      - 0.5 * np.sin(2 * v2) * np.cos(phi2) * np.sqrt(eps2)
                      + np.sqrt((1 - eps1 - eps2) * (1 - y)))


def G_func(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u):
    gamma_val = gamma(eps1, eps2, y, v1, phi1, v2, phi2)

    base = (
        np.cos(v1) * np.cos(phi1) *
        np.cos(v2) * np.cos(phi2) *
        (gamma_val * alpha * u * 0.5)
    )**2 / (omega1 * omega2)

    def safe_factor(x):
        ax = np.abs(x)
        if ax < 1e-4:
            return 1.0
        return (x / np.sinh(x)) ** 2

    if np.isscalar(ksi):
        factor = safe_factor(ksi)
        return base * factor
    else:
        ksi = np.asarray(ksi)
        result = np.zeros_like(ksi, dtype=float)
        for idx in np.ndindex(ksi.shape):
            result[idx] = base[idx] * safe_factor(ksi[idx])
        return result


def p_vv(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2):
    """Вероятность V-V перехода (7) для фиксированных переменных. E — Дж."""
    s_ = abs(i1 - f1)
    s  = abs(i2 - f2)
    if s_ != s:
        raise ValueError("Число переданных квантов должно совпадать: |i1-f1| == |i2-f2|")

    m_r = (m1.mass * m2.mass) / (m1.mass + m2.mass)
    el_lvl = 0

    e1_1 = m1.ev_i[el_lvl][i1]
    e1_2 = m1.ev_i[el_lvl][f1]
    e2_1 = m2.ev_i[el_lvl][i2]
    e2_2 = m2.ev_i[el_lvl][f2]

    omega1 = abs(e1_1 - e1_2) / (s * h_red)
    omega2 = abs(e2_1 - e2_2) / (s * h_red)

    u = np.sqrt(2 * E / m_r)
    ksi = (pi * (omega1 - omega2)) / (alpha * u)

    ns1 = (factorial(max(i1, f1)) / factorial(min(i1, f1))) ** (1.0 / s)
    ns2 = (factorial(max(i2, f2)) / factorial(min(i2, f2))) ** (1.0 / s)
    ns  = np.sqrt(ns1 * ns2)

    G = G_func(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u)

    term1 = (ns1 * ns2 * G) ** s / (factorial(s) ** 2)
    exponent = -(2 * ns * G) / (s + 1) - (ns * G / (s + 1)) ** 2 / (s + 2)
    return term1 * np.exp(exponent)


# ─── p_vv_mm.py ───────────────────────────────────────────

def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapez'):
    """Интегрирование p_vv по 7 угловым/энергетическим переменным."""
    e_in_J = h * c * 100
    E = E * e_in_J   # 1/cm --> J

    if method == 'trapez':
        maxdiv = 9

        eps1 = np.linspace(0, 1, maxdiv)
        eps2 = np.linspace(0, 1, maxdiv)
        y = np.linspace(0, 1, maxdiv)
        v1 = np.linspace(0, np.pi, maxdiv)
        phi1 = np.linspace(0, np.pi, maxdiv)
        v2 = np.linspace(0, np.pi, maxdiv)
        phi2 = np.linspace(0, np.pi, maxdiv)

        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(
            eps1, eps2, y, v1, phi1, v2, phi2, indexing='ij'
        )

        mask = (EPS1 + EPS2) <= 1/2

        EPS1_filtered = EPS1[mask]
        EPS2_filtered = EPS2[mask]
        Y_filtered = Y[mask]
        V1_filtered = V1[mask]
        PHI1_filtered = PHI1[mask]
        V2_filtered = V2[mask]
        PHI2_filtered = PHI2[mask]

        F = p_vv(m1, m2, i1, f1, i2, f2, E,
                 EPS1_filtered, EPS2_filtered, Y_filtered,
                 V1_filtered, PHI1_filtered, V2_filtered, PHI2_filtered)

        F_ = np.zeros_like(EPS1)
        F_[mask] = F
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


# ─── vv_rate.py ───────────────────────────────────────────

def k_vv_mm(m1, m2, i1, f1, i2, f2, T):
    """Численное k_VV методом матричных элементов (см³/с)."""
    T_K = T
    T_inv_cm = T_K * k / (h * c * 100)

    m_red = (m1.mass * m2.mass) / (m1.mass + m2.mass)

    A_eV = 1730.0
    A_J = A_eV * 1.602176634e-19
    alpha_m = 4e10

    U0 = 4.0 * A_J
    R0 = (1.0 / alpha_m) * np.log(U0 / (k * T_K))

    mean_u = np.sqrt(8 * k * T_K / (np.pi * m_red))

    el_lvl = 0
    e_i1 = m1.ev_i[el_lvl][i1]
    e_f1 = m1.ev_i[el_lvl][f1]
    e_i2 = m2.ev_i[el_lvl][i2]
    e_f2 = m2.ev_i[el_lvl][f2]
    delta_E_J = (e_i1 + e_i2) - (e_f1 + e_f2)
    delta_E_cm = delta_E_J / (h * c * 100)

    def integrand(x, *args):
        m1_, m2_, i1_, f1_, i2_, f2_, T_inv_cm_, dE_cm = args
        E_bar_cm = x * T_inv_cm_
        prob = p_vv_int(m1_, m2_, i1_, f1_, i2_, f2_, E_bar_cm, 'trapez')
        return prob * np.exp(-x) * x**3

    args = (m1, m2, i1, f1, i2, f2, T_inv_cm, delta_E_cm)
    result, _ = quad(integrand, 0, np.inf, args=args,
                     epsabs=1e-15, epsrel=1e-15, limit=10000)

    k_vv = np.pi * R0**2 * mean_u * result * 1e6
    return k_vv


def k_vv_mm_for_v(m1_name: str, m2_name: str, v: int, T: float,
                   *, i1: int = 1, f1: int = 0) -> float:
    """Picklable wrapper для k_vv_mm. Принимает имена частиц (строки)."""
    from particles import N2, O2, CO
    _REGISTRY = {"N2": N2, "O2": O2, "CO": CO}
    m1 = _REGISTRY[m1_name]
    m2 = _REGISTRY[m2_name]
    return k_vv_mm(m1, m2, i1=i1, i2=v - 1, f1=f1, f2=v, T=T)
