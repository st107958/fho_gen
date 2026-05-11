from constants import alpha, h_red, pi

import numpy as np
from scipy.special import factorial


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

    # --- вспомогательная функция для безопасного (ksi/sinh(ksi))^2 ---
    def safe_factor(x):
        ax = np.abs(x)
        # При малых x -> 1
        if ax < 1e-4:
            return 1.0
        # При больших x используем асимптотику: (2x e^{-x})^2
        # if ax > 10:
        #     return (2 * ax * np.exp(-ax)) ** 2
        # Основной диапазон
        return (x / np.sinh(x)) ** 2

    # векторизуем, если ksi – массив
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
    """
    Вероятность V‑V перехода (7) для фиксированных переменных.
    E – энергия в Джоулях.
    Возвращает число или массив.
    """
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

    omega1 = abs(e1_1 - e1_2) / (s * h_red)   # 1/с
    omega2 = abs(e2_1 - e2_2) / (s * h_red)

    u = np.sqrt(2 * E / m_r)
    ksi = (pi * (omega1 - omega2)) / (alpha * u)

    ns1 = (factorial(max(i1, f1)) / factorial(min(i1, f1))) ** (1.0 / s)
    ns2 = (factorial(max(i2, f2)) / factorial(min(i2, f2))) ** (1.0 / s)
    ns  = np.sqrt(ns1 * ns2)                     # <-- ИСПРАВЛЕНО
    # ns  = ns1 * ns2

    G = G_func(y, eps1, eps2, v1, phi1, v2, phi2, ksi, omega1, omega2, u)

    term1 = (ns1 * ns2 * G) ** s / (factorial(s) ** 2)
    exponent = -(2 * ns * G) / (s + 1) - (ns * G / (s + 1)) ** 2 / (s + 2)   # <-- ИСПРАВЛЕНО
    # exponent = -(2 * ns * G) / (s + 1) - (ns * G / (s + 1)) ** 2 / (s + 2)
    return term1 * np.exp(exponent)





