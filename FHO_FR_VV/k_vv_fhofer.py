"""
Аналитическая скорость V–V по FHO-FR (порт fhovv из Fortran).

Оптимизации без изменения формул: таблица факториалов, кэш уровней в см⁻¹,
кэш stfac(s), скалярный math вместо scipy/numpy в горячем пути.
"""
from __future__ import annotations

import math
from functools import lru_cache

# Физические константы (как в ноутбуке / fhovv)
K_B = 1.380649e-23
H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
HC_K = 1.4387769
J_TO_CM = 1.0 / (H_PLANCK * C_LIGHT * 100.0)

ALPHA_SI = 4.0 * 1e10
R0_FORT = 2.52e-10
C_CM_PER_S = C_LIGHT * 100.0
_TWO_PI_C_CM = 2.0 * math.pi * C_CM_PER_S

_MAX_FACTORIAL_N = 80
_FACTORIAL = tuple(float(math.factorial(n)) for n in range(_MAX_FACTORIAL_N))

_E_CM_CACHE: dict[int, tuple[float, ...]] = {}


def _factorial(n: int) -> float:
    if n < _MAX_FACTORIAL_N:
        return _FACTORIAL[n]
    return float(math.factorial(n))


def _ns_ratio(v1: int, v2: int, s: int) -> float:
    hi, lo = (v1, v2) if v1 >= v2 else (v2, v1)
    return (_factorial(hi) / _factorial(lo)) ** (1.0 / s)


@lru_cache(maxsize=16)
def _stfac(s: int) -> float:
    stfac = (0.25 * (1.0 + 1.0 / 2.0 ** (s - 1))) ** 4
    stfac /= (s + 1) ** 2 * (2.0**s)
    stfac *= _factorial(s + 3) / 6.0
    return stfac


def _E_cm_levels(mol) -> tuple[float, ...]:
    key = id(mol)
    cached = _E_CM_CACHE.get(key)
    if cached is None:
        cached = tuple(e * J_TO_CM for e in mol.ev_i[0])
        _E_CM_CACHE[key] = cached
    return cached


def _k_vv_fhofer_exo(m1, m2, i1, f1, i2, f2, T, s, mu):
    E1 = _E_cm_levels(m1)
    E2 = _E_cm_levels(m2)

    E_i1, E_f1 = E1[i1], E1[f1]
    E_i2, E_f2 = E2[i2], E2[f2]

    delta_E_cm = (E_i1 + E_i2) - (E_f1 + E_f2)
    delta_E_K = delta_E_cm * HC_K

    if i1 != f1:
        e1 = abs(E_i1 - E_f1) / s
    else:
        e1 = m1.we[0] - 2.0 * m1.wexe[0] * i1
    if i2 != f2:
        e2 = abs(E_i2 - E_f2) / s
    else:
        e2 = m2.we[0] - 2.0 * m2.wexe[0] * i2
    Evib_cm = 0.5 * (e1 + e2)

    omega_rad = _TWO_PI_C_CM * Evib_cm
    ro2 = ALPHA_SI**2 * K_B * T / (2.0 * mu * omega_rad**2)

    ns1 = _ns_ratio(i1, f1, s)
    ns2 = _ns_ratio(i2, f2, s)
    stfac = _stfac(s)

    theta_prime = 4.0 * math.pi**2 * mu * omega_rad**2 / (ALPHA_SI**2 * K_B)
    arg = math.sqrt(theta_prime / T) * abs(delta_E_K / s) / (Evib_cm * HC_K) / math.sqrt(8.0)
    popc = math.exp(-4.0 / 9.0 * arg)
    popc = 0.5 * (3.0 - popc) * popc

    fact_s = _factorial(s)
    zompl0 = ns1 * ns2
    zompl1 = zompl0 * ro2
    zompl2 = stfac * zompl1**s / (fact_s * fact_s)
    zompl3 = (1.0 + 2.0 * zompl1 / (s + 1.0) * stfac) ** (s + 4.0)
    rate_dimless = zompl2 / zompl3 * popc

    mean_v = math.sqrt(8.0 * K_B * T / (math.pi * mu))
    Z = 3.0 * math.pi * R0_FORT**2 * mean_v

    return Z * 1e6 * rate_dimless * math.exp(0.5 * delta_E_K / T)


def k_vv_fhofer(m1, m2, i1, f1, i2, f2, T):
    """
    Аналитическая скорость V-V перехода по FHO-FR (реализация fhovv).
    Возвращает коэффициент скорости в см³/с для экзотермического направления.
    """
    s = abs(i1 - f1)
    if s != abs(i2 - f2) or s == 0:
        raise ValueError("Переход не является V-V обменом одного порядка")

    mu = (m1.mass * m2.mass) / (m1.mass + m2.mass)

    E1 = _E_cm_levels(m1)
    E2 = _E_cm_levels(m2)
    delta_E_cm = (E1[i1] + E2[i2]) - (E1[f1] + E2[f2])

    if delta_E_cm < -1e-5:
        delta_E_K = delta_E_cm * HC_K
        rate_exo = _k_vv_fhofer_exo(m1, m2, f1, i1, f2, i2, T, s, mu)
        return rate_exo * math.exp(-delta_E_K / T)

    return _k_vv_fhofer_exo(m1, m2, i1, f1, i2, f2, T, s, mu)
