
import numpy as np

from scipy.integrate import quad

from p_vv_mm import p_vv_int
from constants import k, h, c


def k_vv_mm_n2_n2_for_v(v: int, T: float = 300.0) -> float:
    """
    N2–N2 V–V: (v1, i1)=(1)→(0), (v2)=(v−1)→(v). См³/с.
    Вынесено в модуль для ProcessPoolExecutor.map (pickle).
    """
    from particles_data import N2

    return k_vv_mm(N2, N2, i1=1, i2=v - 1, f1=0, f2=v, T=T)


def k_vv_mm_o2_o2_for_v(v: int, T: float = 300.0) -> float:
    """
    O2–O2 V–V: i1=v, f1=v−1, i2=0, f2=1. См³/с.
    Вынесено в модуль для ProcessPoolExecutor.map (pickle).
    """
    from particles_data import O2

    return k_vv_mm(O2, O2, i1=v, f1=v - 1, i2=0, f2=1, T=T)


def k_vv_mm(m1, m2, i1, f1, i2, f2, T):
    T_K = T
    T_inv_cm = T_K * k / (h * c * 100)

    m_red = (m1.mass * m2.mass) / (m1.mass + m2.mass)

    # ---------- ИСПРАВЛЕНИЕ: вычисление R0 ----------
    # Параметры потенциала (одинаковые для N2-N2, O2-O2)
    A_eV = 1730.0                     # эВ
    A_J = A_eV * 1.602176634e-19      # Дж
    alpha_m = 4e10                  # 1/м

    # Потенциал U(R) = 4 * A * exp(-alpha * R)
    U0 = 4.0 * A_J

    # R0 из условия U(R0) = kT
    R0 = (1.0 / alpha_m) * np.log(U0 / (k * T_K))
    # ------------------------------------------------

    # Средняя скорость (как раньше)
    mean_u = np.sqrt(8 * k * T_K / (np.pi * m_red))

    # Энергетический дефект ΔE (Дж → см⁻¹)
    el_lvl = 0
    e_i1 = m1.ev_i[el_lvl][i1]
    e_f1 = m1.ev_i[el_lvl][f1]
    e_i2 = m2.ev_i[el_lvl][i2]
    e_f2 = m2.ev_i[el_lvl][f2]
    delta_E_J = (e_i1 + e_i2) - (e_f1 + e_f2)
    delta_E_cm = delta_E_J / (h * c * 100)   # ΔE в см⁻¹

    def integrand(x, *args):
        # x = E_bar / (kT)
        m1_, m2_, i1_, f1_, i2_, f2_, T_inv_cm_, dE_cm = args

        # Симметризованная энергия (см⁻¹)
        E_bar_cm = x * T_inv_cm_

        # Вероятность считается от СИММЕТРИЗОВАННОЙ энергии,
        # как того требует формула (20)
        prob = p_vv_int(m1_, m2_, i1_, f1_, i2_, f2_, E_bar_cm, 'trapez')

        return prob * np.exp(-x) * x**3

    args = (m1, m2, i1, f1, i2, f2, T_inv_cm, delta_E_cm)
    result, _ = quad(integrand, 0, np.inf, args=args, epsabs=1e-15, epsrel=1e-15, limit=10000)

    # Константа скорости в см³/с
    k_vv = np.pi * R0**2 * mean_u * result * 1e6
    return k_vv