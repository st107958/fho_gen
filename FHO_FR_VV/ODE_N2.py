# from R_VV import R_VV
from FHO_FR_VV.particles_data import *
from k_vv_mm import k_vv_mm

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import json
import pandas as pd


with open('coeffs_sparse_N2.json', 'r') as f:
    COEFFS = json.load(f)

def k_lookup(coeffs_dict, v1, v2, v1p, v2p, T):

    key = f"{v1}_{v2}_{v1p}_{v2p}"
    coeffs = coeffs_dict.get(key)

    if coeffs is None:
        return 0.0

    x = T**(-1/3)
    return np.exp(sum(c * x**i for i, c in enumerate(coeffs)))


def R_VV_fast(m1, m2, v, N, v_max, T, coeffs_4d=COEFFS):
    R = 0

    for v_ in range(0, v_max + 1):

        k1 = k2 = k3 = k4 = 0

        if v + 1 <= v_max and v_ + 1 <= v_max:
            k1 = k_lookup(coeffs_4d, v + 1, v, v_, v_ + 1, T)

        if v - 1 >= 0 and v_ - 1 >= 0:
            k2 = k_lookup(coeffs_4d, v - 1, v, v_, v_ - 1, T)

        if v - 1 >= 0 and v_ + 1 <= v_max:
            k3 = k_lookup(coeffs_4d, v, v - 1, v_, v_ + 1, T)
            # print(k3)

        if v + 1 <= v_max and v_ - 1 >= 0:
            k4 = k_lookup(coeffs_4d, v, v + 1, v_, v_ - 1, T)

        Nvp1 = N[v + 1] if v + 1 <= v_max else 0
        Nvm1 = N[v - 1] if v - 1 >= 0 else 0


        R += (k1 * Nvp1 + k2 * Nvm1 - (k3 + k4) * N[v]) * N[v_]


        # print(k1, k2, k3, k4)

    return R # * 8 # * 1000

# def rhs_fast(t, N_flat, r, dr, v_max, m1, m2, T, D):
#     N = N_flat.reshape(v_max + 1, -1)  # [v, i]
#     Nr = len(r)
#     dNdt = np.zeros_like(N)
#
#     for i in range(Nr):  # по радиусу
#         N_local = N[:, i]  # все уровни в точке i
#         for v in range(v_max + 1):
#             vv = R_VV_fast(m1, m2, v, N_local, v_max, T)
#             dNdt[v, i] = vv
#
#     # Диффузия
#     for v in range(v_max + 1):
#         Nv = N[v, :]
#         for i in range(Nr):
#             if i == 0:
#                 diff = D * 2 * (Nv[1] - Nv[0]) / dr ** 2
#             elif i == Nr - 1:
#                 diff = 0
#             else:
#                 diff = D * ((Nv[i + 1] - 2 * Nv[i] + Nv[i - 1]) / dr ** 2 +
#                             (Nv[i + 1] - Nv[i - 1]) / (2 * r[i] * dr))
#             dNdt[v, i] += diff
#
#     return dNdt.flatten()


def rhs_fast(t, N_flat, r, dr, v_max, m1, m2, T, D):
    N = N_flat.reshape(v_max + 1, -1)
    Nr = len(r)
    dNdt = np.zeros_like(N)

    # VV-член
    for i in range(Nr):
        N_local = N[:, i]
        for v in range(v_max + 1):
            vv = R_VV_fast(m1, m2, v, N_local, v_max, T)
            dNdt[v, i] = vv

    # Диффузия
    if D > 0:
        for v in range(v_max + 1):
            Nv = N[v, :]
            for i in range(Nr):
                if i == 0:
                    diff = D * 2 * (Nv[1] - Nv[0]) / dr ** 2
                elif i == Nr - 1:
                    diff = 0
                else:
                    diff = D * ((Nv[i + 1] - 2 * Nv[i] + Nv[i - 1]) / dr ** 2 +
                                (Nv[i + 1] - Nv[i - 1]) / (2 * r[i] * dr))
                dNdt[v, i] += diff

        # Отладка: вывести максимальную диффузию для v=1
        # if int(t * 1e6) % 2 == 0:  # раз в 2 мкс
        #     max_diff = np.max(np.abs(dNdt[1, :]))
        #     print(f"t={t * 1e6:.1f} мкс, max(diff) = {max_diff:.2e}")

    return dNdt.flatten()

#D = 0
D = 0.4    # D ~0.2 см²/с


t_span = (0, 12e-6)  # 65 нс - 5 мкс (как в статье)
t_eval = np.linspace(0, 12e-6, 600)  # больше точек

# r = np.linspace(0, 0.04, 40)
# dr = r[1] - r[0]

v_max = 6

m1 = N2
m2 = N2

T = 300

sigma = 65e-4
Rmax = 5 * sigma
Nr = 50
r = np.linspace(0, Rmax, Nr)
dr = r[1] - r[0]

# Относительные заселённости
N0_rel = 0.63
N1_rel = 0.37

# N_total = 5.14 * 1e18
P_torr = 300  # или 300, 520
P_atm = P_torr / 760
N_total = 2.45e19 * P_atm

# Гауссов профиль
profile = np.exp(-2*r**2 / sigma**2)

# Массив [v, i]
N0 = np.zeros((v_max + 1, Nr))

for v in range(v_max + 1):
    for i in range(Nr):
        if v == 0:
            N0[v, i] = N0_rel * profile[i] * N_total
        elif v == 1:
            N0[v, i] = N1_rel * profile[i] * N_total
        else:
            N0[v, i] = 0.0

# Для solve_ivp
N0_flat = N0.flatten()

sol = solve_ivp(
    lambda t, N: rhs_fast(t, N, r, dr,  v_max, m1, m2, T, D),
    t_span,
    N0_flat,
    t_eval=t_eval,
    method='Radau'          #'BDF'
)

# Восстанавливаем решение и считаем fractional populations
N_sol = sol.y.reshape(v_max + 1, Nr, -1)
full_N = np.array([[2 * np.pi * np.sum(N_sol[v, :, t] * r * dr) for t in range(len(sol.t))] for v in range(v_max + 1)])
fractional = full_N / np.sum(full_N, axis=0)


# Функция для чтения CSV с заменой запятых на точки
def read_csv_with_comma_fix(filename, sep=';'):
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.read()
        content_fixed = content.replace(',', '.')
        from io import StringIO
        return pd.read_csv(StringIO(content_fixed), sep=sep, header=None)
    except FileNotFoundError:
        print(f"Файл {filename} не найден")
        return None


# Один рисунок с несколькими подграфиками
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

# Загрузка экспериментальных данных (только те, что есть)
exp_data = {
    0: read_csv_with_comma_fix('expv0.csv', sep=';'),
    1: read_csv_with_comma_fix('expv1.csv', sep=';'),
    2: read_csv_with_comma_fix('expv2.csv', sep=';'),
    3: read_csv_with_comma_fix('expv3.csv', sep=';'),
}

# Цвет для экспериментальных точек
exp_color = 'red'

for v in range(v_max + 1):
    # Расчётная кривая
    axes[v].plot(sol.t * 1e6, fractional[v, :], linewidth=2, color='black', label='Расчёт (FHO-FR)')

    # Экспериментальные данные (если есть)
    if v in exp_data and exp_data[v] is not None and len(exp_data[v]) > 0:
        time_exp = exp_data[v].iloc[:, 0].values.astype(float)
        pop_exp = exp_data[v].iloc[:, 1].values.astype(float)

        # Точки
        axes[v].scatter(time_exp, pop_exp, marker='o', s=40,
                        color=exp_color, alpha=0.8, zorder=5, label='Эксперимент (Ahn et al.)')


    axes[v].set_xlabel('Время, мкс', fontsize=10)
    axes[v].set_ylabel('Fractional population', fontsize=10)
    axes[v].set_title(f'v = {v}', fontsize=12)
    axes[v].grid(True, alpha=0.3)
    axes[v].legend(fontsize=8)

    # Настройка пределов
    if v >= 2:
        axes[v].set_xlim(0, 6)
    if v >= 3:
        axes[v].set_ylim(0, max(fractional[v, :]) * 1.2)

# Убираем пустые подграфики если v_max < 7
for v in range(v_max + 1, len(axes)):
    axes[v].set_visible(False)

plt.tight_layout()
plt.savefig('N2_VV_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

# Вывод численных значений в консоль
print("\n" + "=" * 70)
print("Fractional populations for N₂ at 300K, 300 Torr")
print("=" * 70)
print(f"{'Time (μs)':^12}", end="")
for v in range(v_max + 1):
    print(f"v={v}:^10", end="")
print()
print("-" * 70)

# Выбираем несколько моментов времени
time_indices = [0, len(sol.t) // 10, len(sol.t) // 5, len(sol.t) // 3,
                len(sol.t) // 2, 2 * len(sol.t) // 3, -1]
for ti in time_indices:
    t_us = sol.t[ti] * 1e6
    if t_us > 12:
        continue
    print(f"{t_us:^12.2f}", end="")
    for v in range(v_max + 1):
        print(f"{fractional[v, ti]:^10.5f}", end="")
    print()

print("=" * 70)