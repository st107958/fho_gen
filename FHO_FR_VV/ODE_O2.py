# from R_VV import R_VV
from FHO_FR_VV.particles_data import *
from k_vv_mm import k_vv_mm

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import json
import pandas as pd


with open('coeffs_sparse2.json', 'r') as f:
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

    return R  * 8 # * 1000


# def rhs(t, N_flat, r, dr, v_max, m1, m2, T, D):
#     N = N_flat.reshape(v_max + 1, -1)  # [v, i]
#     Nr = len(r)
#     dNdt = np.zeros_like(N)
#
#     for i in range(Nr):  # по радиусу
#         N_local = N[:, i]  # все уровни в точке i
#         for v in range(v_max + 1):
#             vv = R_VV(m1, m2, v, N_local, v_max, T)
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
    N = N_flat.reshape(v_max + 1, -1)  # [v, i]
    Nr = len(r)
    dNdt = np.zeros_like(N)

    for i in range(Nr):  # по радиусу
        N_local = N[:, i]  # все уровни в точке i
        for v in range(v_max + 1):
            vv = R_VV_fast(m1, m2, v, N_local, v_max, T)
            dNdt[v, i] = vv

    # Диффузия
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

    return dNdt.flatten()



#D = 0
D = 0.2    # D ~0.2 см²/с

# t_span = (65e-9, 0.7e-6)  # 65 нс - 2 мкс
# t_eval = np.linspace(65e-9, 0.7e-6, 200)

# t_span = (65e-9, 12e-6)  # 65 нс - 5 мкс (как в статье)
# t_eval = np.linspace(65e-9, 12e-6, 600)  # больше точек

t_span = (0, 12e-6)  # 65 нс - 5 мкс (как в статье)
t_eval = np.linspace(0, 12e-6, 600)  # больше точек

# r = np.linspace(0, 0.04, 40)
# dr = r[1] - r[0]

v_max = 5

m1 = O2
m2 = O2

T = 300

a = 80e-4
Rmax = 5 * a
Nr = 50
r = np.linspace(0, Rmax, Nr)
dr = r[1] - r[0]

# Относительные заселённости
N0_rel = 0.63
N1_rel = 0.37

N_total = 5.14 * 1e18

# Гауссов профиль
profile = np.exp(-2*r**2 / a**2)

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
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    # Заменяем запятые на точки (но не разделители)
    content_fixed = content.replace(',', '.')
    from io import StringIO
    return pd.read_csv(StringIO(content_fixed), sep=sep, header=None)


# Один рисунок с несколькими подграфиками
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()

# Загрузка экспериментальных данных с заменой запятых на точки
# Для v=0
data_11 = read_csv_with_comma_fix('v0solid.csv', sep=';')
data_12 = read_csv_with_comma_fix('v0dashed.csv', sep=';')
# Для v=1
data_21 = read_csv_with_comma_fix('v1solid.csv', sep=';')
data_22 = read_csv_with_comma_fix('v1dashed.csv', sep=';')
# Для v=2
data_31 = read_csv_with_comma_fix('v2solid.csv', sep=';')
data_32 = read_csv_with_comma_fix('v2dashed.csv', sep=';')
# Для v=3
data_41 = read_csv_with_comma_fix('v3solid.csv', sep=';')
data_42 = read_csv_with_comma_fix('v3dashed.csv', sep=';')
# Для v=4
data_51 = read_csv_with_comma_fix('v4solid.csv', sep=';')
data_52 = read_csv_with_comma_fix('v4dashed.csv', sep=';')
# Для v=5
data_61 = read_csv_with_comma_fix('v5solid.csv', sep=';')
data_62 = read_csv_with_comma_fix('v5dashed.csv', sep=';')

# Список пар экспериментальных данных для всех v=0..5
exp_data_pairs = [
    [(data_11, 'solid'), (data_12, 'dashed')],  # v=0
    [(data_21, 'solid'), (data_22, 'dashed')],  # v=1
    [(data_31, 'solid'), (data_32, 'dashed')],  # v=2
    [(data_41, 'solid'), (data_42, 'dashed')],  # v=3
    [(data_51, 'solid'), (data_52, 'dashed')],  # v=4
    [(data_61, 'solid'), (data_62, 'dashed')],  # v=5
]

# Цвета для разных экспериментов на одном графике
exp_colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']

for v in range(min(6, v_max + 1)):
    # Расчётная кривая
    axes[v].plot(sol.t * 1e6, fractional[v, :], linewidth=2, color='black', label='Расчёт')

    # Добавляем экспериментальные данные для ВСЕХ v (0..5)
    for idx, (data, style) in enumerate(exp_data_pairs[v]):
        if data is not None and len(data) > 0:
            # Предполагаем, что в CSV: первый столбец - время (мкс), второй - population
            time_exp = data.iloc[:, 0].values.astype(float)
            pop_exp = data.iloc[:, 1].values.astype(float)
            axes[v].scatter(time_exp, pop_exp, marker='o', s=30,
                            color=exp_colors[idx], alpha=0.7,
                            label=f'Статья ({style})')
            # Соединяем точки линией
            axes[v].plot(time_exp, pop_exp, linestyle='--', linewidth=1,
                         color=exp_colors[idx], alpha=0.5)

    axes[v].set_xlabel('Время, мкс')
    axes[v].set_ylabel('Fractional population')
    axes[v].set_title(f'v={v}')
    axes[v].grid(True, alpha=0.3)
    axes[v].legend(fontsize=8)

for v in range(2, 6):
    axes[v].set_xlim(0, 6)      # показывать только до 6 мкс


plt.tight_layout()
plt.show()

# Вывод численных значений
print("\n=== Fractional populations at key times ===")
time_indices = [0, len(sol.t) // 10, len(sol.t) // 5, len(sol.t) // 2, -1]
for ti in time_indices:
    print(f"\nt = {sol.t[ti] * 1e6:.2f} мкс:")
    for v in range(min(6, v_max + 1)):
        print(f"  v={v}: {fractional[v, ti]:.6f}")

#таблицпа
def print_my_table(coeffs_dict, v_max=5, T=300):
    print("\n" + "=" * 60)
    print("O₂ rate coefficients k_{0,1}^{v,v-1} × 10¹³ cm³ s⁻¹")
    print("=" * 60)
    print(f"{'v':^10}", end="")
    for v in range(2, v_max + 1):
        print(f"{v:^12}", end="")
    print()
    print("-" * 60)

    print(f"{'This work':^10}", end="")

    for v in range(2, v_max + 1):
        # Процесс: O₂(v) + O₂(0) → O₂(v-1) + O₂(1)
        # Ключ: v, 0, v-1, 1
        key = f"{0}_{1}_{v}_{v - 1}"

        coeffs = coeffs_dict.get(key)

        if coeffs is not None and coeffs != 0:
            if isinstance(coeffs, list):
                x = T ** (-1 / 3)
                k = np.exp(sum(c * x ** i for i, c in enumerate(coeffs)))
            else:
                k = coeffs
            k_13 = k / 1e-13
            print(f"{k_13:^12.3f}", end="")
        else:
            print(f"{'N/A':^12}", end="")

    print()
    print("=" * 60)


# Запуск
print_my_table(COEFFS, v_max=5, T=300)



