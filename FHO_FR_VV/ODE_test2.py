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
    x = T ** (-1 / 3)
    return np.exp(sum(c * x ** i for i, c in enumerate(coeffs)))


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

        if v + 1 <= v_max and v_ - 1 >= 0:
            k4 = k_lookup(coeffs_4d, v, v + 1, v_, v_ - 1, T)

        Nvp1 = N[v + 1] if v + 1 <= v_max else 0
        Nvm1 = N[v - 1] if v - 1 >= 0 else 0

        R += (k1 * Nvp1 + k2 * Nvm1 - (k3 + k4) * N[v]) * N[v_]

    return R * 8


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


def rhs_fast_with_inflow(t, N_flat, r, dr, v_max, m1, m2, T, D, N_total):
    N = N_flat.reshape(v_max + 1, -1)
    Nr = len(r)
    dNdt = np.zeros_like(N)

    # VV-обмен
    for i in range(Nr):
        N_local = N[:, i]
        for v in range(v_max + 1):
            vv = R_VV_fast(m1, m2, v, N_local, v_max, T)
            dNdt[v, i] = vv

    # Диффузия с притоком из резервуара
    for v in range(v_max + 1):
        Nv = N[v, :]
        for i in range(Nr):
            if i == 0:
                diff = D * 2 * (Nv[1] - Nv[0]) / dr ** 2
            elif i == Nr - 1:
                # На границе: восстанавливаем равновесное распределение
                # Характерное время диффузии на длине ячейки
                L = r[-1]  # размер ячейки
                tau_diff = L ** 2 / D  # ~ время перемешивания во всей ячейке

                if v == 0:
                    # v=0 возвращается к полной плотности
                    N_eq = N_total
                else:
                    # возбужденные уходят в ноль
                    N_eq = 0.0

                # Приток/отток через границу
                diff = (N_eq - Nv[-1]) / tau_diff
            else:
                diff = D * ((Nv[i + 1] - 2 * Nv[i] + Nv[i - 1]) / dr ** 2 +
                            (Nv[i + 1] - Nv[i - 1]) / (2 * r[i] * dr))
            dNdt[v, i] += diff

    return dNdt.flatten()


# Параметры
D = 0.2  # коэффициент диффузии, см²/с

t_span = (0, 12e-6)
t_eval = np.linspace(0, 12e-6, 600)

v_max = 5
m1 = O2
m2 = O2
T = 300

# Геометрия
a = 80e-4  # радиус пучка, см
# Rmax = 5 * a
# Nr = 50
# r = np.linspace(0, Rmax, Nr)
Rmax = 0.5  # 2 мм (вместо 0.03 см = 300 мкм)
Nr = 200
r = np.linspace(0, Rmax, Nr)
dr = r[1] - r[0]

# Начальные условия
N0_rel = 0.63
N1_rel = 0.37
N_total = 5.14e18  # полная плотность O2, см⁻³

# Гауссов профиль
profile = np.exp(-2 * r ** 2 / a ** 2)

N0 = np.zeros((v_max + 1, Nr))
for i in range(Nr):
    N0[0, i] = N0_rel * profile[i] * N_total
    N0[1, i] = N1_rel * profile[i] * N_total

N0_flat = N0.flatten()

# Решение
sol = solve_ivp(
    lambda t, N: rhs_fast(t, N, r, dr, v_max, m1, m2, T, D),
    t_span,
    N0_flat,
    t_eval=t_eval,
    method='Radau'
)


# ========== НОРМИРОВКА С УЧЕТОМ ДИФФУЗИИ ==========
# Используем probe-weighted population (как в эксперименте)
# Это моделирует то, что видит измерительный прибор

def probe_weighted_population(N_sol, r, sigma_probe):
    """
    Усреднение с весом probe beam (Гаусс)
    sigma_probe - радиус probe пучка (1/e^2)
    """
    weight = np.exp(-2.0 * (r / sigma_probe) ** 2)
    denom = np.trapz(weight * r, r)
    nv, nr, nt = N_sol.shape
    obs = np.zeros((nv, nt))
    for v in range(nv):
        for it in range(nt):
            num = np.trapz(N_sol[v, :, it] * weight * r, r)
            obs[v, it] = num / denom if denom > 0 else 0.0
    return obs


# Радиус probe пучка (обычно меньше pump, но для простоты используем тот же)
sigma_probe = a  # 80 мкм

# Восстанавливаем решение
Nr = len(r)
N_sol = sol.y.reshape(v_max + 1, Nr, -1)

# Наблюдаемые (взвешенные) концентрации
obs_N = probe_weighted_population(N_sol, r, sigma_probe)

# ========== КАЛИБРОВКА ПО НАЧАЛЬНОМУ СОСТОЯНИЮ ==========
# При t=0: f(0) = 0.63, f(1) = 0.37 (известно из эксперимента)
# Используем это для калибровки

f0_init_target = 0.63
f1_init_target = 0.37

# Калибровочный множитель (приводим obs_N к правильным начальным долям)
# При t=0: obs_N[0,0] * scale = f0_init_target
scale = f0_init_target / obs_N[0, 0]

# Калиброванные доли (учитывают диффузию)
fractional = obs_N * scale

# Диагностика
print("\n=== ДИАГНОСТИКА НОРМИРОВКИ (с учетом диффузии) ===")
print(f"Калибровочный множитель: {scale:.4f}")
print(f"f(0) при t=0: {fractional[0, 0]:.4f} (цель: {f0_init_target})")
print(f"f(1) при t=0: {fractional[1, 0]:.4f} (цель: {f1_init_target})")
print(f"Сумма f(v) при t=0: {np.sum(fractional[:, 0]):.4f}")
print(f"Сумма f(v) при t=10 мкс: {np.sum(fractional[:, -1]):.4f}")
print(f"Потеря сигнала за 10 мкс: {(1 - np.sum(fractional[:, -1]) / np.sum(fractional[:, 0])) * 100:.2f}%")


# Функция для чтения CSV с заменой запятых на точки
def read_csv_with_comma_fix(filename, sep=';'):
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    content_fixed = content.replace(',', '.')
    from io import StringIO
    return pd.read_csv(StringIO(content_fixed), sep=sep, header=None)


# Графики
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()

# Загрузка экспериментальных данных
data_11 = read_csv_with_comma_fix('v0solid.csv', sep=';')
data_12 = read_csv_with_comma_fix('v0dashed.csv', sep=';')
data_21 = read_csv_with_comma_fix('v1solid.csv', sep=';')
data_22 = read_csv_with_comma_fix('v1dashed.csv', sep=';')
data_31 = read_csv_with_comma_fix('v2solid.csv', sep=';')
data_32 = read_csv_with_comma_fix('v2dashed.csv', sep=';')
data_41 = read_csv_with_comma_fix('v3solid.csv', sep=';')
data_42 = read_csv_with_comma_fix('v3dashed.csv', sep=';')
data_51 = read_csv_with_comma_fix('v4solid.csv', sep=';')
data_52 = read_csv_with_comma_fix('v4dashed.csv', sep=';')
data_61 = read_csv_with_comma_fix('v5solid.csv', sep=';')
data_62 = read_csv_with_comma_fix('v5dashed.csv', sep=';')

exp_data_pairs = [
    [(data_11, 'solid'), (data_12, 'dashed')],
    [(data_21, 'solid'), (data_22, 'dashed')],
    [(data_31, 'solid'), (data_32, 'dashed')],
    [(data_41, 'solid'), (data_42, 'dashed')],
    [(data_51, 'solid'), (data_52, 'dashed')],
    [(data_61, 'solid'), (data_62, 'dashed')],
]

exp_colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']

for v in range(min(6, v_max + 1)):
    axes[v].plot(sol.t * 1e6, fractional[v, :], linewidth=2, color='black', label='Расчёт')

    for idx, (data, style) in enumerate(exp_data_pairs[v]):
        if data is not None and len(data) > 0:
            time_exp = data.iloc[:, 0].values.astype(float)
            pop_exp = data.iloc[:, 1].values.astype(float)
            axes[v].scatter(time_exp, pop_exp, marker='o', s=30,
                            color=exp_colors[idx], alpha=0.7,
                            label=f'Статья ({style})')
            axes[v].plot(time_exp, pop_exp, linestyle='--', linewidth=1,
                         color=exp_colors[idx], alpha=0.5)

    axes[v].set_xlabel('Время, мкс')
    axes[v].set_ylabel('Fractional population')
    axes[v].set_title(f'v={v}')
    axes[v].grid(True, alpha=0.3)
    axes[v].legend(fontsize=8)

for v in range(2, 6):
    axes[v].set_xlim(0, 6)

plt.tight_layout()
plt.show()

# Вывод численных значений
print("\n=== Fractional populations at key times ===")
time_indices = [0, len(sol.t) // 10, len(sol.t) // 5, len(sol.t) // 2, -1]
for ti in time_indices:
    print(f"\nt = {sol.t[ti] * 1e6:.2f} мкс:")
    for v in range(min(6, v_max + 1)):
        print(f"  v={v}: {fractional[v, ti]:.6f}")
    print(f"  сумма: {np.sum(fractional[:, ti]):.6f}")


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


print_my_table(COEFFS, v_max=5, T=300)