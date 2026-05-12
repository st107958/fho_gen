import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp

# -------------------------------------------------------------
# 1. Коэффициенты из regression_coefficients/*.csv (как в regression_new)
# -------------------------------------------------------------
_REG_DIR = Path(__file__).resolve().parent / "regression_coefficients"
PARTICLE1_NAME = "O2"
PARTICLE2_NAME = "O2"


def load_coeffs_from_regression_csv(csv_path: Path) -> dict[str, list[float]]:
    """
    Читает `{p1}_{p2}_regression_coefs.csv` со столбцами i1, f1, i2, f2, a_0, a_1, ...
    и возвращает словарь с ключами `i1_f1_i2_f2` — тот же формат, что у json sparse.
    """
    df = pd.read_csv(csv_path)
    a_cols = sorted(
        (c for c in df.columns if isinstance(c, str) and c.startswith("a_")),
        key=lambda c: int(c.split("_", 1)[1]),
    )
    out: dict[str, list[float]] = {}
    for _, row in df.iterrows():
        i1, f1, i2, f2 = (int(row["i1"]), int(row["f1"]), int(row["i2"]), int(row["f2"]))
        key = f"{i1}_{f1}_{i2}_{f2}"
        coeffs = [float(row[c]) for c in a_cols if pd.notna(row[c])]
        if coeffs:
            out[key] = coeffs
    return out


_csv_reg = _REG_DIR / f"{PARTICLE1_NAME}_{PARTICLE2_NAME}_regression_coefs.csv"
if not _csv_reg.is_file():
    raise FileNotFoundError(
        f"Ожидался файл {_csv_reg} (создаётся в regression_new в папке regression_coefficients)."
    )
COEFFS = load_coeffs_from_regression_csv(_csv_reg)


def k_lookup(coeffs_dict, v1, v2, v1p, v2p, T):
    """Как раньше: ключ v1_v2_v1p_v2p, ln k = sum_i a_i T^(-i/3)."""
    key = f"{v1}_{v2}_{v1p}_{v2p}"
    coeffs = coeffs_dict.get(key)
    if coeffs is None or len(coeffs) == 0:
        return 0.0
    x = T ** (-1.0 / 3.0)
    return float(np.exp(sum(c * x**i for i, c in enumerate(coeffs))))


# -------------------------------------------------------------
# 2. Диффузионный оператор (цилиндрическая симметрия)
# -------------------------------------------------------------
def diffusion_term(n, r, D, bc_outer="dirichlet"):
    """
    dn/dt = (1/r) d/dr (r D dn/dr) в консервативной форме:
    dn/dt = -1/r * d(rJ)/dr, J = -D dn/dr.
    """
    dr = r[1] - r[0]
    nr = len(r)
    res = np.zeros_like(n)
    r_face = 0.5 * (r[:-1] + r[1:])  # r_{i+1/2}
    j_face = -D * (n[1:] - n[:-1]) / dr

    # ось симметрии: dn/dr|0 = 0
    res[0] = 2.0 * D * (n[1] - n[0]) / dr**2

    # внутренние узлы
    for i in range(1, nr - 1):
        div_j = (r_face[i] * j_face[i] - r_face[i - 1] * j_face[i - 1]) / (r[i] * dr)
        res[i] = -div_j

    # внешняя граница
    if bc_outer == "neumann":
        j_out = 0.0
        div_j = (r[-1] * j_out - r_face[-1] * j_face[-1]) / (r[-1] * dr)
        res[-1] = -div_j
    elif bc_outer == "dirichlet":
        res[-1] = 0.0
    else:
        raise ValueError(f"Unknown bc_outer={bc_outer}")

    return res


# -------------------------------------------------------------
# 3. VV-член (локальный)
# -------------------------------------------------------------
def vv_term(n, v_max, T, coeffs_dict):
    """
    n[v] – населённости одного радиального узла (0..v_max)
    возвращает dn_v/dt от VV-обмена
    """
    n_v = len(n)
    vv = np.zeros(n_v)

    for v in range(n_v):
        R = 0.0
        for vp in range(n_v):
            k1 = k2 = k3 = k4 = 0.0

            if v + 1 <= v_max and vp + 1 <= v_max:
                k1 = k_lookup(coeffs_dict, v + 1, v, vp, vp + 1, T)
            if v - 1 >= 0 and vp - 1 >= 0:
                k2 = k_lookup(coeffs_dict, v - 1, v, vp, vp - 1, T)
            if v - 1 >= 0 and vp + 1 <= v_max:
                k3 = k_lookup(coeffs_dict, v, v - 1, vp, vp + 1, T)
            if v + 1 <= v_max and vp - 1 >= 0:
                k4 = k_lookup(coeffs_dict, v, v + 1, vp, vp - 1, T)

            n_vp1 = n[v + 1] if v + 1 <= v_max else 0.0
            n_vm1 = n[v - 1] if v - 1 >= 0 else 0.0

            R += (k1 * n_vp1 + k2 * n_vm1 - (k3 + k4) * n[v]) * n[vp]

        vv[v] =  R * 8.0 #* 2.7

    return vv

# def vv_term(n, v_max, T, coeffs_dict):
#     """
#     n[v] – населённости одного радиального узла (0..v_max)
#     возвращает dn_v/dt от VV-обмена
#     """
#     n_v = len(n)
#     vv = np.zeros(n_v)
#
#     # Множитель 1.0 (убрать 8.0*2.7 после отладки)
#     scale = 8.0  # Временно, для проверки знаков
#
#     for v in range(n_v):
#         R = 0.0
#         for vp in range(n_v):
#             # Приход в v из v-1: (v-1, vp) → (v, vp-1)
#             k_in_from_vm1 = 0.0
#             if v - 1 >= 0 and vp - 1 >= 0:
#                 k_in_from_vm1 = k_lookup(coeffs_dict, v - 1, v, vp, vp - 1, T)
#
#             # Приход в v из v+1: (v+1, vp) → (v, vp+1)
#             k_in_from_vp1 = 0.0
#             if v + 1 <= v_max and vp + 1 <= v_max:
#                 k_in_from_vp1 = k_lookup(coeffs_dict, v + 1, v, vp, vp + 1, T)
#
#             # Уход из v в v+1: (v, vp) → (v+1, vp-1)
#             k_out_to_vp1 = 0.0
#             if v + 1 <= v_max and vp - 1 >= 0:
#                 k_out_to_vp1 = k_lookup(coeffs_dict, v, v + 1, vp, vp - 1, T)
#
#             # Уход из v в v-1: (v, vp) → (v-1, vp+1)
#             k_out_to_vm1 = 0.0
#             if v - 1 >= 0 and vp + 1 <= v_max:
#                 k_out_to_vm1 = k_lookup(coeffs_dict, v, v - 1, vp, vp + 1, T)
#
#             # Населенности для приходных членов
#             n_vm1 = n[v - 1] if v - 1 >= 0 else 0.0
#             n_vp1 = n[v + 1] if v + 1 <= v_max else 0.0
#
#             # Ток: приход - уход
#             # Член с k_in_from_vm1 требует n[v-1] и n[vp]
#             # Член с k_in_from_vp1 требует n[v+1] и n[vp]
#             # Член с k_out_to_vp1 требует n[v] и n[vp]
#             # Член с k_out_to_vm1 требует n[v] и n[vp]
#
#             term = (k_in_from_vm1 * n_vm1 + k_in_from_vp1 * n_vp1 -
#                     (k_out_to_vp1 + k_out_to_vm1) * n[v]) * n[vp]
#
#             R += term
#
#         vv[v] = R * scale
#
#     return vv

# -------------------------------------------------------------
# 4. Правая часть для solve_ivp (диффузия + VV)
# -------------------------------------------------------------
def rhs(t, y, r, D, v_max, T, coeffs_dict, bc_outer="dirichlet"):
    nv = v_max + 1
    nr = len(r)
    N = y.reshape(nv, nr)  # [v, i]
    if bc_outer == "dirichlet":
        N = N.copy()
        N[:, -1] = 0.0
    dNdt = np.zeros_like(N)

    # VV-обмен (локально по каждому радиальному узлу)
    for ir in range(nr):
        n_local = N[:, ir]
        dNdt[:, ir] = vv_term(n_local, v_max, T, coeffs_dict)

    # Диффузия по радиусу для каждого вибрационного уровня
    for v in range(nv):
        diff_profile = diffusion_term(N[v, :], r, D, bc_outer=bc_outer)
        dNdt[v, :] += diff_profile

    if bc_outer == "dirichlet":
        dNdt[:, -1] = 0.0

    return dNdt.flatten()


# -------------------------------------------------------------
# 5. Начальные условия (гауссов пучок, заселение v=1)
# -------------------------------------------------------------
def initial_condition(r, v_max, frac_excited=0.33, w0=46e-4, total_density=5.14e18):
    """
    frac_excited: доля молекул в v=1 на оси пучка
    w0: радиус пучка (см) – в статье ~80 мкм = 80e-4 см
    """
    nr = len(r)
    nv = v_max + 1
    y0 = np.zeros((nv, nr))

    n1_rel = frac_excited
    n0_rel = 1.0 - n1_rel
    profile = np.exp(-2 * r**2 / w0**2)

    for ir in range(nr):
        y0[0, ir] = n0_rel * profile[ir] * total_density
        if nv > 1:
            y0[1, ir] = n1_rel * profile[ir] * total_density

    return y0.flatten()


def probe_weighted_population(N_sol, r, sigma_probe):
    weight = np.exp(-2.0 * (r / sigma_probe) ** 2)
    denom = np.trapz(weight * r, r)
    nv, _, nt = N_sol.shape
    obs = np.zeros((nv, nt))
    for v in range(nv):
        for it in range(nt):
            num = np.trapz(N_sol[v, :, it] * weight * r, r)
            obs[v, it] = num / denom if denom > 0 else 0.0
    return obs


def print_diagnostics(sol, N_sol, r, title=""):
    nv, _, nt = N_sol.shape
    dr = r[1] - r[0]
    min_all = float(np.min(N_sol))
    print(f"\n=== Diagnostics {title} ===")
    print(f"min(N): {min_all:.6e}")

    neg_mask = N_sol < 0
    if np.any(neg_mask):
        first_it = int(np.where(np.any(neg_mask, axis=(0, 1)))[0][0])
        print(f"first negative at t = {sol.t[first_it] * 1e6:.4f} us")
    else:
        print("negative populations: none")

    total = np.zeros(nt)
    for it in range(nt):
        total[it] = np.sum([2 * np.pi * np.sum(N_sol[v, :, it] * r * dr) for v in range(nv)])
    print(f"total integrated N: min={np.min(total):.6e}, max={np.max(total):.6e}")


def run_diffusion_only_check(y0, r, D, v_max, T, bc_outer, t_max=3e-6):
    """
    Диагностика: VV выключен (coeffs_dict={}).
    Показывает стабильность и баланс интегральной плотности для чистой диффузии.
    """
    t_eval = np.linspace(0, t_max, 80)
    sol_d = solve_ivp(
        rhs,
        (0, t_max),
        y0,
        t_eval=t_eval,
        method='Radau',
        args=(r, D, v_max, T, {}, bc_outer),
        rtol=1e-6,
        atol=1e-10,
    )
    N_sol_d = sol_d.y.reshape(v_max + 1, len(r), -1)
    print_diagnostics(sol_d, N_sol_d, r, title=f"[diffusion-only, bc={bc_outer}]")


# -------------------------------------------------------------
# 6. Усреднение по радиусу (сравнение с экспериментом)
# -------------------------------------------------------------
def radial_average(y, r, v_max, t_eval):
    """
    y: решение от solve_ivp (формат: (nr*nv, nt))
    возвращает матрицу <n_v>(t) размером (nt, nv)
    """
    nr = len(r)
    nv = v_max + 1
    nt = len(t_eval)

    # переупаковываем в (nr, nv, nt)
    y_reshaped = y.reshape((nr, nv, -1), order='F')  # Fortran-стиль, т.к. solve_ivp возвращает C

    avg = np.zeros((nt, nv))
    dr = r[1] - r[0]

    for it in range(nt):
        for v in range(nv):
            prof = y_reshaped[:, v, it]
            # интеграл ∫ n_v(r) * r dr
            integral = np.trapz(prof * r, r)
            # усреднение по площади πR^2
            avg[it, v] = 2 * integral / (r[-1] ** 2)
    return avg


# -------------------------------------------------------------
# 7. Параметры расчёта (из статьи)
# -------------------------------------------------------------
v_max = 5          # уровни 0..5
T = 300.0          # комнатная температура
D = 0.2            # коэффициент диффузии, см²/с (N2/O2 ~ при 1 бар)
BC_OUTER = "dirichlet"   # "dirichlet" or "neumann"
SIGMA_PUMP = 46e-4       # 46 um
SIGMA_PROBE = 46e-4      # 46 um (или 65e-4 как sigma_meas)
RUN_DIFFUSION_ONLY_CHECK = False

# # радиальная сетка
# R_max = 0.03       # см = 300 мкм (чуть больше пучка)
# nr = 51
# r = np.linspace(0, R_max, nr)

# радиальная сетка
R_max = 0.5       # см = 300 мкм (чуть больше пучка)
nr = 200
r = np.linspace(0, R_max, nr)

# временная сетка (как в эксперименте)
t_max = 10e-6
t_eval = np.linspace(0, t_max, 200)

# начальные условия
y0 = initial_condition(r, v_max, frac_excited=0.37, w0=SIGMA_PUMP, total_density=5.14e18)
if RUN_DIFFUSION_ONLY_CHECK:
    run_diffusion_only_check(y0, r, D, v_max, T, BC_OUTER)

# -------------------------------------------------------------
# 8. Решение системы
# -------------------------------------------------------------
print("Решаем систему с диффузией...")
sol = solve_ivp(
    rhs,
    (0, t_max),
    y0,
    t_eval=t_eval,
    method='Radau',
    args=(r, D, v_max, T, COEFFS, BC_OUTER),
    rtol=1e-6,
    atol=1e-10
)

print(f"Успешно: {sol.success}")

# # Восстанавливаем решение и считаем fractional populations (как в ODE_O2)
# Nr = len(r)
# nv = v_max + 1
# N_sol = sol.y.reshape(v_max + 1, Nr, -1)
# dr = r[1] - r[0]
# full_N = np.array(
#     [[2 * np.pi * np.sum(N_sol[v, :, t] * r * dr) for t in range(len(sol.t))] for v in range(v_max + 1)]
# )
# obs_N = probe_weighted_population(N_sol, r, sigma_probe=SIGMA_PROBE)
#
# # ИСПРАВЛЕНИЕ: Используем абсолютные концентрации как в статье
# # В статье O2 (2005) калибровка по室温ному спектру, где все молекулы в v=0
# # Поэтому f(v) = N_v / N_total, где N_total - полная плотность O2
# total_density_initial = 5.14e18  # см⁻³, полная плотность O2 при 760 Torr, 300K
# fractional = obs_N / total_density_initial  # абсолютные доли (без перенормировки)
#
# print_diagnostics(sol, N_sol, r, title=f"(bc={BC_OUTER}, D={D})")
#
# # Диагностика: проверка суммы f(v) (должна быть <1 из-за диффузии)
# sum_f = np.sum(fractional, axis=0)
# print(f"\n=== Диагностика нормировки ===")
# print(f"Сумма f(v) при t=0: {sum_f[0]:.4f} (должна быть ~1.0)")
# print(f"Сумма f(v) при t=10 мкс: {sum_f[-1]:.4f} (должна быть <1 из-за диффузии)")
# print(f"Потеря сигнала за 10 мкс: {(1 - sum_f[-1]/sum_f[0])*100:.2f}%")

# Восстанавливаем решение и считаем fractional populations
Nr = len(r)
nv = v_max + 1
N_sol = sol.y.reshape(v_max + 1, Nr, -1)
dr = r[1] - r[0]
full_N = np.array(
    [[2 * np.pi * np.sum(N_sol[v, :, t] * r * dr) for t in range(len(sol.t))] for v in range(v_max + 1)]
)
obs_N = probe_weighted_population(N_sol, r, sigma_probe=SIGMA_PROBE)

# ========== ПРАВИЛЬНАЯ НОРМИРОВКА (калибровка по начальному сигналу) ==========
# В эксперименте: калибровка по室温ному спектру
# При t=0: f(0)_expected = 0.63, f(1)_expected = 0.37
f0_expected = 0.63
f1_expected = 0.37

# Калибровочный множитель (приводим obs_N к абсолютным долям)
# При t=0: obs_N[0,0] * scale = f0_expected
scale = f0_expected / obs_N[0, 0]

fractional = obs_N * scale

print_diagnostics(sol, N_sol, r, title=f"(bc={BC_OUTER}, D={D})")

# Диагностика
print(f"\n=== Диагностика нормировки (калибровка по теоретическим долям) ===")
print(f"Scale factor: {scale:.4f}")
print(f"f(0) при t=0: {fractional[0,0]:.4f} (цель: {f0_expected})")
print(f"f(1) при t=0: {fractional[1,0]:.4f} (цель: {f1_expected})")
print(f"Отношение f(0)/f(1): {fractional[0,0]/fractional[1,0]:.3f}")

sum_f = np.sum(fractional, axis=0)
print(f"Сумма f(v) при t=0: {sum_f[0]:.4f}")
print(f"Сумма f(v) при t=10 мкс: {sum_f[-1]:.4f}")
print(f"Потеря сигнала за 10 мкс: {(1 - sum_f[-1]/sum_f[0])*100:.2f}%")


def diffusion_impact_report(t, frac_with, frac_without, v_levels=[0, 1, 2, 3]):
    """
    Печатает таблицу влияния диффузии
    """
    print("\n" + "=" * 70)
    print("ВЛИЯНИЕ ДИФФУЗИИ НА FRACTIONAL POPULATIONS")
    print("=" * 70)
    print(f"{'t, мкс':^10}", end="")
    for v in v_levels:
        print(f"{f'v={v}':^15}", end="")
    print()
    print("-" * 70)

    # Несколько временных точек
    times_us = [1, 2, 5, 10]
    for t_us in times_us:
        idx = np.argmin(np.abs(t * 1e6 - t_us))
        print(f"{t_us:^10.1f}", end="")
        for v in v_levels:
            with_diff = frac_with[v, idx]
            without_diff = frac_without[v, idx]
            rel_diff = (with_diff - without_diff) / (without_diff + 1e-12) * 100
            print(f"{rel_diff:^+15.2f}%", end="")
        print()
    print("=" * 70)

    # Вывод рекомендации
    max_effect = np.max(np.abs(frac_with - frac_without) / (frac_without + 1e-12))
    if max_effect < 0.01:
        print("⚠️ Влияние диффузии < 1% — можно пренебречь для большинства задач")
    elif max_effect < 0.05:
        print("📊 Влияние диффузии 1-5% — желательно учитывать для точных сравнений")
    else:
        print("🔴 Влияние диффузии >5% — обязательно учитывать диффузию!")


# Если у вас уже есть решение sol и N_sol
# Создаем решение "без диффузии" на лету

def quick_diffusion_check(sol, N_sol, r, t_eval, y0, v_max, T, coeffs_dict):
    """
    Быстрая проверка: пересчитывает без диффузии и сравнивает
    """
    print("Пересчитываем без диффузии для сравнения...")

    nr = len(r)
    nv = v_max + 1

    # Решение без диффузии
    sol_no_diff = solve_ivp(
        rhs, (0, t_max), y0, t_eval=t_eval,
        method='Radau', args=(r, 0.0, v_max, T, coeffs_dict, "dirichlet"),
        rtol=1e-6, atol=1e-10
    )

    N_no_diff = sol_no_diff.y.reshape(nv, nr, -1)

    # Усреднение
    obs_with = probe_weighted_population(N_sol, r, SIGMA_PROBE)
    obs_without = probe_weighted_population(N_no_diff, r, SIGMA_PROBE)

    # Используем абсолютные концентрации для обоих решений
    total_density_initial = 5.14e18
    frac_with = obs_with / total_density_initial
    frac_without = obs_without / total_density_initial

    # Отчет
    diffusion_impact_report(sol.t, frac_with, frac_without)

    return frac_with, frac_without


# Использование (после вашего основного расчета):
frac_with, frac_without = quick_diffusion_check(sol, N_sol, r, t_eval, y0, v_max, T, COEFFS)


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

# Добавляем график суммы f(v) для диагностики диффузии
fig2, ax2 = plt.subplots(1, 1, figsize=(8, 5))
ax2.plot(sol.t * 1e6, sum_f, 'b-', linewidth=2, label='Сумма f(v)')
ax2.axhline(y=1.0, color='r', linestyle='--', label='Начальное значение (1.0)')
ax2.set_xlabel('Время, мкс')
ax2.set_ylabel('Сумма f(v)')
ax2.set_title('Потеря сигнала из-за диффузии')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0.9, 1.05)

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
        # Процесс: O₂(v) + O₂(0) → O₂(v-1) + O₂(1)  → i1=v, f1=v-1, i2=0, f2=1
        key = f"{v}_{v - 1}_{0}_{1}"

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