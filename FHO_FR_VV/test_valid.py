from p_vv_mm import p_vv_int
from k_vv_mm import k_vv_mm
from particles_data import *
from constants import *
from approx_k_vv import *


from FHO_FR_VV.particles_data import *
from k_vv_mm import k_vv_mm

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd




# with open('val_coefs_O2.json', 'r') as f:
#     COEFFS = json.load(f)
#
# def k_lookup(coeffs_dict, v1, v2, v1p, v2p, T):
#
#     key = f"{v1}_{v2}_{v1p}_{v2p}"
#     coeffs = coeffs_dict.get(key)
#
#     if coeffs is None:
#         return 0.0
#
#     x = T**(-1/3)
#     return np.exp(sum(c * x**i for i, c in enumerate(coeffs)))



BASE_DIR = Path('.').resolve()
PARTICLE_A, PARTICLE_B = "O2", "O2"
# REGRESSION_COEFFS_PATH = BASE_DIR / f"{PARTICLE_A}_{PARTICLE_B}_regression_coefs12312.csv"
# REGRESSION_COEFFS_PATH = BASE_DIR / f"{PARTICLE_A}_{PARTICLE_B}_regression_coefs.csv"
REGRESSION_COEFFS_PATH = BASE_DIR / "regression_coefficients" / f"{PARTICLE_A}_{PARTICLE_B}_regression_coefs.csv"
FIGURES_PATH = BASE_DIR / "figures" / f"{PARTICLE_A}_{PARTICLE_B}_coefs_comparison.png"

def coeffs_dict_from_regression_csv(path: Path):
    df = pd.read_csv(path)
    a_cols = sorted(
        (c for c in df.columns if c.startswith("a_")),
        key=lambda s: int(s.split("_")[1]),
    )
    out = {}
    for _, row in df.iterrows():
        key = "_".join(str(int(row[c])) for c in ("i1", "f1", "i2", "f2"))
        coeffs = []
        for cname in a_cols:
            val = row[cname]
            if pd.isna(val):
                continue
            coeffs.append(float(val))
        out[key] = coeffs
    return out

if REGRESSION_COEFFS_PATH.exists():
    COEFFS = coeffs_dict_from_regression_csv(REGRESSION_COEFFS_PATH)
    print(f"k_VV: {REGRESSION_COEFFS_PATH.name} ({len(COEFFS)} переходов)")
else:
    with open(BASE_DIR / "coeffs_sparse2.json", "r") as f:
        COEFFS = json.load(f)
    print("k_VV: fallback coeffs_sparse2.json")

def k_lookup(coeffs_dict, i1, f1, i2, f2, T, rate_multiplier=1.0):
    key = f"{i1}_{f1}_{i2}_{f2}"
    coeffs = coeffs_dict.get(key)
    if coeffs is None:
        return 0.0
    x = T ** (-1/3)
    return rate_multiplier * np.exp(sum(c * x ** i for i, c in enumerate(coeffs)))

s2 = 2
s3 = 3
E = 1000000  # 1/m
# e_in_J = 1.98e-23 # E: 1/cm --> J (h * c * 100)
e_in_J = h * c * 100


# Параметры для O₂ (оценочно)
omega_O2 = 1580.0 * 2 * np.pi * 1e10  # 1580 cm^-1
m_reduced_O2 = 16.0 * 1.660539e-27     # кг
alpha_O2 = 4.3e10                      # м^-1
T = 300.0



# Для нерезонансного перехода Delta_E ≠ 0
# Например, (1,0) -> (0,1) с учётом ангармоничности
# Энергия первого кванта ~ 1580 см^-1 → 1580 * 1.986e-23 = 3.14e-20 J
Delta_E = 0.0  # или вычисли из твоих ev_i




data_1 = pd.read_csv('fho.csv', sep=';', header=None)
data_2 = pd.read_csv('billing.csv', sep=';', header=None)
data_3 = pd.read_csv('work.csv', sep=';', header=None)

data_reg = pd.read_csv('reg.csv', sep=';', header=None)

data_matlab_3000 = pd.read_csv('matlab_3000.csv', sep=';', header=None)
#data_matlab_s3 = pd.read_csv('matlab_s3.csv', sep=';', header=None)

x1 =[]
y_1 = []
y_reg = []
y_reg1 = []
y_reg2 = []
y_reg3 = []
y_reg4 = []
k_rates = []

for i in range(1, 12, 1):
    x1.append(i)
    print('x1:', i)
    y = k_vv_mm(O2, O2, i, i-1, 0, 1, 300)

    k_rate = k_vv_adamovich(i, 0, T, omega_O2, m_reduced_O2, alpha_O2, Delta_E)

    y_reg_ = k_lookup(COEFFS, i, i-1, 0, 1, 300)
    y_reg1_ = 1.8 * k_lookup(COEFFS, i, i - 1, 0, 1, 300)
    y_reg2_ = 2 * k_lookup(COEFFS, i, i - 1, 0, 1, 300)
    y_reg3_ = 2.5 * k_lookup(COEFFS, i, i - 1, 0, 1, 300)
    y_reg4_ = 2.7 * k_lookup(COEFFS, i, i - 1, 0, 1, 300)

    y_1.append(y)
    print('y1:', y)
    k_rates.append(k_rate)
    print('y1:', k_rate)

    y_reg.append(y_reg_)
    y_reg1.append(y_reg1_)
    y_reg2.append(y_reg2_)
    y_reg3.append(y_reg3_)
    y_reg4.append(y_reg4_)

    print('yreg:', y_reg)



fig, ax = plt.subplots()

fho_s2, = ax.plot(x1, y_1, '-o', markersize=5)
fho_s2.set_label('FHO-FR, code (i, i-1 --> 0, 1), T = 300K')

# fho_s_reg, = ax.plot(x1, k_rates, '-^', markersize=5)
# fho_s_reg.set_label('FHO-FR, approx (i, i-1 --> 0, 1), T = 300K')

fho_s_reg, = ax.plot(x1, y_reg, '-^', markersize=5)
fho_s_reg.set_label('FHO-FR, regression (i, i-1 --> 0, 1), T = 300K')

# fho_s_reg, = ax.plot(x1, y_reg1, '--', markersize=5)
# fho_s_reg.set_label('FHO-FR, 1.8*regression (i, i-1 --> 0, 1), T = 300K')
# fho_s_reg, = ax.plot(x1, y_reg2, '--', markersize=5)
# fho_s_reg.set_label('FHO-FR, 2*regression (i, i-1 --> 0, 1), T = 300K')
# fho_s_reg, = ax.plot(x1, y_reg3, '--', markersize=5)
# fho_s_reg.set_label('FHO-FR, 2.5*regression (i, i-1 --> 0, 1), T = 300K')
fho_s_reg, = ax.plot(x1, y_reg4, '--', markersize=5)
fho_s_reg.set_label('FHO-FR, 2.7 * regression (i, i-1 --> 0, 1), T = 300K')



fho_comp_s2, = ax.plot(data_1[0], data_1[1], '-*')
fho_comp_s2.set_label('FHO-FR, article (i, i-1 --> 0, 1), T = 300K')

fho_comp_s3, = ax.plot(data_2[0], data_2[1], '-*')
fho_comp_s3.set_label('Billing, article (i, i-1 --> 0, 1), T = 300K')

fho_comp_s4, = ax.plot(data_3[0], data_3[1], '-s')
fho_comp_s4.set_label('Ahn, Adamovich, article (i, i-1 --> 0, 1), T = 300K')

# iii = [i for i in range(0, 40)]
# fho_reg, = ax.plot(iii, data_reg, '-s')
# fho_reg.set_label('FHOreg')

# fho_matlab_3000, = ax.plot(data_matlab_3000[0], data_matlab_3000[1], '-^')


# fho_matlab_3000.set_label('FHO-FR_matlab_s2')

# fho_matlab_s3, = ax.plot(data_matlab_s3[0], data_matlab_s3[1], '-^')
# fho_matlab_s3.set_label('FHO-FR_matlab_s3')

ax.set_yscale('log')
plt.legend(frameon=True, framealpha=0.5, fontsize='x-small')

# ax.set_xlim(0.2e-7, 1e-2)
# ax.set_ylim(1e-9, 1)
ax.set_xlabel(r'$\mathrm{i, \ vibrational \ quantum \ number}$')
ax.set_ylabel(r'$\mathrm{k_{vv} (T), \ cm^3/s}$')

plt.grid(True, linestyle='--')

fig.savefig(FIGURES_PATH, dpi=600)
plt.show()