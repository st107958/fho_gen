from FHO_FR_VV.p_vv_mm import p_vv_int
from FHO_FR_VV.k_vv_mm import k_vv_mm
from FHO_FR_VV.particles_data import *
from FHO_FR_VV.constants import *

import matplotlib.pyplot as plt
import pandas as pd

import time

E = 1000000  # 1/m
e_in_J = h * c * 100


n_points = 24
T_min, T_max = 500, 16000  # Диапазон
T_data = np.linspace(T_min**(-1/3), T_max**(-1/3), n_points) ** (-3)


df_k = pd.DataFrame()

start_time = time.time()

for j in range(1, 10):
    column_values = []
    for i in range(len(T_data)):
        k = k_vv_mm(N2, N2, 41, 40, j - 1, j, T_data[i])
        column_values.append(k)

    df_k[f'j_{j}'] = column_values

end_time = time.time()
execution_time = end_time - start_time
print(f"Время выполнения: {execution_time:.4f} секунд")

df_k.insert(0, 'Temperature', T_data)

# print(df_k.head())

df_k.to_csv('N2.csv', index=False, header=False)

