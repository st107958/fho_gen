from FHO_FR_VV.p_vv_mm import p_vv_int
from FHO_FR_VV.k_vv_mm import k_vv_mm
from FHO_FR_VV.particles_data import *
from FHO_FR_VV.constants import *

import matplotlib.pyplot as plt
import pandas as pd

import time

E = 1000000  # 1/m
e_in_J = h * c * 100


n_points = 20
T_min, T_max = 1000, 15000  # Диапазон
# T_data = np.linspace(T_min**(-1/3), T_max**(-1/3), n_points) ** (-3)
T_data = np.linspace(T_min, T_max, n_points)

# T_data = np.random.uniform(T_min, T_max, n_points)
T_data = np.sort(T_data)  # если нужно отсортировать


df_k = pd.DataFrame()

start_time = time.time()

column_values = []

for i in range(1, 10):
    k = k_vv_mm(N2, N2, i-1, i, 41, 40,  T_data)
    column_values.append(k)

df_k['k'] = column_values



# df_long = pd.melt(
#     df_k,
#     id_vars=['Temperature'],
#     value_vars=[f'j_{j}' for j in range(1, 2)],
#     var_name='i',  # название для столбца с j
#     value_name='kVV'
# )

# # Преобразуем 'i' из 'j_1', 'j_2', ... в числовые значения
# df_long['i'] = df_long['i'].str.replace('j_', '').astype(int)
#
# # Переименовываем столбцы в нужные имена
# df_long = df_long.rename(columns={'Temperature': 'T'})
#
# # Сортируем для удобства
# df_long = df_long.sort_values(['i', 'T']).reset_index(drop=True)
#
# print("Длинный формат:")
# print(df_long.head())
# print(f"\nРазмер: {df_long.shape}")
#
# # Сохраняем длинный формат
# df_long.to_csv('O2_v1.csv', index=False)



end_time = time.time()
execution_time = end_time - start_time
print(f"Время выполнения: {execution_time:.4f} секунд")



# print(df_k.head())

df_k.to_csv('N2.csv', index=False, header=False)

