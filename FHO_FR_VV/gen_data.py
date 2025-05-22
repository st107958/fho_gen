from p_vv_mm import p_vv_int
from k_vv_mm import k_vv_mm
from particles_data import *
from constants import *

import matplotlib.pyplot as plt
import pandas as pd

E = 1000000  # 1/m
# e_in_J = 1.98e-23 # E: 1/cm --> J (h * c * 100)
e_in_J = h * c * 100


n_points = 24
T_min, T_max = 300, 30000  # Диапазон 300–32000 K
T_data = np.linspace(T_min**(-1/3), T_max**(-1/3), n_points) ** (-3)


K_data = []
for i in range(len(T_data)):
    k = k_vv_mm(N2, N2, 1, 0, 0, 1, T_data[i])
    K_data.append(k)
    print('k:', k)

df = pd.DataFrame({'T': T_data})
df.to_csv('T.csv', index=False, header=False)

df = pd.DataFrame({'K': K_data})
df.to_csv('K.csv', index=False, header=False)

# with open('data/K.csv', 'w', encoding='utf-8') as file:
#     for num in K_data:
#         file.write(f"{num}\n")
#
# with open('data/T.csv', 'w', encoding='utf-8') as file:
#     for num in T_data:
#         file.write(f"{num}\n")