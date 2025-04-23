from p_vv_mm import p_vv_int
from k_vv_mm import k_vv_mm
from particles_data import *
from constants import *

import matplotlib.pyplot as plt
import pandas as pd
import time

s2 = 2
s3 = 3
E = 1000000  # 1/m
# e_in_J = 1.98e-23 # E: 1/cm --> J (h * c * 100)
e_in_J = h * c * 100

#data_matlab1 = pd.read_csv('matlab_kvv.csv', sep=';', header=None)

x1 =[]
y_1 = []
for i in range(4, 40):
    start = time.perf_counter()  # Более точный таймер

    x1.append(i)
    print('x1:', i)
    y = k_vv_mm(N2, N2, 41, 40, i - 1, i, 1000) #* 1e6
    y_1.append(y)
    print('y1:', y)
    end = time.perf_counter()
    print(f"Время выполнения: {end - start:.6f} секунд")

fig, ax = plt.subplots()

fho_1, = ax.plot(x1, y_1, '-P')
fho_1.set_label('FHO-FR_x1')

#data_m_1, = ax.plot(data_matlab1[0], data_matlab1[1], '-s')
#data_m_1.set_label('FHO-FR_comp_2')

ax.set_yscale('log')
plt.legend(frameon=False, framealpha=0, fontsize='small')

# ax.set_xlim(0.2e-7, 1e-2)
# ax.set_ylim(1e-9, 1)
ax.set_xlabel(r'$\mathrm{i}$')
ax.set_ylabel(r'$\mathrm{k_{vv}}$')
plt.show()