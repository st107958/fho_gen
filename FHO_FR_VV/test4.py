from p_vv_mm import p_vv_int
from particles_data import *
from constants import *

import matplotlib.pyplot as plt
import pandas as pd

s2 = 2
s3 = 3
E = 1000000  # 1/m
# e_in_J = 1.98e-23 # E: 1/cm --> J (h * c * 100)
e_in_J = h * c * 100


# data_s2 = pd.read_csv('fho_s2.csv', sep=';', header=None)
# data_s3 = pd.read_csv('fho_s3.csv', sep=';', header=None)

x1 =[]
y_s2 = []
for i in range(100, 10000, 10):
    x1.append(i)
    y_s2.append(p_vv_int(N2, N2, s2+1, 1, 0, s2, i, 'trapez'))

x2 =[]
y_s3 = []
for i in range(100, 10000, 100):
    x2.append(i)
    y_s3.append(p_vv_int(N2, N2, s3+1, 1, 0, s3, i, 'trapez'))

fig, ax = plt.subplots()

fho_s2, = ax.plot(x1, y_s2, '-P')
fho_s2.set_label('FHO-FR_s2')

fho_s3, = ax.plot(x2, y_s3, '-P')
fho_s3.set_label('FHO-FR_s3')

# fho_comp_s2, = ax.plot(data_s2[0], data_s2[1], '-s')
# fho_comp_s2.set_label('FHO-FR_comp_s2')
#
# fho_comp_s3, = ax.plot(data_s3[0], data_s3[1], '-s')
# fho_comp_s3.set_label('FHO-FR_comp_s3')

ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(frameon=False, framealpha=0, fontsize='small')

# ax.set_ylim(0.2e-7, 1e-2)
ax.set_xlabel(r'$\mathrm{E, 1/см}$')
ax.set_ylabel(r'$\mathrm{P}$')
plt.show()