from p_vv_mm import p_vv_int
from k_vv_mm import k_vv_mm
from particles_data import *
from constants import *

import matplotlib.pyplot as plt
import pandas as pd

s2 = 2
s3 = 3
E = 1000000  # 1/m
# e_in_J = 1.98e-23 # E: 1/cm --> J (h * c * 100)
e_in_J = h * c * 100


data_1 = pd.read_csv('3000k.csv', sep=';', header=None)
data_2 = pd.read_csv('300k.csv', sep=';', header=None)
data_3 = pd.read_csv('fho10_3000.csv', sep=';', header=None)

data_matlab_3000 = pd.read_csv('matlab_3000.csv', sep=';', header=None)
#data_matlab_s3 = pd.read_csv('matlab_s3.csv', sep=';', header=None)

x1 =[]
y_1 = []
for i in range(28, 42):
    x1.append(i)
    print('x1:', i)
    y = k_vv_mm(N2, N2, 41, 40, i-1, i, 3000)
    y_1.append(y)
    print('y1:', y)
#
# x2 =[]
# y_2 = []
# for i in range(10, 40, 4):
#     x2.append(i)
#     print('x2:', i)
#     y = k_vv_mm(N2, N2, 41, 40, i-1, i, 300)
#     y_2.append(y)
#     print('y2:', y)

x3 =[]
y_3 = []
for i in range(1, 40, 3):
    x3.append(i)
    print('x3:', i)
    y = k_vv_mm(N2, N2, 1, 0, i-1, i, 3000)
    y_3.append(y)
    print('y3:', y)

fig, ax = plt.subplots()

fho_s2, = ax.plot(x1, y_1, '-P')
fho_s2.set_label('FHO-FR_3000')
#
# fho_s3, = ax.plot(x2, y_2, '-P')
# fho_s3.set_label('FHO-FR_300')

fho_s4, = ax.plot(x3, y_3, '-P')
fho_s4.set_label('FHO-FR_10_3000')

fho_comp_s2, = ax.plot(data_1[0], data_1[1], '-s')
fho_comp_s2.set_label('FHO-FR_comp_1')

fho_comp_s3, = ax.plot(data_2[0], data_2[1], '-s')
fho_comp_s3.set_label('FHO-FR_comp_2')

fho_comp_s4, = ax.plot(data_3[0], data_3[1], '-s')
fho_comp_s4.set_label('FHO-FR_comp_3')



# fho_matlab_3000, = ax.plot(data_matlab_3000[0], data_matlab_3000[1], '-^')
# fho_matlab_3000.set_label('FHO-FR_matlab_s2')

# fho_matlab_s3, = ax.plot(data_matlab_s3[0], data_matlab_s3[1], '-^')
# fho_matlab_s3.set_label('FHO-FR_matlab_s3')

ax.set_yscale('log')
plt.legend(frameon=False, framealpha=0, fontsize='small')

# ax.set_xlim(0.2e-7, 1e-2)
# ax.set_ylim(1e-9, 1)
ax.set_xlabel(r'$\mathrm{i}$')
ax.set_ylabel(r'$\mathrm{k_{vv}}$')
plt.show()