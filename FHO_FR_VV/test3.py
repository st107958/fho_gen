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


data_s2 = pd.read_csv('fho_s2.csv', sep=';', header=None)
data_s3 = pd.read_csv('fho_s3.csv', sep=';', header=None)
data_matlab_s2 = pd.read_csv('matlab_s2.csv', sep=';', header=None)
data_matlab_s3 = pd.read_csv('matlab_s3.csv', sep=';', header=None)

x1 =[]
y_s2 = []
for i in range(2, 41):
    x1.append(i)
    print('x1:', i)
    y = p_vv_int(N2, N2, i, i-s2, 0, s2, 10000, 'trapez')
    y_s2.append(y)
    print('y1:', y)

x2 =[]
y_s3 = []
for i in range(3, 42):
    x2.append(i)
    print('x2:', i)
    y = p_vv_int(N2, N2, i, i-s3, 0, s3, 10000, 'trapez')
    y_s3.append(y)
    print('y2:', y)


# x1 =[]
# y_s2 = []
# for i in range(2, 5):
#     x1.append(i)
#     print('x1:', i)
#     y = p_vv_int(N2, N2, i, i-s2, 0, s2, 10000, 'trapez')
#     y_s2.append(y)
#     print('y1:', y)
#
# x2 =[]
# y_s3 = []
# for i in range(3, 5):
#     x2.append(i)
#     print('x2:', i)
#     y = p_vv_int(N2, N2, i, i-s3, 0, s3, 10000, 'trapez')
#     y_s3.append(y)
#     print('y2:', y)

fig, ax = plt.subplots()

fho_s2, = ax.plot(x1, y_s2, '-^', markersize=5)
fho_s2.set_label(r'$\mathrm{P_{VV}}((i, 0 \to i - s, s);\ s=2;\ E=10^{-4}\ \mathrm{cm}^{-1})$, numerical calculation')

fho_s3, = ax.plot(x2, y_s3, '-p', markersize=5)
fho_s3.set_label(r'$\mathrm{P_{VV}}((i, 0 \to i - s, s);\ s=3;\ E=10^{-4}\ \mathrm{cm}^{-1})$, numerical calculation')

fho_comp_s2, = ax.plot(data_s2[0], data_s2[1], '--')
fho_comp_s2.set_label(r'$\mathrm{P_{VV}}((i, 0 \to i - s, s);\ s=2;\ E=10^{-4}\ \mathrm{cm}^{-1})$, analytic model')

fho_comp_s3, = ax.plot(data_s3[0], data_s3[1], '--')
fho_comp_s3.set_label(r'$\mathrm{P_{VV}}((i, 0 \to i - s, s);\ s=3;\ E=10^{-4}\ \mathrm{cm}^{-1})$, analytic model')

# fho_matlab_s2, = ax.plot(data_matlab_s2[0], data_matlab_s2[1], '-^')
# fho_matlab_s2.set_label('FHO-FR_matlab_s2')
#
# fho_matlab_s3, = ax.plot(data_matlab_s3[0], data_matlab_s3[1], '-^')
# fho_matlab_s3.set_label('FHO-FR_matlab_s3')

ax.set_yscale('log')
plt.legend(frameon=True, framealpha=0.5, fontsize='x-small')

# ax.set_xlim(0.2e-7, 1e-2)
# ax.set_ylim(1e-9, 1)
ax.set_xlabel(r'$\mathrm{Vibrational \ quantum \ number, i}$')
ax.set_ylabel(r'Transition probability, $\mathrm{P_{VV}}$')
plt.grid(True, linestyle='--')

fig.savefig('plot_VV_prob.png', dpi=600)

plt.show()