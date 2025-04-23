import matplotlib.pyplot as plt
import pandas as pd

data_1 = pd.read_csv('k_vt_40_39.csv', sep=';', header=None)
data_2 = pd.read_csv('matlab_kvt.csv', sep=' ', header=None)

fig, ax = plt.subplots()

fho_comp_s2, = ax.plot(data_1[0], data_1[1], '-s')
fho_comp_s2.set_label('FHO-FR_comp_1')

fho_comp_s3, = ax.plot(data_2[0], data_2[1]*1e6, '-s')
fho_comp_s3.set_label('FHO-FR_matlab_2')

ax.set_yscale('log')
plt.legend(frameon=False, framealpha=0, fontsize='small')

# ax.set_xlim(0.2e-7, 1e-2)
# ax.set_ylim(1e-9, 1)
ax.set_xlabel(r'$\mathrm{T}$')
ax.set_ylabel(r'$\mathrm{k_{vt}}$')
plt.show()