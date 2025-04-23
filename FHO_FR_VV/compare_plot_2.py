import matplotlib.pyplot as plt
import pandas as pd

data_1 = pd.read_csv('pvt.csv', sep=';', header=None)
data_2 = pd.read_csv('matlab_pvt.csv', sep=';', header=None)

fig, ax = plt.subplots()

fho_comp_s2, = ax.plot(data_1[0], data_1[1], '-s')
fho_comp_s2.set_label('FHO-FR_comp_1')

fho_comp_s3, = ax.plot(data_2[0], data_2[1], '-s')
fho_comp_s3.set_label('FHO-FR_matlab_2')

ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(frameon=False, framealpha=0, fontsize='small')

# ax.set_xlim(0.2e-7, 1e-2)
# ax.set_ylim(1e-9, 1)
ax.set_xlabel(r'$\mathrm{E}$')
ax.set_ylabel(r'$\mathrm{p_{vt}}$')
plt.show()