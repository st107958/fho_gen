import math
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

c = 299792458
h = 6.6261e-34
k = 1.3806e-23
spectr_constants = {
    'N2': 235857,
    'O2': 158019,
    'NO': 190420
}
pa_to_atm = 9.86923e-6

molecule = 'N2' # тут менять
particle = 'N2' # тут

lines_to_skip = 1
def rel_times(temperature: np.array, coefficient: np.array, molecule: str):
    t_array = []
    for i in range(len(coefficient)):
        T = temperature[i]
        kvt = coefficient[i] / 1e6

        t = (k * T) / (kvt * (1 - math.exp(-(spectr_constants[particle] * c * h) / (k * T)))) * pa_to_atm

        # t_array.append([float(T**(-1/3)), float(t)])
        t_array.append([float(T), float(t)])

    return t_array


T_data = np.loadtxt('T1.csv')
k_data1 = np.loadtxt('K1.csv')
k_data2 = np.loadtxt('K2.csv')
pt_T_array_fho1 = np.array(rel_times(T_data, k_data1, molecule))
pt_T_array_fho2 = np.array(rel_times(T_data, k_data2, molecule))



filename = 'FHO_FR1.csv' #тут
np.savetxt(filename, pt_T_array_fho1)

filename = 'FHO_FR2.csv' #тут
np.savetxt(filename, pt_T_array_fho2)

