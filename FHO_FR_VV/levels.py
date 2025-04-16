from constants import *
import numpy as np


def levels_e_ex(M, elvl, ind=2):
    e_i = []

    for i in range(0, M.num_vibr_levels[elvl]):

        if ind == 1:  # harmonic oscillator
            e_i.append((M.we[elvl]*(i + 0.5))*h*c)

        elif ind == 2:  # anharmonic oscillator, by default
            e_i.append((M.we[elvl]*(i + 0.5) - M.wexe[elvl]*np.power((i + 0.5), 2))*h*c)  # J

        elif ind == 3:  # advanced anharmonic oscillator
            e_i.append((M.we[elvl]*(i + 0.5) - M.wexe[elvl]*np.power((i + 0.5), 2) +
                        M.weye[elvl]*np.power((i + 0.5), 3))*h*c)

    return np.array(e_i)

