import numpy as np

from constants import h, c


class Particle:
    """Molecular particle with vibrational energy levels.

    Parameters
    ----------
    name : str
        Molecular formula, e.g. "N2".
    mass : float
        Molecular mass [kg].
    diameter : float
        Kinetic diameter [m].
    num_elex_levels : int
        Number of electronic states.
    num_vibr_levels : list[int]
        Number of vibrational levels per electronic state.
    we : list[float]
        Harmonic vibrational constant per electronic state [1/m].
    wexe : list[float]
        Anharmonicity constant per electronic state [1/m].
    weye : list[float]
        Third-order anharmonicity constant per electronic state [1/m].

    Attributes
    ----------
    ev_i : list[np.ndarray]
        Vibrational energy levels relative to ground state v=0 [J].
    """
    def __init__(
        self,
        name: str,
        mass: float,
        diameter: float,
        num_elex_levels: int,
        num_vibr_levels: list[int],
        we: list[float],
        wexe: list[float],
        weye: list[float],
    ):
        self.name = name
        self.mass = mass
        self.diameter = diameter
        self.num_elex_levels = num_elex_levels
        self.num_vibr_levels = num_vibr_levels
        self.we = we
        self.wexe = wexe
        self.weye = weye

        self.ev_i: list[np.ndarray] = []
        for elvl in range(num_elex_levels):
            e_i = self.vibrational_energies(elvl)
            self.ev_i.append(e_i - e_i[0])

    def vibrational_energies(self, elvl: int, model: str = "anharmonic") -> np.ndarray:
        v = np.arange(self.num_vibr_levels[elvl]) + 0.5
        we = self.we[elvl]
        wexe = self.wexe[elvl]
        weye = self.weye[elvl]

        if model == "harmonic":
            e_i = we * v
        elif model == "anharmonic":
            e_i = we * v - wexe * v**2
        elif model == "anharmonic_higher":
            e_i = we * v - wexe * v**2 + weye * v**3
        else:
            raise ValueError(f"Unknown model: {model!r}")

        return np.array(e_i) * h * c
