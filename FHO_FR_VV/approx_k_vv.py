import numpy as np
from scipy.constants import k, hbar


def k_vv_adamovich(v, w, T, omega, m, alpha, Delta_E):
    """
    V–V rate from Adamovich FHO-FR, formula (7) from the paper.

    Parameters
    ----------
    v, w : int
        Vibrational quantum numbers (transition v+1,w -> v,w+1)
    T : float
        Temperature (K)
    omega : float
        Oscillator frequency (rad/s) — average vibrational quantum
    m : float
        Collision reduced mass (kg)
    alpha : float
        Repulsive potential parameter (m^-1)
    Delta_E : float
        Vibrational energy defect (J)

    Returns
    -------
    k : float
        Rate coefficient (cm^3/s)
    """
    # theta' formula
    theta_prime = (4.0 * np.pi ** 2 * omega ** 2 * m) / (alpha ** 2 * k)

    # lambda formula
    lambda_val = (1.0 / (3.0 * np.sqrt(2.0))) * (theta_prime / T) * (abs(Delta_E) / (omega * hbar))

    # Prefactor
    prefactor = (1.0 / 16.0) * (alpha ** 2 * k * T) / (2.0 * omega ** 2 * m)

    # Quantum number factor
    quant_factor = (v + 1.0) * (w + 1.0)

    # Lambda exponential factor
    exp_factor = (3.0 - np.exp(-2.0 * lambda_val / 3.0)) * np.exp(-2.0 * lambda_val / 3.0)

    # Energy defect factor (exothermic direction)
    energy_factor = np.exp(Delta_E / (2.0 * k * T))

    # Final rate
    k_rate = prefactor * quant_factor * exp_factor * energy_factor

    return k_rate * 1e-6