"""
k_VV — вычисление коэффициентов скорости V-V переходов.

Модули:
  k_VV          — численное интегрирование (matrix method)
  k_VV_regr     — регрессионные коэффициенты
  k_VV_fortran  — аналитическая формула FHO-FR
"""
from k_VV.k_VV import k_vv_mm, k_vv_mm_for_v
from k_VV.k_VV_regr import k_lookup, load_coeffs, load_coeffs_from_path, coeffs_dict_from_regression_csv
from k_VV.k_VV_fortran import k_vv_fhofer

__all__ = [
    "k_vv_mm", "k_vv_mm_for_v",
    "k_lookup", "load_coeffs", "load_coeffs_from_path", "coeffs_dict_from_regression_csv",
    "k_vv_fhofer",
]
