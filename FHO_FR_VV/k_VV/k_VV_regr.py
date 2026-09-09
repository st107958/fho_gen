"""
Регрессионные коэффициенты k_VV: загрузка и полиномиальная оценка.

k_VV(T) = exp(sum(a_j * T^(-j/3)))
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def coeffs_dict_from_regression_csv(path: Path) -> dict[str, list[float]]:
    """Парсит CSV с колонками i1,f1,i2,f2,a_0,...,a_n → dict["i1_f1_i2_f2"] = [a_j]."""
    df = pd.read_csv(path)
    a_cols = sorted(
        (c for c in df.columns if c.startswith("a_")),
        key=lambda s: int(s.split("_")[1]),
    )
    out: dict[str, list[float]] = {}
    for _, row in df.iterrows():
        key = "_".join(str(int(row[c])) for c in ("i1", "f1", "i2", "f2"))
        coeffs = [float(row[cname]) for cname in a_cols if not pd.isna(row[cname])]
        out[key] = coeffs
    return out


def load_coeffs(base_dir: Path, particle_a: str, particle_b: str) -> dict[str, list[float]]:
    """Загружает коэффициенты регрессии из CSV (или fallback JSON)."""
    path = base_dir / f"{particle_a}_{particle_b}_regression_coefs.csv"
    if path.exists():
        coeffs = coeffs_dict_from_regression_csv(path)
        print(f"k_VV: {path.name} ({len(coeffs)} переходов)")
        return coeffs
    with open(base_dir / "coeffs_sparse2.json", "r", encoding="utf-8") as f:
        print("k_VV: fallback coeffs_sparse2.json")
        return json.load(f)


def load_coeffs_from_path(path: str | Path) -> dict[str, list[float]]:
    """Загружает коэффициенты регрессии из указанного CSV-файла."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"CSV с коэффициентами не найден: {p}")
    coeffs = coeffs_dict_from_regression_csv(p)
    print(f"k_VV: {p.name} ({len(coeffs)} переходов)")
    return coeffs


def k_lookup(
    coeffs_dict: dict,
    i1: int,
    f1: int,
    i2: int,
    f2: int,
    T: float,
    rate_multiplier: float = 1.0,
) -> float:
    """Полиномиальная оценка k_VV по регрессионным коэффициентам."""
    coeffs = coeffs_dict.get(f"{i1}_{f1}_{i2}_{f2}")
    if coeffs is None:
        return 0.0
    x = T ** (-1 / 3)
    return rate_multiplier * float(np.exp(sum(c * x**i for i, c in enumerate(coeffs))))
