from __future__ import annotations

import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from particles_data import N2
from k_vv_mm import k_vv_mm

# Имена — только для имён файлов (объект O2 задаёт данные модели)
PARTICLE1_NAME = "N2"
PARTICLE2_NAME = "N2"

M1, M2 = N2, N2
# Максимальные колебательные уровни 0 … V_MAX включительно (как в R_VV / coeffs_sparse)
V_MAX = 6

# Сетка температур для датасета и регрессии [K]
TEMPERATURES_K = np.linspace(300.0, 6000.0, 10)

WORKDIR = Path(".")
DATASET_PATH = WORKDIR / f"{PARTICLE1_NAME}_{PARTICLE2_NAME}_dataset.csv"
# Переходы, которые вызывает k_lookup внутри R_VV_fast (все четыре ветки, v и v_ по 0…V_MAX)
DATASET_RVV_FAST_PATH = (
    WORKDIR / f"{PARTICLE1_NAME}_{PARTICLE2_NAME}_dataset_rvv_fast.csv"
)
COEFS_PATH = WORKDIR / f"{PARTICLE1_NAME}_{PARTICLE2_NAME}_regression_coefs.csv"

# Подбор степени полинома для ln k (число коэффицициентов = степень + 1)
POLY_TERM_COUNTS = range(4, 8)

# Параллельная генерация dataset_rvv_fast (процессы); None = все доступные ядра
RVV_FAST_WORKERS = os.cpu_count()

# Переход для визуальной проверки графика (начальное → конечное по молекулам 1 и 2)


def transitions_for_r_vv_fast(v_max: int) -> set[tuple[int, int, int, int]]:
    """Уникальные (i1, f1, i2, f2), как в R_VV_fast при переборе v, v_ по 0…v_max."""
    out: set[tuple[int, int, int, int]] = set()
    for v in range(v_max + 1):
        for v_ in range(v_max + 1):
            if v + 1 <= v_max and v_ + 1 <= v_max:
                out.add((v_, v_ + 1, v + 1, v))
            if v - 1 >= 0 and v_ - 1 >= 0:
                out.add((v_, v_ - 1, v - 1, v))
            if v - 1 >= 0 and v_ + 1 <= v_max:
                out.add((v_, v_ + 1, v, v - 1))
            if v + 1 <= v_max and v_ - 1 >= 0:
                out.add((v_, v_ - 1, v, v + 1))
    return out


_M1_RVV = None
_M2_RVV = None


def _rvv_worker_init(m1, m2) -> None:
    global _M1_RVV, _M2_RVV
    _M1_RVV, _M2_RVV = m1, m2


def _compute_rvv_row(task: tuple[int, int, int, int, float]) -> dict:
    i1, f1, i2, f2, Tf = task
    from k_vv_mm import k_vv_mm

    k = float(k_vv_mm(_M1_RVV, _M2_RVV, i1, f1, i2, f2, Tf))
    return {"T": float(Tf), "i1": i1, "f1": f1, "i2": i2, "f2": f2, "k_VV": k}


quad_rvv = sorted(transitions_for_r_vv_fast(V_MAX))
tasks_rvv = [
    (i1, f1, i2, f2, float(T))
    for T in TEMPERATURES_K
    for i1, f1, i2, f2 in quad_rvv
]


def _rvv_executor(max_workers: int | None, m1, m2):
    exec_kw = {"initializer": _rvv_worker_init, "initargs": (m1, m2)}
    return ProcessPoolExecutor(max_workers=max_workers, **exec_kw)


if __name__ == "__main__":
    print(f"Доступно ядер CPU: {os.cpu_count()}")
    print(f"Запрошено workers: {RVV_FAST_WORKERS}")
    print(f"Максимально возможное: {mp.cpu_count()}")
    t0_rvv = time.time()
    rows_rvv: list[dict] = []
    with _rvv_executor(RVV_FAST_WORKERS, M1, M2) as ex:
        fmap = {ex.submit(_compute_rvv_row, t): t for t in tasks_rvv}
        for fut in as_completed(fmap):
            row = fut.result()
            i1, f1, i2, f2 = row["i1"], row["f1"], row["i2"], row["f2"]
            print(
                f"готово: T={row['T']:.6g} K, i1={i1} f1={f1} i2={i2} f2={f2}, k_VV={row['k_VV']:.6g}",
                flush=True,
            )
            rows_rvv.append(row)

    df_dataset_rvv_fast = pd.DataFrame(rows_rvv)
    df_dataset_rvv_fast = df_dataset_rvv_fast.sort_values(
        ["T", "i1", "f1", "i2", "f2"], kind="mergesort"
    ).reset_index(drop=True)
    elapsed_rvv = time.time() - t0_rvv
    print(
        f"R_VV_fast: уникальных переходов {len(quad_rvv)}, строк {len(df_dataset_rvv_fast)}, "
        f"время {elapsed_rvv:.2f} s, workers={RVV_FAST_WORKERS!r}"
    )
    # df_dataset_rvv_fast.to_csv(DATASET_RVV_FAST_PATH, index=False)
    # print(f"Сохранено: {DATASET_RVV_FAST_PATH.resolve()}")
    df_dataset_rvv_fast.head()