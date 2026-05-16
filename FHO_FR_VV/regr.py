from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import warnings

# ============== КОНФИГ ==============
PARTICLE1_NAME = "N2"
PARTICLE2_NAME = "N2"
WORKDIR = Path(".")
DATASET_RVV_FAST_PATH = WORKDIR / f"{PARTICLE1_NAME}_{PARTICLE2_NAME}_dataset_rvv_fast.csv"
COEFS_PATH = WORKDIR / f"{PARTICLE1_NAME}_{PARTICLE2_NAME}_regression_coefs.csv"


# ============== ТВОИ ФУНКЦИИ (БЕЗ ИЗМЕНЕНИЙ) ==============
def mape_loss(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return float(np.mean(np.abs((y_true - y_pred) / y_true)) * 100.0)


def k_model_poly(T, *coeffs):
    x = np.power(np.asarray(T, dtype=float), -1.0 / 3.0)
    out = np.zeros_like(x, dtype=float)
    for k, c in enumerate(coeffs):
        out += float(c) * np.power(x, k)
    return out


def fit_poly(df: pd.DataFrame):
    """Как в `regression.ipynb`: для каждой числовой колонки после `T` подбирает степень 3…6 параметров LM + MAPE."""
    params_arr = []
    pred_df = pd.DataFrame(index=df.index)
    for col_i in range(1, df.shape[1]):
        coeffs = None
        pred = np.zeros(len(df))
        best_mape = np.inf

        series = pd.to_numeric(df.iloc[:, col_i], errors="coerce")
        mask = np.isfinite(series.values) & (series.values > 0)
        T_col = pd.to_numeric(df.iloc[:, 0], errors="coerce")
        T_valid = T_col.values[mask]
        y_valid = series.values[mask]

        if len(T_valid) < 4:
            raise ValueError("Недостаточно положительных точек для полинома.")

        y_log = np.log(y_valid)
        norm_factor = float(np.max(np.abs(y_log)))
        if norm_factor == 0:
            norm_factor = 1.0
        y_norm = y_log / norm_factor

        for coeff_num in range(4, 8):
            p0 = np.ones(coeff_num)
            params, _ = curve_fit(
                k_model_poly,
                T_valid,
                y_norm,
                p0=p0,
                method="lm",
                maxfev=40000,
            )

            params1 = params * norm_factor
            y_pred_log_full = k_model_poly(T_col.values.astype(float), *params1)
            k_pred_full = np.exp(y_pred_log_full)

            valid_pred = np.isfinite(k_pred_full) & mask
            if not np.any(valid_pred):
                continue
            mape_poly = mape_loss(series.values[valid_pred], k_pred_full[valid_pred])

            if mape_poly < best_mape:
                best_mape = mape_poly
                coeffs = params1.astype(float)
                pred = k_pred_full

        if coeffs is None:
            raise RuntimeError("curve_fit не сработала ни для одного порядка.")

        params_arr.append(coeffs)
        pred_df[col_i - 1] = pred

    return params_arr, pred_df


# ============== ЗАГРУЗКА И ПРЕОБРАЗОВАНИЕ ==============
print("Загрузка данных...")
df_long = pd.read_csv(DATASET_RVV_FAST_PATH)

# Приводим типы
df_long['T'] = df_long['T'].astype(float)
df_long['i1'] = df_long['i1'].astype(int)
df_long['f1'] = df_long['f1'].astype(int)
df_long['i2'] = df_long['i2'].astype(int)
df_long['f2'] = df_long['f2'].astype(int)
df_long['k_VV'] = df_long['k_VV'].astype(float)

print(f"Загружено строк: {len(df_long)}")

# Преобразуем в wide формат (как в твоем оригинале)
print("Преобразование в wide формат...")
df_long['transition'] = df_long.apply(
    lambda row: f"{row['i1']}_{row['f1']}_{row['i2']}_{row['f2']}", axis=1
)

# Pivot
df_wide = df_long.pivot(index='T', columns='transition', values='k_VV')

print(f"Wide формат: {df_wide.shape[0]} температур × {df_wide.shape[1]} переходов")

# Сбрасываем индекс чтобы T стала колонкой
df_wide = df_wide.reset_index()

# Переставляем колонки: первая T, остальные - переходы
cols = ['T'] + [c for c in df_wide.columns if c != 'T']
df_wide = df_wide[cols]

print(f"Готово для fit_poly: {df_wide.shape}")

# ============== ЗАПУСК РЕГРЕССИИ ==============
print("\nЗапуск регрессии (это может занять время)...")
print("=" * 60)

try:
    params_arr, pred_df = fit_poly(df_wide)
    print(f"Регрессия завершена! Обработано переходов: {len(params_arr)}")

    # ============== СОХРАНЕНИЕ КОЭФФИЦИЕНТОВ (ИСПРАВЛЕННЫЙ) ==============
    # Получаем названия переходов из колонок (пропуская T)
    transitions = df_wide.columns[1:]

    coeffs_rows = []
    for col_name, coeffs in zip(transitions, params_arr):
        # Парсим переход - убираем возможные .0
        parts = col_name.split('_')
        # Очищаем от .0 и преобразуем в int
        i1 = int(float(parts[0]))
        f1 = int(float(parts[1]))
        i2 = int(float(parts[2]))
        f2 = int(float(parts[3]))

        row = {
            'i1': i1, 'f1': f1,
            'i2': i2, 'f2': f2
        }

        # Добавляем коэффициенты
        for j, c in enumerate(coeffs):
            row[f'a_{j}'] = c

        coeffs_rows.append(row)

    df_coeffs = pd.DataFrame(coeffs_rows)
    df_coeffs.to_csv(COEFS_PATH, index=False)

    print(f"\n✅ Сохранено: {COEFS_PATH.resolve()}")
    print(f"Всего переходов: {len(df_coeffs)}")
    print("\nПервые 5 результатов:")
    print(df_coeffs.head())

except Exception as e:
    print(f"Ошибка: {e}")
    import traceback

    traceback.print_exc()