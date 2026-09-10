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


    # ============== ПОСЛЕ РЕГРЕССИИ - ПОИСК МАКСИМАЛЬНОГО MAPE ==============
    print("\n" + "=" * 60)
    print("ПОИСК МАКСИМАЛЬНОЙ ОШИБКИ ПО ВСЕМ ПЕРЕХОДАМ")
    print("=" * 60)

    from scipy.interpolate import interp1d

    max_mape = 0
    max_mape_transition = None
    all_mapes = []

    for idx, col_name in enumerate(df_wide.columns[1:]):
        # Получаем коэффициенты для этого перехода
        coeffs = params_arr[idx]

        # Предсказываем при тех же температурах
        T_values = df_wide['T'].values
        k_true = df_wide[col_name].values

        # Предсказание через exp(полином)
        pred = np.exp(k_model_poly(T_values, *coeffs))

        # Убираем нулевые/отрицательные значения
        mask = (k_true > 0) & np.isfinite(k_true) & np.isfinite(pred)

        if mask.any():
            mape = np.mean(np.abs((k_true[mask] - pred[mask]) / k_true[mask])) * 100

            if mape > max_mape:
                max_mape = mape
                max_mape_transition = col_name

            all_mapes.append(mape)

    print(f"📊 МАКСИМАЛЬНЫЙ MAPE: {max_mape:.6f}%")
    print(f"🔀 Переход: {max_mape_transition}")
    print(f"📈 Средний MAPE: {np.mean(all_mapes):.6f}%")
    print(f"📉 Медианный MAPE: {np.median(all_mapes):.6f}%")
    print(f"🎯 95-й перцентиль: {np.percentile(all_mapes, 95):.6f}%")
    print("=" * 60)

    # ============== ПОСТРОЕНИЕ ГРАФИКОВ ==============
    print("\nПостроение графиков...")
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d


    # Функция для предсказания по коэффициентам
    def k_model_for_pred(T, *coeffs):
        x = np.power(np.asarray(T, dtype=float), -1.0 / 3.0)
        out = np.zeros_like(x, dtype=float)
        for k, c in enumerate(coeffs):
            out += float(c) * np.power(x, k)
        return np.exp(out)


    # Функция для форматирования подписи перехода в LaTeX
    def format_transition_label(i1, f1, i2, f2):
        """Форматирует переход в вид: $k_{41,0 \\rightarrow 40,1}^{VV}$"""
        return f"$k_{{{i1},{f1} \\rightarrow {i2},{f2}}}^{{VV}}$"


    # Создаём словарь с предсказаниями для каждого перехода
    predictions = {}
    # Создаём более плотную сетку температур для гладкой аппроксимации
    T_min = df_wide['T'].min()
    T_max = df_wide['T'].max()
    T_dense = np.linspace(T_min, T_max, 200)  # 200 точек вместо 10

    for _, row in df_coeffs.iterrows():
        i1, f1, i2, f2 = int(row['i1']), int(row['f1']), int(row['i2']), int(row['f2'])
        # Формируем ключ перехода
        trans_key = f"{i1}_{f1}_{i2}_{f2}"
        # Формируем LaTeX подпись
        trans_label = format_transition_label(i1, f1, i2, f2)
        # Собираем коэффициенты (все колонки кроме i1,f1,i2,f2)
        coeff_cols = [col for col in df_coeffs.columns if col.startswith('a_')]
        coeffs = [row[col] for col in coeff_cols]
        # Предсказываем на плотной сетке
        predictions[trans_key] = {'label': trans_label, 'pred': k_model_for_pred(T_dense, *coeffs)}

    # Выбираем 3 перехода для отображения
    transition_keys = list(predictions.keys())
    selected_transitions = [
        transition_keys[0],  # первый
        transition_keys[len(transition_keys) // 2],  # средний
        transition_keys[-1]  # последний
    ]

    print(f"Выбраны переходы для отображения:")
    for t in selected_transitions:
        print(f"  - {predictions[t]['label']}")

    # Строим ОДИН график со всеми тремя переходами
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Цвета и стили для разных переходов
    colors = ['blue', 'green', 'red']
    markers = ['o', 's', '^']

    for idx, trans_key in enumerate(selected_transitions):
        color = colors[idx % len(colors)]
        marker = markers[idx % len(markers)]

        # Получаем данные перехода
        trans_data = predictions[trans_key]
        trans_label = trans_data['label']

        # Исходные данные из длинного формата
        parts = trans_key.split('_')
        mask = (df_long['i1'] == int(parts[0])) & \
               (df_long['f1'] == int(parts[1])) & \
               (df_long['i2'] == int(parts[2])) & \
               (df_long['f2'] == int(parts[3]))

        df_trans = df_long[mask].copy()
        df_trans = df_trans.sort_values('T')

        T_orig = df_trans['T'].values
        k_orig = df_trans['k_VV'].values

        # Предсказанные значения на плотной сетке
        k_pred = trans_data['pred']

        # Строим исходные точки
        ax.scatter(T_orig, k_orig, color=color, s=50, alpha=0.7,
                   marker=marker, label=f'{trans_label} (Точные значения FHO-FR)', zorder=3)

        # Строим аппроксимацию пунктирной линией на плотной сетке
        ax.plot(T_dense, k_pred, color=color, linestyle='--', linewidth=2.5,
                label=f'{trans_label} (Аппроксимация)', zorder=2)

        # # Строим исходные точки
        # ax.scatter(T_orig, k_orig, color=color, s=50, alpha=0.7,
        #            marker=marker, label=f'{trans_label} (FHO-FR exact values)', zorder=3)
        #
        # # Строим аппроксимацию пунктирной линией на плотной сетке
        # ax.plot(T_dense, k_pred, color=color, linestyle='--', linewidth=2.5,
        #         label=f'{trans_label} (approximation)', zorder=2)

        # Вычисляем и выводим MAPE в консоль
        interp_func = interp1d(T_dense, k_pred, kind='linear',
                               bounds_error=False, fill_value=np.nan)
        k_pred_at_orig = interp_func(T_orig)

        valid_mask = np.isfinite(k_orig) & (k_orig > 0) & np.isfinite(k_pred_at_orig)
        if valid_mask.any():
            mape = np.mean(np.abs((k_orig[valid_mask] - k_pred_at_orig[valid_mask]) /
                                  k_orig[valid_mask])) * 100
            print(f"  {trans_label}: MAPE = {mape:.10f}%")

    # Настройки графика
    ax.set_xlabel('Температура T, K', fontsize=18)
    ax.set_ylabel(r'Коэффициент скорости энергообмена $k^{VV}$, см³/с', fontsize=18)
    # ax.set_xlabel('Temperature T, K', fontsize=18)
    # ax.set_ylabel(r'Rate Coefficients $k^{VV}$, cm³/s', fontsize=18)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_yscale('log')  # Логарифмическая шкала
    ax.legend(loc='best', fontsize=14, framealpha=0.9)
    ax.tick_params(axis='both', labelsize=14)  # размер делений

    plt.tight_layout()
    plot_path = WORKDIR / f'{PARTICLE1_NAME}_{PARTICLE2_NAME}_regression_fits.png'
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.show()

    print(f"\n✅ График сохранён: {plot_path.resolve()}")

    # ============== ДОБАВЛЯЕМ НА ГРАФИК k_vv_fhofer ==============
    print("\n" + "=" * 60)
    print("ДОБАВЛЯЕМ k_vv_fhofer НА ГРАФИК")
    print("=" * 60)

    from k_VV import k_vv_fhofer
    from particles import N2

    # Температуры для расчёта k_vv_fhofer
    T_fhofer = 300  # или какая у тебя

    # Выбираем переходы для сравнения (как в ноутбуке)
    # Переходы: (i1=1, f1=0, i2=v-1, f2=v) для v от 1 до 9
    v_levels = np.arange(1, 10)
    k_fhofer_values = []

    for v in v_levels:
        try:
            k = k_vv_fhofer(N2, N2, i1=1, f1=0, i2=v - 1, f2=v, T=T_fhofer)
            k_fhofer_values.append(k)
            print(f"v={v}: k_vv_fhofer = {k:.2e}")
        except Exception as e:
            print(f"Ошибка для v={v}: {e}")
            k_fhofer_values.append(np.nan)

    k_fhofer_values = np.array(k_fhofer_values)

    # ============== ДОБАВЛЯЕМ НА СУЩЕСТВУЮЩИЙ ГРАФИК ==============
    # Находим или создаём figure
    import matplotlib.pyplot as plt

    # Если есть открытая фигура - используем её, иначе создаём новую
    fig = plt.gcf()
    ax = plt.gca()

    # Строим k_vv_fhofer поверх существующего графика
    ax.plot(v_levels, k_fhofer_values, 'o-', color='red', linewidth=2,
            markersize=8, label=f'k_vv_fhofer (T={T_fhofer} K)', zorder=10)

    # Обновляем легенду
    ax.legend(loc='best', fontsize=14, framealpha=0.9)

    # Обновляем график
    plt.tight_layout()

    # Сохраняем обновлённый график
    plot_path_updated = WORKDIR / f'{PARTICLE1_NAME}_{PARTICLE2_NAME}_regression_fits_with_fhofer.png'
    plt.savefig(plot_path_updated, dpi=150, bbox_inches='tight')
    print(f"\n✅ Обновлённый график сохранён: {plot_path_updated.resolve()}")

    plt.show()

except Exception as e:
    print(f"Ошибка: {e}")
    import traceback

    traceback.print_exc()


