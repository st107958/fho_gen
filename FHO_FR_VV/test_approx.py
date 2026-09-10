import time
from pathlib import Path
import numpy as np
import pandas as pd

# ============== КОНФИГ ==============
PARTICLE1_NAME = "N2"
PARTICLE2_NAME = "N2"
WORKDIR = Path(".")
DATASET_RVV_FAST_PATH = WORKDIR / f"{PARTICLE1_NAME}_{PARTICLE2_NAME}_dataset_rvv_fast.csv"

# ============== ЗАМЕР ТОЛЬКО ЧТЕНИЯ И ПРЕОБРАЗОВАНИЯ ==============
print("=" * 60)
print("ЗАМЕР: чтение датасета и преобразование в wide формат")
print("=" * 60)

# 1. Чтение CSV
t0_read = time.time()
df_long = pd.read_csv(DATASET_RVV_FAST_PATH)
t1_read = time.time()
print(f"1. Чтение CSV: {t1_read - t0_read:.3f} сек")
print(f"   Размер: {len(df_long)} строк × {len(df_long.columns)} колонок")
print(f"   Память (приблизительно): {df_long.memory_usage(deep=True).sum() / 1024**2:.1f} MB")

# 2. Приведение типов
t0_convert = time.time()
df_long['T'] = df_long['T'].astype(float)
df_long['i1'] = df_long['i1'].astype(int)
df_long['f1'] = df_long['f1'].astype(int)
df_long['i2'] = df_long['i2'].astype(int)
df_long['f2'] = df_long['f2'].astype(int)
df_long['k_VV'] = df_long['k_VV'].astype(float)
t1_convert = time.time()
print(f"2. Приведение типов: {t1_convert - t0_convert:.3f} сек")

# 3. Создание колонки transition (строковая операция)
t0_transition = time.time()
df_long['transition'] = df_long.apply(
    lambda row: f"{row['i1']}_{row['f1']}_{row['i2']}_{row['f2']}", axis=1
)
t1_transition = time.time()
print(f"3. Создание transition (apply): {t1_transition - t0_transition:.3f} сек")

# 4. Pivot (wide format) - самая тяжелая операция
t0_pivot = time.time()
df_wide = df_long.pivot(index='T', columns='transition', values='k_VV')
t1_pivot = time.time()
print(f"4. Pivot (wide format): {t1_pivot - t0_pivot:.3f} сек")

# 5. Сброс индекса
t0_reset = time.time()
df_wide = df_wide.reset_index()
t1_reset = time.time()
print(f"5. Сброс индекса: {t1_reset - t0_reset:.3f} сек")

# 6. Перестановка колонок
t0_reorder = time.time()
cols = ['T'] + [c for c in df_wide.columns if c != 'T']
df_wide = df_wide[cols]
t1_reorder = time.time()
print(f"6. Перестановка колонок: {t1_reorder - t0_reorder:.3f} сек")

# ============== ИТОГО ==============
total_time = t1_reorder - t0_read
print(f"\n" + "=" * 60)
print(f"ИТОГО ВРЕМЯ: {total_time:.3f} сек")
print(f"  из них чтение CSV: {t1_read - t0_read:.3f} сек ({100*(t1_read-t0_read)/total_time:.1f}%)")
print(f"  из них pivot: {t1_pivot - t0_pivot:.3f} сек ({100*(t1_pivot-t0_pivot)/total_time:.1f}%)")
print("=" * 60)

# Дополнительная информация
print(f"\n📊 Итоговый wide формат: {df_wide.shape[0]} температур × {df_wide.shape[1]} переходов")
print(f"   (температур обычно {len(df_long['T'].unique())}, переходов ~{df_wide.shape[1]-1})")

# Проверка что все корректно
print("\n✅ Первые 3 строки (первые 5 колонок):")
print(df_wide.iloc[:3, :5].to_string())