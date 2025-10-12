import time

import cupy as cp

# Проверка доступности GPU
print("Доступно GPU:", cp.cuda.runtime.getDeviceCount())  # Должно быть > 0
print("Текущий GPU:", cp.cuda.Device().id)  # Должно быть 0 или 1
print("Массив создан на GPU?", cp.array([1, 2, 3]).device)  # Должно быть <CUDA Device 0>


# Тест производительности GPU через time.time()
def test_gpu_performance():
    a = cp.random.rand(10000, 10000)
    b = cp.random.rand(10000, 10000)

    start_time = time.time()
    c = a @ b  # Матричное умножение
    cp.cuda.Stream.null.synchronize()  # Синхронизация GPU
    elapsed_time = time.time() - start_time

    print(f"Время умножения матриц 10000x10000: {elapsed_time:.3f} сек")
    print(f"Загрузка GPU можно проверить через 'nvidia-smi'")


test_gpu_performance()