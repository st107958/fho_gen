from scipy.integrate import nquad, trapezoid
from scipy.integrate import simpson
import numpy as np

from numba import cuda

#np.set_printoptions(threshold=np.inf)
import cupy as cp

from constants import *
from p_vv_mm_ij import p_vv, g, gamma
from particles_data import *

# threads_per_block = 256
# blocks_per_grid = 65535
# @cuda.jit
# def trapz_cuda_kernel(F, eps1, eps2, y, v1, phi1, v2, phi2, result):
#     idx = cuda.grid(1)
#     n0, n1, n2, n3, n4, n5, n6 = F.shape
#
#     if idx < n0 * n1 * n2 * n3 * n4 * n5 * n6:
#         # Преобразование линейного индекса в 7D координаты
#         i0 = idx // (n1 * n2 * n3 * n4 * n5 * n6)
#         remainder = idx % (n1 * n2 * n3 * n4 * n5 * n6)
#         i1 = remainder // (n2 * n3 * n4 * n5 * n6)
#         remainder %= (n2 * n3 * n4 * n5 * n6)
#         i2 = remainder // (n3 * n4 * n5 * n6)
#         remainder %= (n3 * n4 * n5 * n6)
#         i3 = remainder // (n4 * n5 * n6)
#         remainder %= (n4 * n5 * n6)
#         i4 = remainder // (n5 * n6)
#         remainder %= (n5 * n6)
#         i5 = remainder // n6
#         i6 = remainder % n6
#
#         # Вычисление весов
#         weight = 1.0
#
#         # Вес для eps1 (ось 6)
#         if i6 == 0 or i6 == n6 - 1:
#             weight *= (eps1[1] - eps1[0]) / 2
#         else:
#             weight *= (eps1[i6 + 1] - eps1[i6 - 1]) / 2
#
#         # Вес для eps2 (ось 5)
#         if i5 == 0 or i5 == n5 - 1:
#             weight *= (eps2[1] - eps2[0]) / 2
#         else:
#             weight *= (eps2[i5 + 1] - eps2[i5 - 1]) / 2
#
#         # Вес для y (ось 4)
#         if i4 == 0 or i4 == n4 - 1:
#             weight *= (y[1] - y[0]) / 2
#         else:
#             weight *= (y[i4 + 1] - y[i4 - 1]) / 2
#
#         # Вес для v1 (ось 3)
#         if i3 == 0 or i3 == n3 - 1:
#             weight *= (v1[1] - v1[0]) / 2
#         else:
#             weight *= (v1[i3 + 1] - v1[i3 - 1]) / 2
#
#         # Вес для phi1 (ось 2)
#         if i2 == 0 or i2 == n2 - 1:
#             weight *= (phi1[1] - phi1[0]) / 2
#         else:
#             weight *= (phi1[i2 + 1] - phi1[i2 - 1]) / 2
#
#         # Вес для v2 (ось 1)
#         if i1 == 0 or i1 == n1 - 1:
#             weight *= (v2[1] - v2[0]) / 2
#         else:
#             weight *= (v2[i1 + 1] - v2[i1 - 1]) / 2
#
#         # Вес для phi2 (ось 0)
#         if i0 == 0 or i0 == n0 - 1:
#             weight *= (phi2[1] - phi2[0]) / 2
#         else:
#             weight *= (phi2[i0 + 1] - phi2[i0 - 1]) / 2
#
#         # Атомарное суммирование
#         cuda.atomic.add(result, 0, F[i0, i1, i2, i3, i4, i5, i6] * weight)
#
#
# def compute_7d_integral(F, eps1, eps2, y, v1, phi1, v2, phi2):
#     # Перенос данных на GPU
#     F_gpu = cuda.to_device(F)
#     eps1_gpu = cuda.to_device(eps1)
#     eps2_gpu = cuda.to_device(eps2)
#     y_gpu = cuda.to_device(y)
#     v1_gpu = cuda.to_device(v1)
#     phi1_gpu = cuda.to_device(phi1)
#     v2_gpu = cuda.to_device(v2)
#     phi2_gpu = cuda.to_device(phi2)
#     result_gpu = cuda.device_array(1, dtype=np.float64)
#
#     # Инициализация результата
#     result_gpu[0] = 0.0
#
#     # Запуск ядра
#     trapz_cuda_kernel[blocks_per_grid, threads_per_block](
#         F_gpu, eps1_gpu, eps2_gpu, y_gpu, v1_gpu, phi1_gpu, v2_gpu, phi2_gpu, result_gpu
#     )
#
#     # Возвращаем результат
#     return result_gpu.copy_to_host()[0]

def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapez'):
    e_in_J = h * c * 100
    E = E * e_in_J   # 1/cm --> J

    # print('E', E)

    if m1 == m2 and i1 == f2 and i2 == f1:
        raise ValueError("ksi = 0, resonance process")

    if method == 'trapez':
        maxdiv = 5  # макс. кол-во делений по координате: 18 (больше - много памяти)

        # пределы интегрирования
        eps1 = np.linspace(0, 1, maxdiv)
        eps2 = np.linspace(0, 1, maxdiv)
        y = np.linspace(0, 1, maxdiv)
        v1 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        phi1 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        v2 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        phi2 = np.linspace(-np.pi/2, np.pi/2, maxdiv)

        EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2, indexing='ij')


        # print('eps1 grid shape', EPS1.shape)

        # mask = (EPS1 + EPS2) <= (1-np.power(Y, 2))/(2*(1-np.power(Y, 2)/2))
        mask = (EPS1 + EPS2) <= 1/2
        #############################################

        EPS1_filtered = EPS1[mask]
        EPS2_filtered = EPS2[mask]
        Y_filtered = Y[mask]
        V1_filtered = V1[mask]
        PHI1_filtered = PHI1[mask]
        V2_filtered = V2[mask]
        PHI2_filtered = PHI2[mask]


        # print('eps1_filtered grid shape', EPS1_restored.shape)

        F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1_filtered, EPS2_filtered, Y_filtered,
                 V1_filtered, PHI1_filtered, V2_filtered, PHI2_filtered)

        # F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1_restored, EPS2_restored, Y_restored,
        #          V1_restored, PHI1_restored, V2_restored, PHI2_restored)

        # print('F', F.shape)

        F_ = np.zeros_like(EPS1)  # Исходная форма (maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv)
        F_[mask] = F  # допустимые точки
        F_ = np.nan_to_num(F_)

        # print('F_reshaped', F_.shape)

        #1

        # result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
        #     F_, eps1, axis=6), eps2, axis=5), y, axis=4), v1, axis=3),
        #     phi1, axis=2), v2, axis=1), phi2, axis=0)

        #2



        #3

        result = np.trapz(
            np.trapz(
                np.trapz(
                    np.trapz(
                        np.trapz(
                            np.trapz(
                                np.trapz(F_, eps1, axis=6),
                                eps2, axis=5),
                            y, axis=4),
                        v1, axis=3),
                    phi1, axis=2),
                v2, axis=1),
            phi2, axis=0)

        # print(result)

        #############################################################

        # F = np.nan_to_num(p_vv(m1, m2, i1, f1, i2, f2, E, EPS1, EPS2, Y, V1, PHI1, V2, PHI2) * mask, nan=0.0)
        #
        # result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
        #     F, eps1, axis=6), eps2, axis=5), y, axis=4), v1, axis=3),
        #     phi1, axis=2), v2, axis=1), phi2, axis=0)
        #
        # print(result)

        result = result / (np.pi ** 4)

    return result

# print(p_vv_int(CO, CO, 2, 0, 1, 3, 10000, 'trapez'))
#
#
#
# print(p_vv_int(N2, N2, 1, 0, 2, 3, 10000, 'trapez'))
#
# print(p_vv_int(N2, N2, 4, 2, 0, 3, 10000, 'trapez'))
#
# print(p_vv_int(N2, N2, 5, 2, 0, 3, 10000, 'trapez'))
# print(p_vv_int(N2, N2, 5, 1, 0, 4, 10000, 'trapez'))
# print(p_vv_int(N2, N2, 6, 1, 0, 5, 10000, 'trapez'))



