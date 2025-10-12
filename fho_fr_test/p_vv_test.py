from scipy.integrate import nquad, trapezoid
from scipy.integrate import simpson
import numpy as np
import time

from numba import cuda

#np.set_printoptions(threshold=np.inf)
import cupy as cp

from fho_fr_test.constants import *
from fho_fr_test.p_vv_ij_test import p_vv
# from fho_fr_test.p_vv_mm_ij import p_vv
from fho_fr_test.particles_data import *


print(cp.cuda.runtime.getDeviceCount())  # Число доступных GPU
current_gpu = cp.cuda.Device().id
print("Текущий GPU:", current_gpu)
print("CuPy использует CUDA?", cp.cuda.is_available())


def p_vv_int(m1, m2, i1, f1, i2, f2, E, method='trapez'):
    e_in_J = h * c * 100
    E = E * e_in_J   # 1/cm --> J
    # print(f"i1: {i1}, f1: {f1}")


    # print('E', E)

    if m1 == m2 and i1 == f2 and i2 == f1:
        raise ValueError("ksi = 0, resonance process")

    if method == 'trapez':
        # maxdiv = 5  # макс. кол-во делений по координате: 18 (больше - много памяти)
        #
        # # пределы интегрирования
        # eps1 = np.linspace(0, 1, maxdiv)
        # eps2 = np.linspace(0, 1, maxdiv)
        # y = np.linspace(0, 1, maxdiv)
        # v1 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        # phi1 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        # v2 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        # phi2 = np.linspace(-np.pi/2, np.pi/2, maxdiv)
        #
        # EPS1, EPS2, Y, V1, PHI1, V2, PHI2 = np.meshgrid(eps1, eps2, y, v1, phi1, v2, phi2, indexing='ij')
        #
        #
        # # print('eps1 grid shape', EPS1.shape)
        #
        # # mask = (EPS1 + EPS2) <= (1-np.power(Y, 2))/(2*(1-np.power(Y, 2)/2))
        # mask = (EPS1 + EPS2) <= 1/2
        # #############################################
        #
        # EPS1_filtered = EPS1[mask]
        # EPS2_filtered = EPS2[mask]
        # Y_filtered = Y[mask]
        # V1_filtered = V1[mask]
        # PHI1_filtered = PHI1[mask]
        # V2_filtered = V2[mask]
        # PHI2_filtered = PHI2[mask]
        #
        #
        # # print('eps1_filtered grid shape', EPS1_restored.shape)
        #
        # F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1_filtered, EPS2_filtered, Y_filtered,
        #          V1_filtered, PHI1_filtered, V2_filtered, PHI2_filtered)
        #
        # # F = p_vv(m1, m2, i1, f1, i2, f2, E, EPS1_restored, EPS2_restored, Y_restored,
        # #          V1_restored, PHI1_restored, V2_restored, PHI2_restored)
        #
        # # print('F', F.shape)
        #
        # F_ = np.zeros_like(EPS1)  # Исходная форма (maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv, maxdiv)
        # F_[mask] = F  # допустимые точки
        # F_ = np.nan_to_num(F_)

        # print('F_reshaped', F_.shape)

        #1

        # result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
        #     F_, eps1, axis=6), eps2, axis=5), y, axis=4), v1, axis=3),
        #     phi1, axis=2), v2, axis=1), phi2, axis=0)

        #2
        # print('ВЫчисление')
        # start_time = time.time_ns()
        #
        # count = 0
        # threads_per_block = 256
        # blocks_per_grid = 65535
        # @cuda.jit(device=True)
        # def integrand(eps1, eps2, y, v1, phi1, v2, phi2):
        #     nonlocal count
        #     count += 1
        #     if count % 1000 == 0:
        #         print(f"Progress: {count} evaluations")
        #
        #     return p_vv(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2)
        #
        #
        # ranges = [
        #     (0, 0.5),  # eps1 (верхний предел 0.5 из-за условия eps1 + eps2 <= 0.5)
        #     (0, 0.5),  # eps2
        #     (0, 1),  # y
        #     (-np.pi / 2, np.pi / 2),  # v1
        #     (-np.pi / 2, np.pi / 2),  # phi1
        #     (-np.pi / 2, np.pi / 2),  # v2
        #     (-np.pi / 2, np.pi / 2)  # phi2
        # ]
        #
        # # Вычисляем интеграл с адаптацией
        # result, error = nquad(integrand, ranges)
        #
        # end_time = time.time_ns()
        # print((end_time - start_time) * 1e-9)

        #3

        maxdiv = 4

        eps1 = cp.linspace(0, 1, maxdiv)
        eps2 = cp.linspace(0, 1, maxdiv)
        y = cp.linspace(0, 1, maxdiv)
        v1 = cp.linspace(-np.pi/2, np.pi/2, maxdiv)
        phi1 = cp.linspace(-np.pi/2, np.pi/2, maxdiv)
        v2 = cp.linspace(-np.pi/2, np.pi/2, maxdiv)
        phi2 = cp.linspace(-np.pi/2, np.pi/2, maxdiv)

        deps1 = eps1[1] - eps1[0]
        deps2 = eps2[1] - eps2[0]
        dy = y[1] - y[0]
        dv1 = v1[1] - v1[0]
        dv2 = v2[1] - v2[0]
        dphi1 = phi1[1] - phi1[0]
        dphi2 = phi2[1] - phi2[0]

        start_time = time.time_ns()

        count = 0
        threads_per_block = 256
        blocks_per_grid = 65535

        integrand_kernel = cp.ElementwiseKernel(
            'float64 eps1_val, float64 eps2_val, float64 y_val, float64 v1, float64 phi1, float64 v2, float64 phi2',
            'float64 val',
            '''
            // Здесь должна быть ваша функция p_vv()
            // Пока используем заглушку
            val = eps1_val * eps2_val * y_val * v1 * phi1 * v2 * phi2;
            ''',
            'integrand_kernel'
        )

        def integrand(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2):
            nonlocal count
            count += 1
            if count % 1000 == 0:
                print(f"Progress: {count} evaluations")

            return p_vv(m1, m2, i1, f1, i2, f2, E, eps1, eps2, y, v1, phi1, v2, phi2)


        result = 0.0
        for k1 in range(maxdiv):
            for k2 in range(maxdiv):
                for k3 in range(maxdiv):
                    eps1_val = eps1[k1]
                    eps2_val = eps2[k2]
                    y_val = y[k3]

                    # GPU батч по остальным 4 переменным
                    V1, PHI1, V2, PHI2 = cp.meshgrid(v1, phi1, v2, phi2, indexing='ij')

                    val = integrand(m1, m2, i1, f1, i2, f2, E, eps1_val, eps2_val, y_val, V1, PHI1, V2, PHI2)
                    result += cp.sum(val)
                    result = result * deps1 * deps2 * dy * dv1 * dv2 * dphi1 * dphi2
                    print(result)

        # result = result * deps1 * deps2 * dy * dv1 * dv2 * dphi1 * dphi2

        end_time = time.time_ns()
        print((end_time - start_time) * 1e-9)

        #4

        # result = np.trapz(
        #     np.trapz(
        #         np.trapz(
        #             np.trapz(
        #                 np.trapz(
        #                     np.trapz(
        #                         np.trapz(F_, eps1, axis=6),
        #                         eps2, axis=5),
        #                     y, axis=4),
        #                 v1, axis=3),
        #             phi1, axis=2),
        #         v2, axis=1),
        #     phi2, axis=0)

        # print(result)

        #############################################################

        # F = np.nan_to_num(p_vv(m1, m2, i1, f1, i2, f2, E, EPS1, EPS2, Y, V1, PHI1, V2, PHI2) * mask, nan=0.0)
        #
        # result = trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(trapezoid(
        #     F, eps1, axis=6), eps2, axis=5), y, axis=4), v1, axis=3),
        #     phi1, axis=2), v2, axis=1), phi2, axis=0)
        #
        # print(result)

        result = result / (cp.pi ** 4)

    return result

print(p_vv_int(N2, N2, 1, 0, 20, 21, 10000))
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

# def gamma(..)
# def g(...)
# def p_vv(...)
#     k = l/(2*np.pi)
#     return g(..)*gamma(..)*k
# def integrand()
#     return p_vv(...)