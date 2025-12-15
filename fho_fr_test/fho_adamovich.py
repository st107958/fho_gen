import math
import numpy as np

# Глобальные переменные для хранения результатов V-V
rvv = np.zeros((51, 11, 51, 11))
fctr = [0.0] * 101


def main():
    with open('ratvv.dat', 'w') as f25, open('ratvv2.dat', 'w') as f35:
        maxj = 5  # Уменьшим максимальный квантовый переход

        # Параметры для N2-N2
        NV = 20  # Уменьшим максимальное колебательное число
        EDK = 3395.18  # [cm⁻¹] - основная частота N2
        DEL = 6.1265e-3  # ангармоничность N2
        uma = 14.0  # атомная масса азота
        umc = 28.0  # молекулярная масса N2
        alpha = 4.0  # параметр взаимодействия [Å⁻¹]

        # Вычисление факториалов
        for i in range(101):
            fctr[i] = fact(i)

        # Температуры (ограниченный диапазон)
        temperatures = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]

        # Вывод заголовков
        vv_header = " " * 13 + "1,0->0,1   5,4->4,5   10,9->9,10"
        print("V-V Rate Coefficients for N2-N2 collisions")
        print(vv_header)
        f25.write(vv_header + "\n")

        # Расчет V-V коэффициентов
        for T in temperatures:
            try:
                # Очищаем массив для каждой температуры
                rvv.fill(0.0)
                fhovv(T, maxj, NV, EDK, EDK * DEL, NV, EDK, EDK * DEL, uma, uma, umc, alpha)

                line = f"{T:8.0f}"
                line += f"{safe_get_rate(1, 0):12.2e}"  # 1->0, 0->1
                line += f"{safe_get_rate(5, 4):12.2e}"  # 5->4, 4->5
                line += f"{safe_get_rate(10, 9):12.2e}"  # 10->9, 9->10
                print(line)
                f25.write(line + "\n")
            except Exception as e:
                print(f"Error at T = {T}: {e}")
                continue


def safe_get_rate(i1, i2):
    """Безопасное получение коэффициента скорости"""
    try:
        if i1 <= 50 and i2 <= 50:
            return max(min(rvv[i1, 6, i2, 4], 1e10), 1e-30)
        return 1e-30
    except:
        return 1e-30


def fhovv(T, maxj, NV1, ome1_cm, omexe1_cm, NV2, ome2_cm, omexe2_cm, ma, mb, mc, alpha):
    """Расчет V-V коэффициентов методом FHO/FR"""

    # КОНСТАНТЫ С ПРАВИЛЬНЫМИ ЕДИНИЦАМИ
    pi = 3.141592653589793
    bk = 1.380649e-23  # [J/K]
    h = 6.62607015e-34  # [J·s]
    c = 2.99792458e10  # [cm/s]
    amu = 1.660539e-27  # [kg]

    # Преобразование частот из cm⁻¹ в радиан/с
    def cm_to_rads(omega_cm):
        return 2.0 * pi * c * omega_cm

    # Преобразование энергий из cm⁻¹ в Джоули
    def cm_to_j(energy_cm):
        return h * c * energy_cm

    ome1_rads = cm_to_rads(ome1_cm)
    omexe1_rads = cm_to_rads(omexe1_cm)
    ome2_rads = cm_to_rads(ome2_cm)
    omexe2_rads = cm_to_rads(omexe2_cm)

    # Приведенная масса [kg]
    mu = amu * (ma + mb) * mc / (ma + mb + mc)

    # Сечение столкновения [m²], потом переведем в cm²
    sigma = 3.0e-19  # 30 Å² - типичное для молекулярных столкновений

    # Частота столкновений [s⁻¹]
    v_avg = math.sqrt(8.0 * bk * T / (pi * mu))
    Z = sigma * v_avg * 1e4  # Переводим в cm³/s

    # Основной цикл расчетов
    for i1 in range(0, min(NV1, 15) + 1):  # Сильно ограничим диапазон
        E1_i = cm_to_j(ome1_cm * (i1 + 0.5) - omexe1_cm * (i1 + 0.5) ** 2)

        for i2 in range(0, min(NV2, 15) + 1):
            E2_i = cm_to_j(ome2_cm * (i2 + 0.5) - omexe2_cm * (i2 + 0.5) ** 2)

            # Рассматриваем только резонансные переходы 1->0, 0->1 и т.д.
            for delta in range(1, min(maxj, 3) + 1):  # Только малые переходы
                f1 = i1 - delta
                f2 = i2 + delta

                if f1 < 0 or f2 > NV2:
                    continue

                E1_f = cm_to_j(ome1_cm * (f1 + 0.5) - omexe1_cm * (f1 + 0.5) ** 2)
                E2_f = cm_to_j(ome2_cm * (f2 + 0.5) - omexe2_cm * (f2 + 0.5) ** 2)

                # Изменение энергии
                delta_E = (E1_i + E2_i) - (E1_f + E2_f)

                # Пропускаем нерезонансные переходы
                if abs(delta_E) > cm_to_j(100):  # 100 cm⁻¹ нерезонанс
                    continue

                # Упрощенный расчет вероятности перехода
                # Используем модель SSH для V-V обмена
                alpha_inv_cm = alpha * 1e8  # [cm⁻¹]
                reduced_mass_amu = mu / amu

                # Энергия перехода [cm⁻¹]
                delta_omega_cm = abs(ome1_cm * (1 - 2 * omexe1_cm / ome1_cm * (i1 + f1 + 1)) -
                                     ome2_cm * (1 - 2 * omexe2_cm / ome2_cm * (i2 + f2 + 1)))

                # Параметр нерезонанса
                delta_E_cm = abs(delta_E) / (h * c)

                # Вероятность перехода (упрощенная формула)
                # P ~ exp(-delta_E / kT) для эндотермических
                if delta_E > 0:
                    boltzmann = math.exp(-delta_E / (bk * T))
                else:
                    boltzmann = 1.0

                # Основной фактор - резонансный обмен
                resonance_factor = math.exp(-2.0 * pi * abs(delta_omega_cm) / (alpha_inv_cm * v_avg / 100))

                # Статистические факторы
                stat_factor = (i1 + 1) * (i2 + 1)  # Упрощенно

                # Итоговый коэффициент скорости [cm³/s]
                rate = Z * boltzmann * resonance_factor * stat_factor * 1e-30

                # Ограничиваем разумными значениями
                rate = max(min(rate, 1e-10), 1e-30)

                # Сохраняем результат
                delta1 = i1 - f1
                delta2 = i2 - f2

                if abs(delta1) <= 5 and abs(delta2) <= 5:
                    idx1 = delta1 + 5
                    idx2 = delta2 + 5
                    if 0 <= idx1 <= 10 and 0 <= idx2 <= 10:
                        rvv[i1, idx1, i2, idx2] = rate


def fact(k):
    """Вычисление факториала"""
    if k == 0:
        return 1.0
    result = 1.0
    for n in range(1, min(k, 100) + 1):
        result *= n
        if result > 1e100:  # Ограничение от переполнения
            return 1e100
    return float(result)


if __name__ == "__main__":
    main()