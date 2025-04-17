import numpy as np
from scipy.integrate import nquad, trapezoid

# Функция, которую интегрируем
def f(x, y, z):
    return np.sqrt(np.maximum(0, (1 - x - y)*(1 - z)))  # Замени на свою функцию

def f_1(x, y, z):
    return np.where(x+y < 1, np.sqrt(np.maximum(0, (1 - x - y)*(1 - z))), 0)

# Параметры сетки
Nx, Ny = 100, 100  # Количество точек по x и y
x_vals = np.linspace(0, 10, Nx)
y_vals = np.linspace(0, 1, Ny)
z_vals = np.linspace(0, 1, Ny)

# Создаём сетку
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals)

# Маска для области x + y <= 1
mask = X + Y <= 1

# Применяем функцию к сетке
F = f(X, Y, Z) * mask
F_1 = f_1(X, Y, Z)

# Метод трапеций (по X и Y)
integral = np.trapz(np.trapz(np.trapz(F, x_vals, axis=2), y_vals, axis=1), z_vals)
integral2 = trapezoid(trapezoid(trapezoid(F, x_vals, axis=2), y_vals, axis=1), z_vals)
integral3 = trapezoid(trapezoid(trapezoid(F_1, x_vals, axis=2), y_vals, axis=1), z_vals)

print(f"Интеграл: {integral:.6f}")
print(integral2)
print(integral3)

x = [1.0, 2.0, -1.0, 0.0]

x_array1 = np.array(x)

x_array2 = np.array(x)

mask = x_array2 >= 0
x_array2 = x_array2 * mask
print(x_array2)

print(np.where(x_array2 > 0, np.power(x_array2, 1/2), 0.0))
# print(np.sqrt(0))

