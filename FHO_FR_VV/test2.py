import numpy as np
from scipy.integrate import nquad, trapezoid

# Функция, которую интегрируем
def f(x, y, z):
    return np.sqrt(np.maximum(0, 1 - x - y))  # Замени на свою функцию

# Параметры сетки
Nx, Ny = 100, 100  # Количество точек по x и y
x_vals = np.linspace(0, 1, Nx)
y_vals = np.linspace(0, 1, Ny)
z_vals = np.linspace(0, 1, Ny)

# Создаём сетку
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals)

# Маска для области x + y <= 1
mask = X + Y <= 1

# Применяем функцию к сетке
F = f(X, Y, Z) * mask

# Метод трапеций (по X и Y)
integral = np.trapz(np.trapz(np.trapz(F, x_vals, axis=2), y_vals, axis=1), z_vals)
integral2 = trapezoid(trapezoid(trapezoid(F, x_vals, axis=2), y_vals, axis=1), z_vals)

print(f"Интеграл: {integral:.6f}")
print(integral2)