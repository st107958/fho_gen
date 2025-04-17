import numpy as np
from scipy.integrate import nquad, trapezoid

def f(x, y, z):
    return np.sqrt((1 - x - y)*(1 - z))

def f_1(x, y, z):
    return np.where(x+y < 1, np.sqrt(np.maximum(0, (1 - x - y)*(1 - z))), 0)

# Параметры сетки
Nx, Ny, Nz = 10, 10, 10  # Количество точек по x и y
x_vals = np.linspace(0, 1, Nx)
y_vals = np.linspace(0, 1, Ny)
z_vals = np.linspace(0, 1, Nz)

# Создаём сетку
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals)

# Маска для области x + y <= 1
mask = X + Y <= 1

X_filtered = X[mask]  # 1D-массив допустимых X
Y_filtered = Y[mask]  # 1D-массив допустимых Y
Z_filtered = Z[mask]  # 1D-массив допустимых Z

F = f(X_filtered, Y_filtered, Z_filtered)

F_ = np.zeros_like(X)  # Исходная форма (10,10,10)
F_[mask] = F    # Заполняем только допустимые точки
F_ = np.nan_to_num(F_)  # На всякий случай заменяем оставшиеся NaN на 0

# Применяем функцию к сетке
# F = np.nan_to_num(f(X, Y, Z) * mask, nan=0.0)
# print(F.shape, F)
# F_1 = f_1(X, Y, Z)

# Метод трапеций (по X и Y)
# integral = np.trapz(np.trapz(np.trapz(F, x_vals, axis=2), y_vals, axis=1), z_vals)
# integral2 = trapezoid(trapezoid(trapezoid(F, x_vals, axis=2), y_vals, axis=1), z_vals, axis=0)
# integral3 = trapezoid(trapezoid(trapezoid(F_1, x_vals, axis=2), y_vals, axis=1), z_vals)

# print("F shape:", F.shape)
# print("x_vals shape:", x_vals.shape)
# print("y_vals shape:", y_vals.shape)
# print("z_vals shape:", z_vals.shape)

integral_x = trapezoid(F_, x_vals, axis=0)
integral_xy = trapezoid(integral_x, y_vals, axis=0)
integral_xyz = trapezoid(integral_xy, z_vals, axis=0)

# print(f"Интеграл: {integral:.6f}")
print(integral_xyz)
# print(integral3)

x = [1.0, 2.0, -1.0, 0.0]

x_array1 = np.array(x)

x_array2 = np.array(x)

mask = x_array2 >= 0
x_array2 = x_array2 * mask
# print(x_array2)
#
# print(np.where(x_array2 > 0, np.power(x_array2, 1/2), 0.0))
# print(np.sqrt(0))

