import random
import numpy as np

y_average_mass = []
y_mass = []
teta_mass = []
xn = [[-1, -1], [-1, 1], [1, -1]]
variant = 308

x1_min, x1_max, x2_min, x2_max = -30, 0, -15, 35

y_min, y_max = (20 - variant) * 10, (30 - variant) * 10

for i in range(3):
    y_mass.append([])
    for j in range(7):
        m_a_x = random.randint(y_min, y_max)
        y_mass[i].append(m_a_x)

for i in range(3):
    y_av = (y_mass[i][0] + y_mass[i][1] + y_mass[i][2] + y_mass[i][3] + y_mass[i][4] + y_mass[i][5] + y_mass[i][6]) / 7
    y_average_mass.append(y_av)

for i in range(3):
    teta_sum = 0
    for j in range(7):
        teta_sum += (y_mass[i][j] - y_average_mass[i]) ** 2
    teta_mass.append(teta_sum / 7)

teta = np.sqrt(2 * (2 * 7 - 2) / 7 * (7 - 4))

Fuv1 = teta_mass[0] / teta_mass[1]
Fuv2 = teta_mass[2] / teta_mass[0]
Fuv3 = teta_mass[2] / teta_mass[1]

for i in range(3):
    globals()['Quv' + str(i + 1)] = globals()['Fuv' + str(i + 1)] * ((7 - 2) / 7)

for i in range(3):
    globals()['Ruv' + str(i + 1)] = (globals()['Quv' + str(i + 1)] - 1) / teta

mx1, mx2, my, a1, a2, a3, a11, a22 = 0, 0, 0, 0, 0, 0, 0, 0

for i in range(3):
    mx1 += xn[i][0]
    mx2 += xn[i][1]
    my += y_average_mass[i]
    a1 += xn[i][0] ** 2
    a2 = xn[i][0] * xn[i][1]
    a3 += xn[i][1] ** 2
    a11 += xn[i][0] * y_average_mass[i]
    a22 += xn[i][1] * y_average_mass[i]

mx1, mx2, my, a1, a2, a3, a11, a22 = mx1 / 3, mx2 / 3, my / 3, a1 / 3, a2 / 3, a3 / 3, a11 / 3, a22 / 3

b0 = np.linalg.det([[my, mx1, mx2], [a11, a1, a2], [a22, a2, a3]]) / np.linalg.det(
    [[1, mx1, mx2], [mx1, a1, a2], [mx2, a2, a3]])
b1 = np.linalg.det([[1, my, mx2], [mx1, a11, a2], [mx2, a22, a3]]) / np.linalg.det(
    [[1, mx1, mx2], [mx1, a1, a2], [mx2, a2, a3]])
b2 = np.linalg.det([[1, mx1, my], [mx1, a1, a11], [mx2, a2, a22]]) / np.linalg.det(
    [[1, mx1, mx2], [mx1, a1, a2], [mx2, a2, a3]])

delta_x1 = (x1_max - x1_min) / 2
delta_x2 = (x2_max - x2_min) / 2
x10 = (x1_min + x1_max) / 2
x20 = (x2_min + x2_max) / 2
a_0 = b0 - b1 * (x10 / delta_x1) - b2 * (x20 / delta_x2)
a_1 = b1 / delta_x1
a_2 = b2 / delta_x2

for j in range(3):
    print(xn[j])
for i in range(3):
    print(y_mass[i])
print("Average Y = ", y_average_mass)
print("Teta = ", teta_mass)
print("Відхилення = ", teta)
for i in range(3):
    print("Fuv" + str(i), "=", globals()['Fuv' + str(i + 1)])
for i in range(3):
    print("Quv" + str(i), "=", globals()['Quv' + str(i + 1)])
for i in range(3):
    print("Ruv" + str(i), "=", globals()['Ruv' + str(i + 1)])
print("mx1 = " + str(mx1) + ", ", "mx2 = " + str(mx2) + ", ", "my = " + str(my) + ", ", "a1 = " + str(a1) + ", ",
      "a2 = " + str(a2) + ", ", "a3 = " + str(a3))
print("a11 = " + str(a11) + ",", "a22 = " + str(a22) + ",", "b0 = " + str(b0) + ",", "b1 = " + str(b1) + ",",
      "b2 = " + str(b2))
print("Нормоване рівняння регресії")
print(str(b0) + " " + str(b1) + "*x1" + " " + str(b2) + "*x2")
print("Перевірка")
print(b0 + b1 * (-1) + b2 * (-1), " = ", y_average_mass[0])
print(b0 + b1 * (-1) + b2 * (1), " = ", y_average_mass[1])
print(b0 + b1 * (1) + b2 * (-1), " = ", y_average_mass[2])
print("delta x1 = " + str(delta_x1) + ", ", "delta x2 = " + str(delta_x2) + ", ", "x10 = " + str(x10) + ", ",
      "x20 = " + str(x20) + ", ", "a0 = " + str(a_0) + ", ", "a1 = " + str(a_1), "a2 = " + str(a_2))
print("Натуралызоване рівняння")
print(str(a_0) + " " + str(a_1) + " * x1" + " " + str(a_2) + " * x2")
print("Перевірка по рядках")
print(a_0 + a_1 * (-30) + a_2 * (-15), "=", y_average_mass[0])
print(a_0 + a_1 * (-30) + a_2 * (35), "=", y_average_mass[1])
print(a_0 + a_1 * (0) + a_2 * (-15), "=", y_average_mass[2])
