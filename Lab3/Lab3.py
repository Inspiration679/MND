import numpy as np
import random


def mk_mat(x1, x2):
    tmp_mass = []
    for i in range(4):
        tmp_mass.append(random.randint(x1, x2 + 1))
    return tmp_mass


def sum_riv(x1, x2, x3):
    global b0, b1, b2, b3
    y = b0 + b1 * x1 + b2 * x2 + b3 * x3
    return y


x1_min, x1_max = -30, 0
x2_min, x2_max = -15, 35
x3_min, x3_max = -30, 35
y_min, y_max = int(200 + (x1_min + x2_min + x3_min) / 3), int(200 + (x1_max + x2_max + x3_max) / 3)

plan_matrix = [[1, -1, -1, -1],
               [1, -1, 1, 1],
               [1, 1, -1, 1],
               [1, 1, 1, -1]]
x0_plan = [1, 1, 1, 1]
x1_plan = mk_mat(x1_min, x1_max)
x2_plan = mk_mat(x2_min, x2_max)
x3_plan = mk_mat(x3_min, x3_max)
y1_plan = mk_mat(y_min, y_max)
y2_plan = mk_mat(y_min, y_max)
y3_plan = mk_mat(y_min, y_max)

y_avg = []
mx = []
a = []
a_1 = []
my = 0

for i in range(4):
    t = 0
    for j in range(3):
        t += globals()["y" + str(j + 1) + "_plan"][i]
    y_avg.append(t / 3)

ya1, ya2, ya3, ya4 = y_avg[0], y_avg[1], y_avg[2], y_avg[3]

for i in range(3):
    t = 0
    for j in range(4):
        t += globals()["x" + str(i + 1) + "_plan"][j]
    mx.append(t / 4)

mx1, mx2, mx3 = mx[0], mx[1], mx[2]

for i in range(4):
    my += y_avg[i]

my = my / 4

for i in range(3):
    t = 0
    for j in range(4):
        t += globals()["x" + str(i + 1) + "_plan"][j] * y_avg[j]
    a.append(t / 4)

a1, a2, a3 = a[0], a[1], a[2]

for i in range(3):
    t = 0
    for j in range(4):
        t += globals()["x" + str(i + 1) + "_plan"][j] * globals()["x" + str(i + 1) + "_plan"][j]
    a_1.append(t / 4)

a11, a22, a33 = a_1[0], a_1[1], a_1[2]

a12 = a21 = (x1_plan[0] * x2_plan[0] + x1_plan[1] * x2_plan[1] + x1_plan[2] * x2_plan[2] + x1_plan[3] * x2_plan[3]) / 4
a13 = a31 = (x1_plan[0] * x3_plan[0] + x1_plan[1] * x3_plan[1] + x1_plan[2] * x3_plan[2] + x1_plan[3] * x3_plan[3]) / 4
a23 = a32 = (x2_plan[0] * x3_plan[0] + x2_plan[1] * x3_plan[1] + x2_plan[2] * x3_plan[2] + x2_plan[3] * x3_plan[3]) / 4

b0 = np.linalg.det([[my, mx1, mx2, mx3], [a1, a11, a12, a13], [a2, a12, a22, a32],
                    [a3, a13, a23, a33]]) / np.linalg.det(
    [[1, mx1, mx2, mx3], [mx1, a11, a12, a13], [mx2, a12, a22, a32], [mx3, a13, a23, a33]])

b1 = np.linalg.det([[1, my, mx2, mx3], [mx1, a1, a12, a13], [mx2, a2, a22, a32],
                    [mx3, a3, a23, a33]]) / np.linalg.det(
    [[1, mx1, mx2, mx3], [mx1, a11, a12, a13], [mx2, a12, a22, a32], [mx3, a13, a23, a33]])

b2 = np.linalg.det([[1, mx1, my, mx3], [mx1, a11, a1, a13], [mx2, a12, a2, a32],
                    [mx3, a13, a3, a33]]) / np.linalg.det(
    [[1, mx1, mx2, mx3], [mx1, a11, a12, a13], [mx2, a12, a22, a32], [mx3, a13, a23, a33]])

b3 = np.linalg.det([[1, mx1, mx2, my], [mx1, a11, a12, a1], [mx2, a12, a22, a2],
                    [mx3, a13, a23, a3]]) / np.linalg.det(
    [[1, mx1, mx2, mx3], [mx1, a11, a12, a13], [mx2, a12, a22, a32], [mx3, a13, a23, a33]])

s_mass = []

for i in range(4):
    t = 0
    for j in range(3):
        t += (globals()["y" + str(j + 1) + "_plan"][i] - y_avg[i]) ** 2
    s_mass.append(t / 3)

Gp = max(s_mass) / sum(s_mass)

S_B = sum(s_mass) / 4
S_Bs = S_B / 12
SBs = np.sqrt(S_Bs)
B = []
T = []

for i in range(4):
    t = 0
    for j in range(4):
        t += y_avg[j] * plan_matrix[j][i]
    B.append(t / 4)

for i in range(4):
    T.append(abs(B[i]) / SBs)

koef_T = []
new_T = []

for i in range(len(T)):
    if T[i] > 2.306:
        koef_T.append(i)

for i in range(4):
    t = 0
    for j in koef_T:
        t += globals()["b" + str(j)] * globals()["x" + str(j) + "_plan"][i]
    new_T.append(t)

Sad = 0

for i in range(4):
    Sad += (new_T[i] - y_avg[i]) ** 2

S_ad = Sad * 1.5
Fp = S_ad / S_Bs
Ft = 4.5

print("Матриця планування", "-", plan_matrix)
for i in range(3):
    print("x" + str(i + 1), "-", globals()["x" + str(i + 1) + "_plan"])
for i in range(3):
    print("y" + str(i + 1), "-", globals()["y" + str(i + 1) + "_plan"])
print("Середні значення у - ", y_avg)
print("mx1, mx2, mx3 - ", mx)
print("my - ", my)
print("a1, a2, a3 - ", a)
print("a11, a22, a33 - ", a11, a22, a33)
print("a12, a13, a23, a21, a31, a32 - ", a12, a13, a23, a21, a31, a32)
print("b0, b1, b2, b3 - ", b0, b1, b2, b3)
print("Отримане рівняння регресії: ", b0, "+", b1, "*", "x1", "+", b2, "*", "x2", "+", b3, "*", "x3", )
print(sum_riv(x1_plan[0], x2_plan[0], x3_plan[0]), "=", ya1)
print(sum_riv(x1_plan[1], x2_plan[1], x3_plan[1]), "=", ya2)
print(sum_riv(x1_plan[2], x2_plan[2], x3_plan[2]), "=", ya3)
print(sum_riv(x1_plan[3], x2_plan[3], x3_plan[3]), "=", ya4)
print("Дисперсія - ", s_mass)
print("Gp - ", Gp)
if Gp < 0.7679:
    print("Дисперсія однорідна")
else:
    print("Дисперсія неоднорідна")
print("S^2b - ", S_B)
print("S^2bs - ", S_Bs)
print("Sbs - ", SBs)
print("B0, B1, B2, B3 - ", B)
for i in range(len(T)):
    if T[i] > 2.306:
        print(T[i], "Входить в рівняння")
    else:
        print(T[i], "Виключається з рівнняня")
print("Рівнняня")
for i in koef_T:
    print("b" + str(i), "*", "x" + str(i), end=" ")
    print("+", end=" ")
print(0)
print("y1^, y2^, y3^, y4^ - ", new_T)
print("Ft - ", Ft)
print("Fp - ", Fp)
print("Sad - ", Sad)
if Fp > Ft:
    print("Fp>Fт.Отже, рівняння регресії неадекватно оригіналу при рівні значимості 0.05")
