import numpy as np
import copy
from scipy.stats import t as krit
import random


def cramer(matrix, koef):
    matrix = np.array(matrix)
    det = np.linalg.det(matrix)
    b = []
    for _ in range(8):
        dop = copy.deepcopy(matrix)
        dop[_] = koef
        det_dop = np.linalg.det(dop)
        b.append(det_dop / det)
    return b


def disp(y, y_avg):
    dispersion = []
    for _ in range(8):
        s = 0
        for j in range(3):
            s += (y[_][j] - y_avg[_]) ** 2
        dispersion.append(s / 3)
    return dispersion


def cohren(dispersion):
    Gp = max(dispersion) / sum(dispersion)
    Gt = 0.5157
    if Gp < Gt:
        print("За критерієм Кохрена дисперсія однорідна")


def fisher(y_new, y_avg):
    fisher = [4.5, 3.6, 3.2, 3, 2.9, 2.7, 2.4, 2.2, 2]
    d = len(eq)
    S_ad_1 = (3 / (8 - d))
    S_ad_2 = 0
    for i in range(8):
        S_ad_2 += (y_new[i] - y_avg[i]) ** 2
    S_ad = S_ad_1 * S_ad_2
    Fp = S_ad / S_B

    Ft = fisher[8 - d - 1]
    if Fp > Ft:
        print("Fp=", Fp, ">Fт=", Ft, "Отже, рівняння регресії неадекватно оригіналу при рівні значимості 0.05")


def stud(dispersion, y_avg, plan, b_mass, x_mass):
    global eq, S_B
    S_B = sum(dispersion) / len(dispersion)
    S_Bs = S_B / (3 * 8)
    SBs = np.sqrt(S_Bs)
    eq = []

    B0, B1, B2, B3, B4, B5, B6, B7 = 0, 0, 0, 0, 0, 0, 0, 0
    t = []
    t_ = []

    for _ in range(8):
        Suma = 0
        for j in range(8):
            Suma += y_avg[j] * plan[j][_]
        globals()["B" + str(_)] = Suma / 8

    for _ in range(8):
        t.append(globals()["B" + str(_)] / SBs)

    kriterii = krit.ppf(q=0.975, df=2 * 8)

    print("За критерієм Стьюдента")

    for _ in range(8):
        if t[_] < kriterii:
            print("Коефіцієнт", "b" + str(_), "приймаємо незначним")
        elif t[_] > kriterii:
            eq.append([b_mass[_], _])

    for i in range(8):
        S = 0
        for j in range(len(eq)):
            if eq[j][1] == 0:
                S += eq[j][0]
            else:
                tmp = eq[j][1] - 1
                S += eq[j][0] * x_mass[tmp][i]

        #массив t_ містить отримані функції відгуку зі спрощеними коефіцієнтами
        t_.append(S)
        #----------------------------------------------------------------------

    return t_


plan = [
    [+1, -1, -1, -1, +1, +1, +1, -1],
    [+1, -1, +1, +1, -1, -1, +1, -1],
    [+1, +1, -1, +1, -1, +1, -1, -1],
    [+1, +1, +1, -1, +1, -1, -1, -1],
    [+1, -1, -1, +1, +1, -1, -1, +1],
    [+1, -1, +1, -1, -1, +1, -1, +1],
    [+1, +1, -1, -1, -1, -1, +1, +1],
    [+1, +1, +1, +1, +1, +1, +1, +1]
]

x1_min, x1_max = -5, 15
x2_min, x2_max = 25, 45
x3_min, x3_max = 15, 45
x_min_avg, x_max_avg = 0, 0

for i in range(3):
    x_min_avg += globals()["x" + str(i + 1) + "_min"]
    x_max_avg += globals()["x" + str(i + 1) + "_max"]

x_min_avg, x_max_avg = round(x_min_avg / 3), round(x_max_avg / 3)
y_min_avg, y_max_avg = 200 + x_min_avg, 200 + x_max_avg

y_mass = [[random.randint(y_min_avg, y_max_avg + 1) for i in range(3)] for j in range(8)]

y_sum = [round(sum(i) / 3, 2) for i in y_mass]

N = 8
x_mass = np.array([
    [x1_min, x2_min, x3_min, x1_min * x2_min, x1_min * x3_min, x2_min * x3_min, x1_min * x2_min * x3_min],
    [x1_min, x2_max, x3_max, x1_min * x2_max, x1_min * x3_max, x2_max * x3_max, x1_min * x2_max * x3_max],
    [x1_max, x2_min, x3_max, x1_max * x2_min, x1_max * x3_max, x2_min * x3_max, x1_max * x2_min * x3_max],
    [x1_max, x2_max, x3_min, x1_max * x2_max, x1_max * x3_min, x2_max * x3_min, x1_max * x2_max * x3_min],
    [x1_min, x2_min, x3_max, x1_min * x2_min, x1_min * x3_max, x2_min * x3_max, x1_min * x2_min * x3_max],
    [x1_min, x2_max, x3_min, x1_min * x2_max, x1_min * x3_min, x2_max * x3_min, x1_min * x2_max * x3_min],
    [x1_max, x2_min, x3_min, x1_max * x2_min, x1_max * x3_min, x2_min * x3_min, x1_max * x2_min * x3_min],
    [x1_max, x2_max, x3_max, x1_max * x2_max, x1_max * x3_max, x2_max * x3_max, x1_max * x2_max * x3_max]
]).T

x1 = x_mass[0]
x2 = x_mass[1]
x3 = x_mass[2]
b_delta = [[N, sum(x1), sum(x2), sum(x3), sum(x1 * x2), sum(x1 * x3), sum(x2 * x3), sum(x1 * x2 * x3)],
           [sum(x1), sum(x1 ** 2), sum(x1 * x2), sum(x1 * x3), sum(x1 ** 2 * x2), sum(x1 ** 2 * x3),
            sum(x1 * x2 * x3), sum(x1 ** 2 * x2 * x3)],
           [sum(x2), sum(x1 * x2), sum(x2 ** 2), sum(x2 * x3), sum(x1 * x2 ** 2), sum(x1 * x2 * x3),
            sum(x2 ** 2 * x3), sum(x1 * x2 ** 2 * x3)],
           [sum(x3), sum(x1 * x3), sum(x2 * x3), sum(x3 ** 2), sum(x1 * x2 * x3), sum(x1 * x3 ** 2),
            sum(x2 * x3 ** 2), sum(x1 * x2 * x3 ** 2)],
           [sum(x1 * x2), sum(x1 ** 2 * x2), sum(x1 * x2 ** 2), sum(x1 * x2 * x3), sum(x1 ** 2 * x2 ** 2),
            sum(x1 ** 2 * x2 * x3), sum(x1 * x2 ** 2 * x3), sum(x1 ** 2 * x2 ** 2 * x3)],
           [sum(x1 * x3), sum(x1 ** 2 * x3), sum(x1 * x2 * x3), sum(x1 * x3 ** 2), sum(x1 ** 2 * x2 * x3),
            sum(x1 ** 2 * x3 ** 2), sum(x1 * x2 * x3 ** 2), sum(x1 ** 2 * x2 * x3 ** 2)],
           [sum(x2 * x3), sum(x1 * x2 * x3), sum(x2 ** 2 * x3), sum(x2 * x3 ** 2), sum(x1 * x2 ** 2 * x3),
            sum(x1 * x2 * x3 ** 2), sum(x2 ** 2 * x3 ** 2), sum(x1 * x2 ** 2 * x3 ** 2)],
           [sum(x1 * x2 * x3), sum(x1 ** 2 * x2 * x3), sum(x1 * x2 ** 2 * x3), sum(x1 * x2 * x3 ** 2),
            sum(x1 ** 2 * x2 ** 2 * x3), sum(x1 ** 2 * x2 * x3 ** 2), sum(x1 * x2 ** 2 * x3 ** 2),
            sum(x1 ** 2 * x2 ** 2 * x3 ** 2)]]
print("Матриця планування для m=3")
for i in range(3):
    print("Y" + str(i + 1), "=", np.array(y_mass).T[i])
print("Середні значення функції відгуку за рядками:")
print("Y_R = ", y_sum)
k = [sum(y_sum), sum(y_sum * x1), sum(y_sum * x2), sum(y_sum * x3), sum(y_sum * x2 * x1), sum(y_sum * x3 * x1),
     sum(y_sum * x2 * x3), sum(y_sum * x2 * x1 * x3)]
b_mass = cramer(b_delta, k)
print("Коефіцієнти рівняння регресії:")
for i in range(8):
    print("b" + str(i), "=", b_mass[i])
d = disp(y_mass, y_sum)
cohren(d)
new_t = stud(d, y_sum, plan, b_mass, x_mass)
print("Отримані функції відгуку зі спрощеними коефіцієнтами = ", new_t)
fisher(new_t, y_sum)
