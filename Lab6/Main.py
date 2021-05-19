from copy import deepcopy
from math import sqrt
from random import random
from time import process_time as clock
import numpy as np
from prettytable import PrettyTable
x1_min = -30
x1_max = 0
x2_min = -15
x2_max = 35
x3_min = -30
x3_max = 35
t1,t2,t3=0,0,0
koefs = [9.6, 8.8, 9.2, 6.9, 7.5, 0.7, 6.9, 7.4, 0.1, 5.9, 10.0]  # Koefs for search y
x_average_max = (x1_max + x2_max + x3_max) / 3
x_average_min = (x1_min + x2_min + x3_min) / 3
y_max = 200 + x_average_max
y_min = 200 + x_average_min
T1,T2,T3=0,0,0

def replace_column(list_: list, column, list_replace):
    list_ = deepcopy(list_)
    for i in range(len(list_)):
        list_[i][column] = list_replace[i]
    return list_


def append_to_list_x(x: list, variant: int):
    if variant == 1:
        for i in range(len(x)):
            x[i].append(x[i][1] * x[i][2])
            x[i].append(x[i][1] * x[i][3])
            x[i].append(x[i][2] * x[i][3])
            x[i].append(x[i][1] * x[i][2] * x[i][3])
    if variant == 2:
        for i in range(len(x)):
            x[i].append(x[i][1] * x[i][2])
            x[i].append(x[i][1] * x[i][3])
            x[i].append(x[i][2] * x[i][3])
            x[i].append(x[i][1] * x[i][2] * x[i][3])
            x[i].append(x[i][1] * x[i][1])
            x[i].append(x[i][2] * x[i][2])
            x[i].append(x[i][3] * x[i][3])
    for i in range(len(x)):
        for j in range(len(x[i])):
            if round(x[i][j], 3) == 0:
                x[i][j] = 0
            x[i][j] = round(x[i][j], 3)


def get_value(table: dict, key: int):
    value = table.get(key)
    if value is not None:
        return value
    for i in table:
        if type(i) == range and key in i:
            return table.get(i)


def main(m, n):
    global T1,T2,T3
    if n == 14:
        const_l = 1.73
        print(
            'ŷ = b0 + b1 * x1 + b2 * x2 + b3 * x3 + b12 * x1 * x2 + b13 * x1 * x3 + b23 * x2 * x3 + b123 * x1 * x2 * '
            'x3 + b11 * x1 * x1 + b22 * x2 * x2 + b33 * x3 * x3')
        norm_x = [
            [+1, -1, -1, -1],
            [+1, -1, +1, +1],
            [+1, +1, -1, +1],
            [+1, +1, +1, -1],
            [+1, -1, -1, +1],
            [+1, -1, +1, -1],
            [+1, +1, -1, -1],
            [+1, +1, +1, +1],
            [+1, -const_l, 0, 0],
            [+1, const_l, 0, 0],
            [+1, 0, -const_l, 0],
            [+1, 0, const_l, 0],
            [+1, 0, 0, -const_l],
            [+1, 0, 0, const_l],
        ]

        delta_x1 = (x1_max - x1_min) / 2
        delta_x2 = (x2_max - x2_min) / 2
        delta_x3 = (x2_max - x3_min) / 2
        x01 = (x1_min + x1_max) / 2
        x02 = (x2_min + x2_max) / 2
        x03 = (x3_min + x3_max) / 2

        x = [
            [1, x1_min, x2_min, x3_min],
            [1, x1_min, x2_max, x3_max],
            [1, x1_max, x2_min, x3_max],
            [1, x1_max, x2_max, x3_min],
            [1, x1_min, x2_min, x3_max],
            [1, x1_min, x2_max, x3_min],
            [1, x1_max, x2_min, x3_min],
            [1, x1_max, x2_max, x3_max],
            [1, -const_l * delta_x1 + x01, x02, x03],
            [1, const_l * delta_x1 + x01, x02, x03],
            [1, x01, -const_l * delta_x2 + x02, x03],
            [1, x01, const_l * delta_x2 + x02, x03],
            [1, x01, x02, -const_l * delta_x3 + x03],
            [1, x01, x02, const_l * delta_x3 + x03],
        ]

        append_to_list_x(norm_x, variant=2)
        append_to_list_x(x, variant=2)

    if n == 8:
        print(
            'ŷ = b0 + b1 * x1 + b2 * x2 + b3 * x3 + b12 * x1 * x2 + b13 * x1 * x3 + b23 * x2 * x3 + b123 * x1 * x2 * x3'
        )
        norm_x = [
            [+1, -1, -1, -1],
            [+1, -1, +1, +1],
            [+1, +1, -1, +1],
            [+1, +1, +1, -1],
            [+1, -1, -1, +1],
            [+1, -1, +1, -1],
            [+1, +1, -1, -1],
            [+1, +1, +1, +1]
        ]

        x = [
            [1, x1_min, x2_min, x3_min],
            [1, x1_min, x2_max, x3_max],
            [1, x1_max, x2_min, x3_max],
            [1, x1_max, x2_max, x3_min],
            [1, x1_min, x2_min, x3_max],
            [1, x1_min, x2_max, x3_min],
            [1, x1_max, x2_min, x3_min],
            [1, x1_max, x2_max, x3_max]
        ]

        append_to_list_x(norm_x, variant=1)
        append_to_list_x(x, variant=1)

    if n == 4:
        print('ŷ = b0 + b1 * x1 + b2 * x2 + b3 * x3')
        norm_x = [
            [+1, -1, -1, -1],
            [+1, -1, +1, +1],
            [+1, +1, -1, +1],
            [+1, +1, +1, -1],
        ]
        x = [
            [1, x1_min, x2_min, x3_min],
            [1, x1_min, x2_max, x3_max],
            [1, x1_max, x2_min, x3_max],
            [1, x1_max, x2_max, x3_min],
        ]
    if n == 14:
        y = [[round(sum([koefs[j] * i[j] for j in range(len(koefs))]) + random() * 10 - 5, 3) for k in range(m)] for i
             in x]
    else:
        y = np.random.randint(y_min, y_max, size=(n, m))
    # y = np.random.randint(y_min, y_max, size=(n, m))
    y_av = list(np.average(y, axis=1))

    for i in range(len(y_av)):
        y_av[i] = round(y_av[i], 3)

    if n == 14:
        t = PrettyTable(['N', 'norm_x_0', 'norm_x_1', 'norm_x_2', 'norm_x_3', 'norm_x_1_x_2', 'norm_x_1_x_3',
                         'norm_x_2_x_3', 'norm_x_1_x_2_x_3', 'norm_x_1_x_1', 'norm_x_2_x_2', 'norm_x_3_x_3', 'x_0',
                         'x_1', 'x_2', 'x_3', 'x_1_x_2', 'x_1_x_3', 'x_2_x_3', 'x_1_x_2_x_3', 'x_1_x_1', 'x_2_x_2',
                         'x_3_x_3'] + [f'y_{i + 1}' for i in range(m)] + ['y_av'])

    if n == 8:
        t = PrettyTable(['N', 'norm_x_0', 'norm_x_1', 'norm_x_2', 'norm_x_3', 'norm_x_1_x_2', 'norm_x_1_x_3',
                         'norm_x_2_x_3', 'norm_x_1_x_2_x_3', 'x_0', 'x_1', 'x_2', 'x_3', 'x_1_x_2', 'x_1_x_3',
                         'x_2_x_3', 'x_1_x_2_x_3'] + [f'y_{i + 1}' for i in range(m)] + ['y_av'])
    if n == 4:
        t = PrettyTable(
            ['N', 'norm_x_0', 'norm_x_1', 'norm_x_2', 'norm_x_3', 'x_0', 'x_1', 'x_2', 'x_3'] +
            [f'y_{i + 1}' for i in range(m)] + ['y_av'])

    for i in range(n):
        t.add_row([i + 1] + list(norm_x[i]) + list(x[i]) + list(y[i]) + [y_av[i]])
    print(t)

    m_ij = []
    for i in range(len(x[0])):
        m_ij.append([round(sum([x[k][i] * x[k][j] for k in range(len(x))]) / 14, 3) for j in range(len(x[i]))])

    k_i = []
    for i in range(len(x[0])):
        a = sum(y_av[j] * x[j][i] for j in range(len(x))) / 14
        k_i.append(a)

    det = np.linalg.det(m_ij)
    det_i = [np.linalg.det(replace_column(m_ij, i, k_i)) for i in range(len(k_i))]

    b_i = [round(i / det, 3) for i in det_i]
    if n == 14:
        print(
            f"\nThe naturalized regression equation: "
            f"y = {b_i[0]:.5f} + {b_i[1]:.5f} * x1 + {b_i[2]:.5f} * x2 + "
            f"{b_i[3]:.5f} * x3 + {b_i[4]:.5f} * x1 * x2 + "
            f"{b_i[5]:.5f} * x1 * x3 + {b_i[6]:.5f} * x2 * x3 + {b_i[7]:.5f} * x1 * x2 * x3 + {b_i[8]:.5f} * x1 * x1 + "
            f"{b_i[9]:.5f} * x2 * x2 + {b_i[10]:.5f} * x3 * x3")
    if n == 8:
        print(
            f"\nThe naturalized regression equation: "
            f"y = {b_i[0]:.5f} + {b_i[1]:.5f} * x1 + {b_i[2]:.5f} * x2 + "
            f"{b_i[3]:.5f} * x3 + {b_i[4]:.5f} * x1 * x2 + "
            f"{b_i[5]:.5f} * x1 * x3 + {b_i[6]:.5f} * x2 * x3 + {b_i[7]:.5f} * x1 * x2 * x3")
    if n == 4:
        print(
            f"\nThe naturalized regression equation: "
            f"y = {b_i[0]:.5f} + {b_i[1]:.5f} * x1 + {b_i[2]:.5f} * x2 + {b_i[3]:.5f} * x3\n")

    check_i = [round(sum(b_i[j] * i[j] for j in range(len(b_i))), 3) for i in x]
    for i in range(len(check_i)):
        print(f'ŷ{i + 1} = {check_i[i]}, y_av{i + 1} = {y_av[i]}')

    print("\n[ Kohren's test ]")
    f_1 = m - 1
    f_2 = n
    s_i = [sum([(i - y_av[j]) ** 2 for i in y[j]]) / m for j in range(len(y))]
    g_p = max(s_i) / sum(s_i)

    table = {2: 0.75, 3: 0.6841, 4: 0.6287, 5: 0.5892, 6: 0.5598, 7: 0.5365, 8: 0.5175, 9: 0.5017, 10: 0.4884,
             range(11, 17): 0.4366, range(17, 37): 0.3720, range(37, 2 ** 100): 0.3093}
    g_t = get_value(table, m)

    if g_p < g_t:
        print(f"The variance is homogeneous: Gp = {g_p:.5} < Gt = {g_t}")
    else:
        print(f"The variance is not homogeneous Gp = {g_p:.5} > Gt = {g_t}\nStart again with m = m + 1 = {m + 1}")
        return main(m=m + 1, n=n)
    if timer==9:
        t1_1=clock()
        T1=t1_1-t1

    print("\n[ Student's test ]")
    s2_b = sum(s_i) / n
    s2_beta_s = s2_b / (n * m)
    s_beta_s = sqrt(s2_beta_s)
    beta_i = [sum([norm_x[i][j] * y_av[i] for i in range(len(norm_x))]) / n for j in range(len(norm_x[0]))]
    beta_i = [round(i, 3) for i in beta_i]

    t = [abs(i) / s_beta_s for i in beta_i]
    if n == 14:
        beta_i = b_i
    f_3 = f_1 * f_2
    t_table = {4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160,
               14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086, 21: 2.08, 22: 2.074,
               23: 2.069, 24: 2.064, range(25, 30): 2.06, range(30, 40): 2.042, range(40, 60): 2.021, range(60, 100): 2,
               range(100, 2 ** 100): 1.96}
    d = deepcopy(len(beta_i))
    for i in range(len(t)):
        if get_value(t_table, f_3) > t[i]:
            beta_i[i] = 0
            d -= 1
    if n == d:
        n = 8 if n == 4 else 14
        print(f"n=d\nStart again with n = {n} and m = {m}")
        return main(m=m, n=n)
    if n == 14:
        print(
            f"\nThe naturalized simplified regression equation: "
            f"y = {beta_i[0]:.5f} + {beta_i[1]:.5f} * x1 + "
            f"{beta_i[2]:.5f} * x2 + {beta_i[3]:.5f} * x3 + {beta_i[4]:.5f} * x1 * x2 + "
            f"{beta_i[5]:.5f} * x1 * x3 + {beta_i[6]:.5f} * x2 * x3 + {beta_i[7]:.5f} * x1 * x2 * x3 + "
            f"{beta_i[8]:.5f} * x1 * x1 + {beta_i[9]:.5f} * x2 * x2 + {beta_i[10]:.5f} * x3 * x3")
        check_i = [round(sum(beta_i[j] * i[j] for j in range(len(beta_i))), 3) for i in x]

    if n == 8:
        print(
            f"\nThe normalized regression equation: "
            f"y = {beta_i[0]:.5f} + {beta_i[1]:.5f} * x1 + {beta_i[2]:.5f} * x2 + "
            f"{beta_i[3]:.5f} * x3 + {beta_i[4]:.5f} * x1 * x2 + "
            f"{beta_i[5]:.5f} * x1 * x3 + {beta_i[6]:.5f} * x2 * x3 + {beta_i[7]:.5f} * x1 * x2 * x3")
        check_i = [round(sum(beta_i[j] * i[j] for j in range(len(beta_i))), 3) for i in norm_x]

    if n == 4:
        print(
            f"\nThe normalized regression equation: "
            f"y = {beta_i[0]:.5f} + {beta_i[1]:.5f} * x1 + {beta_i[2]:.5f} * x2 + "
            f"{beta_i[3]:.5f} * x3")
        check_i = [round(sum(beta_i[j] * i[j] for j in range(len(beta_i))), 3) for i in norm_x]

    for i in range(len(check_i)):
        print(f'ŷ{i + 1} = {check_i[i]}, y_av{i + 1} R {y_av[i]}')
    if timer==9:
        t2_2=clock()
        T2=t2_2-t2

    print("\n[ Fisher's test ]")
    f_4 = n - d
    s2_ad = m / f_4 * sum([(check_i[i] - y_av[i]) ** 2 for i in range(len(y_av))])
    f_p = s2_ad / s2_b
    f_t = {
        1: [164.4, 199.5, 215.7, 224.6, 230.2, 234, 235.8, 237.6],
        2: [18.5, 19.2, 19.2, 19.3, 19.3, 19.3, 19.4, 19.4],
        3: [10.1, 9.6, 9.3, 9.1, 9, 8.9, 8.8, 8.8],
        4: [7.7, 6.9, 6.6, 6.4, 6.3, 6.2, 6.1, 6.1],
        5: [6.6, 5.8, 5.4, 5.2, 5.1, 5, 4.9, 4.9],
        6: [6, 5.1, 4.8, 4.5, 4.4, 4.3, 4.2, 4.2],
        7: [5.5, 4.7, 4.4, 4.1, 4, 3.9, 3.8, 3.8],
        8: [5.3, 4.5, 4.1, 3.8, 3.7, 3.6, 3.5, 3.5],
        9: [5.1, 4.3, 3.9, 3.6, 3.5, 3.4, 3.3, 3.3],
        10: [5, 4.1, 3.7, 3.5, 3.3, 3.2, 3.1, 3.1],
        11: [4.8, 4, 3.6, 3.4, 3.2, 3.1, 3, 3],
        12: [4.8, 3.9, 3.5, 3.3, 3.1, 3, 2.9, 2.9],
        13: [4.7, 3.8, 3.4, 3.2, 3, 2.9, 2.8, 2.8],
        14: [4.6, 3.7, 3.3, 3.1, 3, 2.9, 2.8, 2.7],
        15: [4.5, 3.7, 3.3, 3.1, 2.9, 2.8, 2.7, 2.7, 2.7, 2.7, 2.6, 2.6],
        16: [4.5, 3.6, 3.2, 3, 2.9, 2.7, 2.6, 2.6],
        17: [4.5, 3.6, 3.2, 3, 2.8, 2.7, 2.5, 2.3],
        18: [4.4, 3.6, 3.2, 2.9, 2.8, 2.7, 2.5, 2.3],
        19: [4.4, 3.5, 3.1, 2.9, 2.7, 2.7, 2.4, 2.3],
        range(20, 22): [4.4, 3.5, 3.1, 2.8, 2.7, 2.7, 2.4, 2.3],
        range(22, 24): [4.3, 3.4, 3.1, 2.8, 2.7, 2.6, 2.4, 2.3],
        range(24, 26): [4.3, 3.4, 3, 2.8, 2.6, 2.5, 2.3, 2.2],
        range(26, 28): [4.2, 3.4, 3, 2.7, 2.6, 2.5, 2.3, 2.2],
        range(28, 30): [4.2, 3.3, 3, 2.7, 2.6, 2.4, 2.3, 2.1],
        range(30, 40): [4.2, 3.3, 3, 2.7, 2.6, 2.4, 2.3, 2.1, 2, 2, 2, 2],
        range(40, 60): [4.1, 3.2, 2.9, 2.6, 2.5, 2.3, 2.2, 2, 1.9, 1.9, 1.9, 1.9],
        range(60, 120): [4, 3.2, 2.8, 2.5, 2.4, 2.3, 2.1, 1.9, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8],
        range(120, 2 ** 100): [3.8, 3, 2.6, 2.4, 2.2, 2.1, 2, 2, 1.9, 1.9, 1.9, 1.8, 1.8]
    }
    if f_p > get_value(f_t, f_3)[f_4]:
        n = 8 if n == 4 else 14
        print(
            f"fp = {f_p} > ft = {get_value(f_t, f_3)[f_4]}.\n"
            f"The mathematical model is not adequate to the experimental data\n"
            f"Start again with m = {m} and n = {n}")
        return False
    else:
        print(
            f"fP = {f_p} < fT = {get_value(f_t, f_3)[f_4]}.\n"
            f"The mathematical model is adequate to the experimental data\n")
        if timer==9:
            t3_3=clock()
            T3=t3_3-t3
        return True

timer=0
# n = 14 because if you start with 4 then it will not reach 14
t1=clock()
for i in range(10):
    t5=False
    while t5 ==False:
        t5 = main(m=2, n=14)
    timer+=1
print("Кохрен - ",T1/10,"Стьюдент - ",T2/10,"Фішер - ",T3/10)