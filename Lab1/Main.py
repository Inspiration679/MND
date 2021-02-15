import random
a0,a1,a2,a3=0,0,0,0
for i in range(4):
    globals()['a' + str(i)]=random.randint(50,80)


class Main:
    def __init__(self, m, n, a0, a1, a2, a3):
        self.a0, self.a1, self.a2, self.a3, self.row, self.column = a0, a1, a2, a3, m, n
        self.main_func()

    def main_func(self):
        self.make_mass(self.row, self.column)
        self.y_d0_x0()
        self.make_second_mass()
        self.show()

    def make_mass(self, size1, size2):
        self.main_mass = []
        for i in range(size1):
            temp_mass = []
            for j in range(size2):
                temp_mass.append(float(random.randint(0, 20)))
            self.main_mass.append(temp_mass)

    def make_second_mass(self):
        self.second_mass = []
        for i in range(self.row):
            temp_mass = []
            for j in range(self.column):
                temp_mass.append(round((self.main_mass[i][j] - self.x0_mass[j]) / self.dx_mass[j], 2))
            self.second_mass.append(temp_mass)

    def y_d0_x0(self):
        self.dx_mass, self.x0_mass, self.y_mass, self.mass, self.average_y, self.new_mass,self.y_x0 = [], [], [], ["Num", "X1",
                                                                                                         "X2", "X3",
                                                                                                         "Y", "Xn1",
                                                                                                         "Xn2",
                                                                                                         "Xn3"], 0, [],self.a0
        for i in range(self.row):
            self.y_mass.append(
                self.a0 + self.main_mass[i][0] * self.a1 + self.main_mass[i][1] * self.a2 + self.main_mass[i][
                    2] * self.a3)

        for i in range(self.column):
            temp_mass = []
            for j in range(self.row):
                temp_mass.append(self.main_mass[j][i])
            self.x0_mass.append(((max(temp_mass) + min(temp_mass)) / 2))

        for i in range(self.column):
            temp_mass = []
            for j in range(self.row):
                temp_mass.append(self.main_mass[j][i])
            self.dx_mass.append(self.x0_mass[i] - min(temp_mass))

        for i in range(len(self.y_mass)):
            self.average_y += self.y_mass[i]

        self.average_y = self.average_y / len(self.y_mass)

        for i in range(len(self.y_mass)):
            if self.y_mass[i] < self.average_y:
                self.new_mass.append(self.average_y - self.y_mass[i])

        self.var = self.average_y - min(self.new_mass)

        for i in range(3):
            self.y_x0+=globals()['a' + str(i+1)]*self.x0_mass[i]
        self.dx_mass.append(self.average_y)
        self.x0_mass.append(round(self.y_x0,2))


        for i in range(3):
            self.dx_mass.append("---")
            self.x0_mass.append("---")

    def show(self):
        print("{:<2}".format(self.mass[0]), "{:>3}".format(self.mass[1]), "{:>6}".format(self.mass[2]),
              "{:>6.5}".format(self.mass[3]), "{:>5.5}".format(self.mass[4]), "{:>7}".format(self.mass[5]),
              "{:>7}".format(self.mass[6]), "{:>5.5}".format(self.mass[7]))
        for i in range(self.row):
            print("{:<3}".format(i + 1), end=" ")
            for j in range(self.column):
                print("{:<6}".format(self.main_mass[i][j]), end=" ")
            print("{:<6}".format(self.y_mass[i]), end=" ")
            for k in range(self.column):
                print("{:<6}".format(self.second_mass[i][k]), end=" ")
            print("|")

        print("{:<3}".format("x0"), end=" ")
        for i in self.x0_mass:
            print("{:<6}".format(i), end=" ")
        print("{:<5}".format("|\ndx"), end=" ")
        for i in self.dx_mass:
            print("{:<6}".format(i), end=" ")
        print("|")
        for i in range(4):
            print("a" + str(i) + " = ", globals()['a' + str(i)], end=", ")
        print("Середнє Y - ", self.average_y, end=", ")
        print("> Y - " + str(self.var) + ".")


if __name__ == "__main__":
    t = Main(8, 3, a0,a1,a2,a3)
