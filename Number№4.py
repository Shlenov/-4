import json
import os
import scipy.constants as constants
import scipy.special as special
import numpy as np
import csv
import matplotlib.pyplot as plt
import wget

D = None
Fmin = None
Fmax = None

try:
    wget.download('https://jenyay.net/uploads/Student/Modelling/task_rcs.csv')
except:
    print(-1)
    
with open('task_rcs.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in spamreader:
        if row[0] == '11,':
            #print(D, Fmin, Fmax)
            print(row[0], row[1][0:], row[2][0:], row[3][0:])
            D = float(row[1][0:-1])
            Fmin = float(row[2][0:-1])
            Fmax = float(row[3][0:])
    
r = D/2

#значения разброса частот с шагом 10000000
freq_range = range(int(Fmin) , int(Fmax), 10000000)

# Класс отвечает за рассчет ЭПР,
class CalculateEDA:
    def __init__(self, r, frequency):
        # в конструкторе задаем начальные переменные рассчета
        self.radius = r
        self.wave_length = 0
        self.k = 0
        self._freq_range_ = frequency

    def a_n(self, n):
        numerator = np.longdouble(special.spherical_jn(n, self.k * self.radius))
        divider = self.h_n(n, self.k * self.radius)
        return np.divide(numerator, divider)

    def b_n(self, n):
        numerator = self.k * self.radius * np.longdouble(special.spherical_jn(n - 1, self.k * self.radius)) - n * np.longdouble(special.spherical_jn(n, self.k * self.radius))
        divider = self.k * self.radius * self.h_n(n - 1, self.k * self.radius) - n * self.h_n(n, self.k * self.radius)
        return np.divide(numerator, divider)

    def h_n(self, n, arg):
        return np.clongdouble(special.spherical_jn(n, arg) + 1j*special.spherical_yn(n, arg))

    def EDA(self):
        coef = self.wave_length**2 / constants.pi
        partForml = 0
        # оператор суммы в формуле c верхним пределом 50
        for n in range(1, 50):
            partForml += (-1)**n * (n + 0.5) * (self.b_n(n) - self.a_n(n))
        result = coef * abs(partForml) ** 2
        return result

    # функция запускает процесс просчёта для введёных границ частоты
    def calculateData(self):
        self.data = []
        self.Freq = []
        self.Lambda = []
        self.Rcs = []
        for freq in self._freq_range_:
            # обновляем длину волны и волновое число
            self.wave_length = np.longdouble(constants.c / freq)
            self.k = np.longdouble(2 * constants.pi / self.wave_length)
            # получаем значение ЭПР для новых параметров
            temp_eda = self.EDA()
            self.Freq.append(float(freq))
            self.Lambda.append(float(self.wave_length))
            self.Rcs.append(float(temp_eda))
        #self.data.append({"freq": self.Freq, "lambda": self.Lambda, "rcs": self.Rcs})
        self.data = [self.Freq, self.Lambda, self.Rcs]
        return self.data

class Output:
    def __init__(self, data):
        # в данном случае в конструкторе нас интересует только передача массива с данными о точках
        self.data = data
    # функция для создания графика
    def drawPlotData(self):
        freq = self.data[0]
        rcs = self.data[2]
        plt.plot(freq, rcs)
        plt.ylabel('ЭПР')
        plt.xlabel('Частота')
        plt.title('зависимость ЭПР от частоты')
        plt.grid()
        plt.show()
    def saveJson(self, filename):
        with open(filename, 'w') as f:
            json.dump({"freq": self.data[0],"lambda": self.data[1], "rcs": self.data[2]}, f, indent=4)


# Создаем объект рассчёта
calculator = CalculateEDA(r, freq_range)
data = calculator.calculateData()

# Создаём объект отвечающий за вывод
output = Output(data)
output.saveJson('result.json')
output.drawPlotData()