'''
Решение уравнения Бюргерса по схеме Лакса-Вендроффа
'''

import matplotlib.pyplot as plt
import numpy as np
import math


# initial conditions
def initial_conditions(U, X, NX):
    # sinusoidal signal
    for i in range(1, NX):
        U[i] = C0 + C1 * math.sin(K * X[i])

    # impulse signal
'''
    for i in range(1, NX // 5 + 1):
        U[i] = np.array(0.5)
    for i in range(NX // 5, 2 * NX // 5 + 1):
        U[i] = np.array(0.0)
    for i in range(2 * NX // 5, 3 * NX // 5 + 1):
        U[i] = np.array(0.5)
    for i in range(3 * NX // 5, 4 * NX // 5 + 1):
        U[i] = np.array(-0.5)
    for i in range(4 * NX // 5, NX + 1):
        U[i] = np.array(0.5)
'''


# boundary conditions
def boundary_conditions(U, NX):
    U[0] = U[NX - 1]
    U[NX] = U[1]


# saving results
def output_data(output, X, UN, NX):
    for i in range(0, NX + 1):
        output.write(str(X[i]) + ', ' + str(UN[i]) + '\n')


with open('input_4.txt', 'r') as data:
    L = float(data.readline().split('#')[0])  # длина расчетной области
    M = float(data.readline().split('#')[0])  # номер гармоники
    C = float(data.readline().split('#')[0])  # коэффициент переноса
    C0 = int(data.readline().split('#')[0])  # коэфф, задающий начальные усл-я
    C1 = int(data.readline().split('#')[0])  # коэфф, задающий начальные усл-я
    NX = int(data.readline().split('#')[0])  # число узлов по координате
    NT = int(data.readline().split('#')[0])  # число узлов по времени
    CFL = float(data.readline().split('#')[0])  # число Куранта

K = M * math.pi / L

dx = L / (NX - 1)

dt = CFL * dx / C

Time = dt * NT

X = np.array([0.0 for i in range(NX + 1)])

U = np.array([0.0 for i in range(NX + 1)])  # characteristic form

U2 = np.array([0.0 for i in range(NX + 1)])  # divergent form

UN = np.array([0.0 for i in range(NX + 1)])  # characteristic form

UN2 = np.array([0.0 for i in range(NX + 1)])  # divergent form

Z = np.array([0.0 for i in range(NX + 1)])

g = np.array([0.0 for i in range(NX + 1)])

U0 = np.array([0.0 for i in range(NX + 1)])

# nodes by coordinate
X[0] = - dx / 2.0
for i in range(1, NX):
    X[i] = X[i - 1] + dx
X[NX] = L + dx / 2.0

initial_conditions(U, X, NX)

boundary_conditions(U, NX)

initial_conditions(U2, X, NX)

boundary_conditions(U2, NX)

# initial signal
for i in range(0, NX + 1):
    U0[i] = U[i]

# solution of the Burgers equation (characteristic form)
for j in range(1, NT + 1):
    for i in range(0, NX):
        UN[i] = (U[i] - (dt * U[i]) * (U[i + 1] - U[i - 1]) / (2 * dx) + (dt * U[i] / dx) ** 2 * (U[i + 1] - 2 * U[i] +
                                                                                                  U[i - 1]) / 2)
    boundary_conditions(UN, NX)
    U[:] = UN[:]

# solution of the Burgers equation (divergent form)
for j in range(1, NT + 1):
    for i in range(0, NX):
         UN2[i] = (U2[i] - dt * ((U2[i + 1]) ** 2 - (U2[i - 1]) ** 2) / (4 * dx) + (dt / dx) ** 2 * ((U2[i + 1]) ** 2 -
         2 * (U2[i]) ** 2 + (U2[i - 1]) ** 2) / 4)
    boundary_conditions(UN2, NX)
    U2[:] = UN2[:]


# saving results
output = open('output_4.txt', 'w')
output_data(output, X, UN, NX)
output.close()

# drawing graph
plt.rcParams['font.family'] = 'Courier New'

plt.rcParams['font.size'] = '24'

plt.xlabel('x, m', fontweight='bold')

plt.ylabel('u, m/s', fontweight='bold')

plt.plot(X[1:NX], UN[1:NX], 'red', linewidth = 4, label = 'characteristic form')  # characteristic form

plt.plot(X[1:NX], U0[1:NX], 'green', linestyle = '-.', linewidth = 3, label = 'initial signal')  # initial signal

plt.plot(X[1:NX], UN2[1:NX], 'black',  linewidth = 4, label = 'divergent form')  # divergent form

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17), shadow=False, edgecolor='black', ncol=2, fontsize='small')

plt.xlim([0, 1])

plt.ylim([-1, 1])

plt.grid()

plt.show()

np.delete(X, 0)
np.delete(U, 0)
np.delete(UN, 0)


