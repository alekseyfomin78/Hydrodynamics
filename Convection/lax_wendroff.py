'''
Решение задачи конвекции по схеме Лакса-Вендроффа
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


with open('input_2.txt', 'r') as data:
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

U = np.array([0.0 for i in range(NX + 1)])

UN = np.array([0.0 for i in range(NX + 1)])

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

# initial signal
for i in range(0, NX + 1):
    U0[i] = U[i]

# numerical solution
for j in range(1, NT + 1):
    for i in range(0, NX):
        UN[i] = U[i] - (dt * C) * (U[i + 1] - U[i - 1]) / (2 * dx) + (dt * C / dx) ** 2 * (U[i + 1] - 2 * U[i] + U[i - 1]) / 2
    boundary_conditions(UN, NX)
    U[:] = UN[:]


# analytical solution
for i in range(0, NX + 1):
    Z[i] = math.sin(K * X[i] - K * C * NT * dt)


# numerical solution estimate
G = ((1 + CFL ** 2 * (np.cos(K * dx) - 1)) ** 2 + (-CFL * np.sin(K * dx)) ** 2) ** 0.5
FE = - K * dx * CFL
F = np.arctan(-CFL * np.sin(K * dx) / (1 + CFL ** 2 * (np.cos(K * dx) - 1)))
for i in range(0, NX + 1):
    g[i] = (G ** NT) * np.sin(K * X[i] + F * NT)


# saving results
output = open('output_2.txt', 'w')
output_data(output, X, UN, NX)
output.close()


# drawing graph
plt.rcParams['font.family'] = 'Courier New'

plt.rcParams['font.size'] = '24'

plt.xlabel('x, m', fontweight='bold')

plt.ylabel('u, m/s', fontweight='bold')

plt.plot(X[1:NX], UN[1:NX], 'black', linewidth = 4, label = 'numerical solution')  # numerical solution

plt.plot(X[1:NX], Z[1:NX], 'skyblue', linewidth = 4, label = 'analytical solution')  # analytical solution

plt.plot(X[1:NX], U0[1:NX], 'green', linestyle = '-.', linewidth = 3, label = 'initial signal')  # initial signal

plt.plot(X[1:NX], g[1:NX], 'red', linestyle = ':', linewidth = 4, label = 'evaluation')  # numerical
# solution estimate

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17), shadow=False, edgecolor='black', ncol=2, fontsize='small')

plt.xlim([0, 1])

plt.ylim([-1, 1])

plt.grid()

plt.show()



