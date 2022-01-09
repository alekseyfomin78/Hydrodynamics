'''
Решение уравнения Бюргерса по явной противопоточной схеме
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
        U[i] = np.array(1.5)
    for i in range(NX // 5, 2 * NX // 5 + 1):
        U[i] = np.array(1.0)
    for i in range(2 * NX // 5, 3 * NX // 5 + 1):
        U[i] = np.array(1.5)
    for i in range(3 * NX // 5, 4 * NX // 5 + 1):
        U[i] = np.array(0.5)
    for i in range(4 * NX // 5, NX + 1):
        U[i] = np.array(1.5)
'''


# boundary conditions
def boundary_conditions(U, NX):
    U[0] = U[NX - 1]
    U[NX] = U[1]


# saving results
def output_data(output, X, UN, NX):
    for i in range(0, NX + 1):
        output.write(str(X[i]) + ', ' + str(UN[i]) + '\n')


with open('input_3.txt', 'r') as data:
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

U = np.array([0.0 for i in range(NX + 1)])  # divergent form

UN = np.array([0.0 for i in range(NX + 1)])  # divergent form

U2 = np.array([0.0 for i in range(NX + 1)])  # characteristic form

UN2 = np.array([0.0 for i in range(NX + 1)])  # characteristic form

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

# solution of the Burgers equation (divergent form)
for j in range(1, NT + 1):
    for i in range(0, NX + 1):
        # явная противопоточная схема (дивергентная форма)
        if U[i] > 0:
            UN[i] = U[i] - dt / dx * (0.5 * U[i] ** 2 - 0.5 * U[i - 1] ** 2)
        else:
            UN[i] = U[i] - dt / dx * (0.5 * U[i + 1] ** 2 - 0.5 * U[i] ** 2)
    boundary_conditions(UN, NX)
    U[:] = UN[:]

# solution of the Burgers equation (characteristic form)
for j in range(1, NT + 1):
    for i in range(0, NX + 1):
        # явная противопоточная схема (характеристическая форма)
        if U2[i] > 0:
            UN2[i] = U2[i] - U2[i] * dt / dx * (U2[i] - U2[i - 1])
        else:
            UN2[i] = U2[i] - U2[i] * dt / dx * (U2[i + 1] - U2[i])
    boundary_conditions(UN2, NX)
    U2[:] = UN2[:]


# saving results
output = open('output_3.txt', 'w')
output_data(output, X, UN, NX)
output.close()

# drawing graph
plt.rcParams['font.family'] = 'Courier New'

plt.rcParams['font.size'] = '24'

plt.xlabel('x, m', fontweight='bold')

plt.ylabel('u, m/s', fontweight='bold')

plt.plot(X[0:NX], UN[0:NX], 'black', linewidth=3, label='divergent form')  # divergent form

plt.plot(X[0:NX], UN2[0:NX], 'red', linewidth=3, label='characteristic form')  # characteristic form

plt.plot(X[0:NX], U0[0:NX], 'green', linestyle='-.', linewidth=2, label='initial signal')  # initial signal

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17), shadow=False, edgecolor='black', ncol=2, fontsize='small')

plt.xlim([0, 1])

plt.ylim([-1, 1])

plt.grid()

plt.show()

