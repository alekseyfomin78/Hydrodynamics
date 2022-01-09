'''
Решение задачи диффузии по центральной схеме
'''

import matplotlib.pyplot as plt
import numpy as np


# initial conditions
def initial_conditions(T1, T2, NX):
    for i in range(0, NX // 2 + 1):
        T1[i] = T_left
    for i in range(NX // 2, NX):
        T2[i] = T_right


# boundary conditions
def boundary_conditions(T1, T2, T_b, NX):
    T1[0] = T_left
    T_b[NX // 2] = (l1 * T1[NX // 2 - 1] + l2 * T2[NX // 2 + 1]) / (l1 + l2)
    T1[NX // 2] = T_b[NX // 2]
    T2[NX // 2] = T_b[NX // 2]
    T2[NX - 1] = T_right


# saving results
def output_data(output, X, TN1, TN2, NX):
    for i in range(0, NX // 2 + 1):
        output.write(str(X[i]) + ',' + str(TN1[i]) + '\n')
    for i in range(NX // 2 + 1, NX):
        output.write(str(X[i]) + ',' + str(TN2[i]) + '\n')


with open('input_1.txt', 'r') as data:
    L = float(data.readline().split('#')[0])  # длина одного тела
    T_right = float(data.readline().split('#')[0])  # температура на правом конце
    T_left = float(data.readline().split('#')[0])  # температура на левом конце
    NX = int(data.readline().split('#')[0])  # число узлов по координате
    NT = int(data.readline().split('#')[0])  # число узлов по времени
    a1 = float(data.readline().split('#')[0])
    a2 = float(data.readline().split('#')[0])
    l1 = float(data.readline().split('#')[0])
    l2 = float(data.readline().split('#')[0])

# VNM = 0.5

dx = 2 * L / (NX - 1)  # coordinate step

dt = 0.2  # time step (dx ** 2 * VNM / a1)

X = np.array([0.0 for i in range(NX)])  # nodes by coordinate

T1 = np.array([0.0 for i in range(NX)])  # temperature of the first body on the layer "n"

T2 = np.array([0.0 for i in range(NX)])  # temperature of the second body on the layer "n"

TN1 = np.array([0.0 for i in range(NX)])  # temperature of the first body on the layer "n + 1"

TN2 = np.array([0.0 for i in range(NX)])  # temperature of the second body on the layer "n + 1"

T_b = np.array([0.0 for i in range(NX)])  # temperature in the contact zone of two bodies


# nodes by coordinate
X[0] = -L
for i in range(1, NX - 1):
    X[i] = X[i - 1] + dx
X[NX - 1] = L


initial_conditions(T1, T2, NX)

boundary_conditions(T1, T2, T_b, NX)


# numerical solution
for j in range(1, NT + 1):

    for i in range(1, NX // 2):
        TN1[i] = T1[i] + (a1 * dt / dx ** 2) * (T1[i - 1] - 2 * T1[i] + T1[i + 1])

    for i in range(NX // 2 + 1, NX - 1):
        TN2[i] = T2[i] + (a2 * dt / dx ** 2) * (T2[i - 1] - 2 * T2[i] + T2[i + 1])

    boundary_conditions(TN1, TN2, T_b, NX)
    T1[:] = TN1[:]
    T2[:] = TN2[:]


# saving results
output = open('output_1.txt', 'w')
output_data(output, X, TN1, TN2, NX)
output.close()


# drawing graph
plt.rcParams['font.family'] = 'Courier New'

plt.rcParams['font.size'] = '24'

plt.xlabel('Coordinate', fontweight='bold')

plt.ylabel('Temperature', fontweight='bold')

plt.plot(X[0:NX // 2 + 1], TN1[0:NX // 2 + 1], 'black', linewidth=3, label='explicit central scheme')

plt.plot(X[NX // 2:NX], TN2[NX // 2:NX], 'black', linewidth=3)

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), shadow=False, edgecolor='black', ncol=2)

plt.xlim([-1.1, 1.1])

plt.grid()

plt.show()
