import matplotlib.pyplot as plt
from sympy.abc import t, x, y
import numpy as np
import sympy as sy
import math


def f(x, y):
    return -6*(x+y)

#точное решение
def exact_sol(x, y):
    return x**3 + y**3

h1 = 0.1
h2 = 0.1
l1 = 1
l2 = 1
n1 = int(l1/h1)
n2 = int(l2/h2)


phi_f = [[0] * (n1+1) for _ in range(n2+1)]

#заполняем значения функции фи
for i in range(0, n1+1):
    x = i*h1
    for j in range(0, n2+1):
        y = j*h2
        phi_f[i][j] = f(x,y)
        
for i in range(1, n1):
    x = i*h1
    phi_f[i][1] += ((x**3)/(h2**2))
    phi_f[i][n2-1] += ((x**3+1)/(h2**2))

for j in range(1, n2):
    y = j*h2
    phi_f[1][j] += ((y**3)/(h1**2))
    phi_f[n1-1][j] += ((y**3+1)/(h1**2))
     
    
beta = [[0] * (n1+1) for _ in range(n2+1)]
alpha = [[0] * (n2+1) for _ in range(n1+1)]


#находим коэффициенты бета по формуле (коэф-ты разложения фи по по собственным функциям оператора)
for i in range (0, n1+1):
    for k in range(1, n2):
        for j in range(1, n2):
            beta[i][k] += h2*math.sin(k*sy.pi*j*h2/l2)*phi_f[i][j]
        beta[i][k] *= (2/l2)**0.5


#находим коэффициенты альфа методом прогонки (коэф-ты для разностного уравнения)          
q = [0] * (n1+1) 
phi = [0] * (n1+1) 
q[0] = 0
phi[0] = 0
         

for k in range(1, n2):
        lambda_k = (-4/(h2**2))*((math.sin((k*sy.pi*h2)/(2*l2)))**2)
        for i in range(0, n1):
            q[i+1] = -(1/(q[i]+(lambda_k*(h1**2) - 2)))
            phi[i+1] = ((-(h1**2 * beta[i+1][k]))-phi[i])/(q[i]+(lambda_k*(h1**2)-2))
        for t in range(n1-1, 0, -1):
            alpha[k][t] = q[t]*alpha[k][t+1] + phi[t]


v = [[0] * (n1+1) for _ in range(n2+1)]
#заполняем граничные значения v
for i in range(0, n1+1):
    x = i*h1
    v[i][0] = x**3
    v[i][n2] = x**3+1

for j in range(0, n2+1):
    y = j*h2
    v[0][j] = y**3
    v[n1][j] = y**3+1

#вычисляем их во внутренних точках 
for i in range(1, n1):
    for j in range(1, n2):
        y = j*h2
        for k in range(1, n2):
            v[i][j] += alpha[k][i]*((2/l2)**0.5)*math.sin((k*sy.pi*y)/l2)
         

exact = [[0] * (n1+1) for _ in range(n2+1)]
for i in range (0, n1+1):
    for j in range(0, n2+1):
        exact[i][j] = exact_sol(i*h1, j*h2)
print() 

    
for i in range(0, n1+1):
    print('i =', i)
    for j in range(0, n2+1):
        print('\t', round(v[i][j], 6),'\t',round(exact[i][j], 6))

max_t = 0
for i in range(0, n1+1):
    for j in range(0, n2+1):
        res = abs(exact[i][j] - v[i][j])
        if res > max_t:
            max_t = res
            
#print ('max error = ', max_t)
 
x = np.linspace(0, l1, n1+1)
y = np.linspace(0, l2, n2+1)


X, Y = np.meshgrid(x, y)
Z2 = exact_sol(X, Y)
Z = np.array(v)
    
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection='3d')
 
ax.plot_surface(X, Y, Z, cmap='cool', alpha=0.8)
ax.scatter(X, Y, Z2, color='blue', alpha=0.8)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.set_zlabel('z', fontsize=12)
 
plt.show()
