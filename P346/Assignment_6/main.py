import sys
import LabLibrary as mylib
import Eigenvalueandvectors as eigen
import matplotlib.pyplot as plt
import numpy as np

sys.stdin = open("input.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")


# question 1
omega0 = 1.2
m, c, k = 1, 0.1, 1
oscillator = mylib.HarmonicOscillator(m, c, k)

def applied_force(t):
    return omega0 ** 2

f = mylib.rk_derivative_factory(oscillator, applied_force)
x_osc = oscillator.solve(f, [2, -1])

x, dx = x_osc.T
t = np.linspace(0, 10 * oscillator.period(), 1000)

plt.scatter(x=t, y=x, label="$y(t)$")
plt.scatter(x=t, y=dx, label="$y'(t)$")

plt.xlabel("t", size=18)
plt.ylabel("y", size=18)
plt.legend(prop={"size": 14})
plt.savefig('forced harmonic oscillator.png')
plt.clf()

# question 2
guess = list(map(float,input().strip().split()))[:2]
param = list(map(int,input().strip().split()))[:3]
L = float(input())

z, y, x = mylib.boundary_value_problem(guess,param,L)

plt.scatter(x, y, color = 'r', label="$T(x)$")
p = 0
for i in range (0, len(y)):
    if abs(y[i] - param[1]) < 0.1:
        p = i
        print('At {} m the temperature of rod gets T = 100 C'.format(x[p]))
plt.axvline(x = x[p])
plt.axhline(y = 100)
plt.savefig('heat conduction.png')
plt.clf()


# question 3
hp = list(map(int,input().strip().split()))[:4]
mylib.One_D_heat(hp[0], hp[1], hp[2], hp[3])

# question 4
dim =list(map(int,input().strip().split()))[:2]
ele = list(map(float,input().strip().split()))[:(dim[0]*dim[1])]

A = eigen.Eigenvalue(dims=(dim[0],dim[1]),fill = ele)
print(A)
k = A.eigenvalue(x0 = [1,1,1])
print('Founded in {}th iteration'.format(k))
print('Eigenvalue: {}'.format(A.ev1))
A.eigenvectors()
print('Eigenvectors: {}'.format(A.ev))
