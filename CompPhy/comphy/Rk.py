import numpy as np
import matplotlib.pyplot as plt

global N, alpha, Ta, max_iter, eps
N = 1000
alpha = 0.01
Ta = 20
max_iter = 20
eps = 10**(-3)
g = 10

def RungeKutta(f, y0, x):
    y = np.zeros((len(x), len(y0)))
    y[0, :] = np.array(y0)
    h = x[1] - x[0]
    for i in range(0, len(x) - 1):
        # Many slight changes below
        k1 = np.array(f(y[i, :], x[i]))
        k2 = np.array(f(y[i, :] + h * k1 / 2, x[i] + h / 2))
        k3 = np.array(f(y[i, :] + h * k2 / 2, x[i] + h / 2))
        k4 = np.array(f(y[i, :] + h * k3, x[i] + h))
        y[i + 1, :] = y[i, :] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


# A constant function in your case, but this can be any function of `t`
def applied_force(t):
    # Note, you did not provide a value for `omega0`
    return omega0**2*math.cos(omega0*t)


def rk_derivative_factory(osc, F):
    return lambda X, t: np.array([X[1], (F(t) - osc.c * X[1] - osc.k * X[0]) / osc.m])

class HarmonicOscillator:
    def __init__(self, m, c, k):
        if (type(m) in (int, float)) and (m > 0):
            self.m = m
        else:
            raise ValueError("Parameter 'm' must be a positive number")
        if (type(c) in (int, float)) and (c > 0):
            self.c = c
        else:
            raise ValueError("Parameter 'c' must be a positive number")
        if (type(k) in (int, float)) and (k > 0):
            self.k = k
        else:
            raise ValueError("Parameter 'k' must be a positive number")

        self.T = 2 * np.pi * (self.m / self.k)**(0.5)

    def period(self):
        return self.T

    def solve(self, func, y0):
        t = np.linspace(0, 10 * self.period(), 1000)
        return RungeKutta(func, y0, t)

def One_D_heat(L,T,Nx,Nt):
    V0 = []
    V1 = []
    x = []
    hx = L/Nx
    ht = T/Nt
    alpha = ht/(hx*hx)
    if alpha >= 0.5:
        print('Stability compromised!!!')

    for i in range(Nx+1):
        V1.append(0)
        if i*hx == 1:
            V0.append(300)
        else:
            V0.append(0)
        x.append(i*hx)

    plt.plot(x, V0, label = '0')

    for j in range(N):
        for i in range(Nx+1):
            if i == 0:
                V1[i] = alpha*V0[i+1] + (1- 2*alpha)*V0[i]
            elif i == Nx:
                V1[i] = alpha*(V0[i-1]) + (1- 2*alpha)*V0[i]
            else:
                V1[i] = alpha*(V0[i+1] + V0[i-1]) + (1- 2*alpha)*V0[i]

        for i in range(Nx+1):
            V0[i] = V1[i]

        if j == 10 or j == 100 or j == 200 or j == 500 or j == 1000 or j == 2000 or j == 3000 or j == 4000 or j == 5000:
            plt.plot(x, V1, label = str(j)+'s')
            plt.xlabel("$x$", size=18)
            plt.ylabel("$T^{o} \ C$", size=18)
            plt.legend(prop={"size": 14})
    plt.savefig('1-dimensional heat.png')


def dyx(y, z, x):
    f = z
    return f

def dzx(y, z, x):
    f = - g
    return f

def ORK4(dyx,dzx,y0, z0, upper, p):#p determines if i want graph or value for interpolation
    yi = y0
    zi = z0
    x = 0
    yl = [y0]
    zl = [z0]
    xl = [0]
    h = (upper-0)/N
    while yi <= upper:
        k1y = h*dyx(yi, zi, x)
        k1z = h*dzx(yi, zi, x)

        k2y = h*dyx(yi + k1y/2, zi+k1z/2, x+h/2)
        k2z = h*dzx(yi + k1y/2, zi+k1z/2, x+h/2)

        k3y = h*dyx(yi+k2y/2, zi+k2z/2, x+h/2)
        k3z = h*dzx(yi+k2y/2, zi+k2z/2, x+h/2)

        k4y = h*dyx(yi+k3y, zi+k3z, x+h)
        k4z = h*dzx(yi+k3y, zi+k3z, x+h)

        yi += (k1y+2*k2y+2*k3y+k4y)/6
        zi += (k1z+2*k2z+2*k3z+k4z)/6
        yl.append(yi)
        zl.append(zi)
        x += h
        xl.append(x)

    if p == 0:
        return zl, yl, xl
    else:
        return zi, yi

def boundary_value_problem(g1,p1,L):
    z1h, y1 = ORK4(p1[0], g1[0], L, 1)
    z2h, y2 = ORK4(p1[0], g1[1], L, 1)

    iter = 0

    while abs(y1-p1[2]) >= eps and abs(y2-p1[2]) >= eps and iter <= max_iter:
        iter +=1
        z = g1[1] + ((g1[0] - g1[1])*(p1[2] - y2))/(y1 - y2)
        z3h, y3 = ORK4(p1[0], z, L, 1)
        if abs(y3 - p1[2]) < eps:
            zlistf, ylistf, xlistf = ORK4(p1[0], z, L, 0)
            return zlistf, ylistf, xlistf
        else:
            if y3 < p1[2]:
                g1[1] = z
                y2 = y3
            else:
                g1[0] = z
                y1 = y3

def solve(a,b,c,x):
    plt.plot(xlistf, ylistf, color = 'r')

    p = 0
    for i in range (0, len(ylistf)):
        if abs(ylistf[i] - 100) < 0.1:
            p = i
            break
    plt.axvline(x = xlistf[p])
    plt.axhline(y = 100)
    plt.show()
