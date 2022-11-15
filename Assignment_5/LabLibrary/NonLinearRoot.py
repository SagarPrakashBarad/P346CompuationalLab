import matplotlib.pyplot as plt
import math as m

global mi, eps, d
mi = 1000
eps = 1/10**6
d = 0.1

def mf(x,cof):
    y = 0
    for i in range(len(cof)):
        y += cof[i]*x**i
    return y

def df(x,cof):
    y = 0
    t = cof.copy()
    n = len(cof)
    for i in range(len(cof)):
        t[i] *= i
        y += t[i]*x**(i-1)
    t = t[1:n]
    return y,t

def ddf(x,cof):
    y = 0
    n = len(cof)
    t = cof.copy()
    for i in range(len(cof)):
        t[i] *= i*(i-1)
        y += t[i]*x**(i-2)
    col = t[2:n]
    return y,col

def regularfalsi(f,a,b):
    if f(a)*f(b) >= 0:
        return -1
    c = a
    for i in range(mi):
        cb = c
        c = (a*f(b) - b*f(a))/(f(b) - f(a))

        if (f(c) * f(a)) < 0:
            b = c
        elif (f(c) * f(b)) < 0:
            a = c
        print("Iteration: {} | Current Solution: {}".format(i, f(c)))
        if abs(cb -c) < eps:
            return c

def newtonraphson(f,df,p):
    for i in range(mi):
        tb = p
        if abs(df(p)) > eps:
            p -= (f(p)/df(p))
            print("Iteration: {} | Current Solution: {}".format(i, f(p)))
            if abs(tb - p) < eps or abs(f(p) - d) < eps:
                return p

def bisect(f, a, b):
    for i in range(mi):
        c = (a + b)/2
        print("Iteration: {} | Current Solution: {}".format(i, f(c)))
        if abs(a - b) < eps:
            return c
        if (f(c)*f(a)) < 0:
            b = c
        elif (f(c)*f(b) < 0):
            a = c


def bracket(f,a,b):
    for i in range(mi):
        if f(a)*f(b) < 0:
            return a, b
        elif (f(a) * f(b)) > 0:
            if abs(f(a)) > abs(f(b)):
                b += d*(b-a)
            elif abs(f(a)) < abs(f(b)):
                a -= d*(b-a)

def deflation(root,cof):
    t = cof.copy()
    n = len(cof)
    for i in range(n-2,-1,-1):
        t[i] += root*t[i+1]
    t = t[1:n]
    return t

def Laguerre(cof,p):

    for i in range(mi):
        f = mf(p,cof)
        t = p
        l1 = df(p,cof)
        l2 = ddf(p,cof)
        G = l1[0]/f
        H = G*G - l2[0]/f
        n = len(cof)
        if G >= 0:
            a = n/(G+m.sqrt((n-1)*abs((n*H - G*G))))
        else:
            a = n/(G-m.sqrt((n-1)*abs((n*H - G*G))))
        p -= a
        if abs(t - p) < eps or abs(f) < eps:
            return p

def all_roots(cof,p):
    roots = []
    a = Laguerre(cof,p)
    t1 = deflation(a,cof)
    roots.append(a)
    b = Laguerre(t1,a)
    t2 = deflation(b,t1)
    roots.append(b)
    c = Laguerre(t2,b)
    t3 = deflation(c,t2)
    roots.append(c)
    d = Laguerre(t3,c+1)
    t4 = deflation(d,t3)
    roots.append(d)

    return roots
