import matplotlib.pyplot as plt
import math

global mi, eps, d
mi = 1000
eps = 1/10**6
d = 1/10**4

def makefc(x,cof,pow):
    y = 0
    for i in range(len(cof)):
        y += cof[i]*x**pow[i]
    return y

def fc(x):
    return (x-5)*math.exp(x) + 5

def dfc(x):
    return (x-4)*math.exp(x)

def ddfc(x):
    return (6*x)

def regularfalsi(a,b):
    if func(a)*func(b) >= 0:
        return -1

    c = a

    for i in range(mi):
        cb = c
        c = (a*func(b) - b*func(a))/(func(b) - func(a))

        if (func(c) * func(a)) < 0:
            b = c
        elif (func(c) * func(b)) < 0:
            a = c
        if abs(cb -c) < eps:
            print(i)
            return c

def newtonraphson(p):
    t = p
    for i in range(mi):
        tb = t
        if dfc(t) > eps:
            t -= (fc(t)/dfc(t))
            if abs(tb - t) < eps or abs(fc(t) -d) < eps:
                return t


def bisect(a, b):
    t1 = a
    t2 = b

    for i in range(mi):
        if abs(t1 - t2) < eps:
            print(i)
            return c
        c = (t1 + t2)/2
        if (func(c)*func(t1)) < 0:
            t2 = c
        elif (func(c)*func(t2) < 0):
            t1 = c


def bracket1(a,b,k):
    t1 = a
    t2 = b

    for i in range(mi):
        if func(t1)*func(t2) < 0:
            print(i)
            return regularfalsi(t1,t2)
        elif (func(t1) * func(t2)) > 0:
            if abs(func(t1)) > abs(func(t2)):
                t1 += k*(b-a)
            elif abs(func(t1)) < abs(func(t2)):
                t1 += k*(b-a)

def bracket2(a,b,k):
    t1 = a
    t2 = b

    for i in range(mi):
        if func(t1)*func(t2) < 0:
            print(i)
            return bisect(t1,t2)
        elif (func(t1) * func(t2)) > 0:
            if abs(func(t1)) > abs(func(t2)):
                t1 += k*(b-a)
            elif abs(func(t1)) < abs(func(t2)):
                t1 += k*(b-a)

def deflation(root,cof):
    dfcof = cof.copy()
    for i in range(len(root)):
        for j in range(len(root)):
            cof[j+1] += root[i]*cof[j]
            dfcof[j] = cof[j]*(len(root)-j)
        dfcof[len(cof)-1] = 0
        print(cof)
        print(dfcof)

    return cof,dfcof
