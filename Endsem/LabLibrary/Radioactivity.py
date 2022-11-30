import math
global a, c, m, N, T, ta, tb, x0
a = 1103515245
c = 12345
m = 32768
N = 5000
T = 1000
ta = 20
tb = 30
x0 = 0.1

def LCGPRNG():
    global x0
    x0 = (a*x0 + c) % m
    return x0

def BoxProblem(N,dt):
    Na = []
    Nb = []
    w = []
    na = N
    nb = 0
    t = 0
    n0a = N
    while (t<T):
        for i in range(n0a):
            n0a = na
            if LCGPRNG()/m <= n0a/N:
                na -= 1
                nb += 1
            elif LCGPRNG()/m <= (1-n0a/N):
                nb -= 1
                na += 1

        Na.append(na)
        Nb.append(nb)
        t += dt
        w.append(t)

    return Na, Nb, w
