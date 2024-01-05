import math as m
global a, c, m
a = 1103515245
c = 12345
m = 32768

x0 = 0.1

def LCGPRNG():
    global x0
    x0 = (a*x0 + c) % m
    return x0

def MidInt(f,a,b,N):
    h = (b-a)/N
    M = 0
    for i in range(1,N+1):
        M += f(a+(2*i-1)/2*h)

    M *= h
    return M


def trapezInt(f,a,b,N):
    h = (b-a)/N
    M = (f(a) + f(b))/2
    for i in range(1,N):
            M += f(a+(i)*h)
    M *= h
    return M


def SimpsonInt(f,a,b,N):
    h = (b-a)/N
    M = f(a) + f(b)
    for i in range(1,N):
        if i%2 == 0:
            M += 2*f(a+(i)*h)
        else:
            M += 4*f(a+(i)*h)
    M *= h/3
    return M

def monte_carlo_int(f,a,b,N):
    """
    monte carlo interation with more uniform speed
    """
    est = 0
    sigma = 0
    X = a

    for i in range(1,N+1):
        xb = X
        yb = est
        X = a + (b-a)*(LCGPRNG()/m)
        est += ((b-a)/N)*f(X)
        sigma += f(X)**2/N
        #plt.plot([xb,X],[yb,est])

    sigma -= (est/(b-a))**2
    sigma = (sigma)**(1/2)

    #plt.show()

    return est
