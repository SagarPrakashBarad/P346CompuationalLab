import matplotlib.pyplot as plt
global a, c, m
a = 1103515245
c = 12345
m = 32768
x0 = 1

def LCGPRNG():
    global x0
    x0 = (a*x0 + c) % m
    return x0

def area(N):
    inellipse = 0

    for i in range(int(N)):
        x1 = (LCGPRNG()/(m/2))
        y1 = (LCGPRNG()/m)

        dis = (x1*x1/4 + y1*y1)

        if dis <= 1:
            inellipse += 1
    area = 8*inellipse/N # area in the first quadrant times 8
    Rarea = 2*1*3.14
    return area, abs(Rarea -area)/Rarea*100

def RandomWalk(N):
    x = 0
    y = 0
    wlkd = 0

    for i in range(1,N+1):
        xb = x
        yb = y
        x += (2*(LCGPRNG()/m) -1)
        y += (2*(LCGPRNG()/m) -1)
        plt.plot([xb,x],[yb,y])
        wlkd += (x - xb)**2 + (y- yb)**2

    plt.show()

    RmsDist = (wlkd**(1/2))
    displacement = ((x)**2 + (y)**2)**(0.5)


    return x, y, RmsDist, displacement
