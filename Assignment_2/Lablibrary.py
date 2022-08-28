import matplotlib.pyplot as plt
def LCGPRNG(seed):
    global a, c, m
    a = 1103515245
    c = 12345
    m = 32768
    seed = (a*seed + c) % m
    return seed

def volapprox(seed1, seed2, seed3, N):
    x1, y1, z1 = seed1, seed2, seed3
    insphere = 0

    for i in range(int(N)):
        x1 = LCGPRNG(x1)/m
        y1 = LCGPRNG(y1)/m
        z1 = LCGPRNG(z1)/m

        radius = (x1*x1 + y1*y1 + z1*z1)**(0.5)

        if radius <= 1:
            insphere += 1
    vol = insphere/N

    return vol

def RandomWalk(seed1 ,seed2 ,N):
    x = y = 0
    x_value = []
    y_value = []
    x_value.append(x)
    y_value.append(y)
    x1 = seed1
    y1 = seed2
    wlkd = 0
    N = int(N)

    for i in range(1,N+1):
        xl = x
        yl = y
        x1 = LCGPRNG(x1)
        y1 = LCGPRNG(y1)
        x += (2*(x1/m) -1)
        y += (2*(y1/m) -1)
        x_value.append(x)
        y_value.append(y)
        wlkd += (x - xl)**2 + (y - yl)**2

    RmsDist = (wlkd**(1/2)/N)
    displacement = (x_value[N] - x_value[0])**2 + (y_value[N] - y_value[0])**2 


    return x_value, y_value, RmsDist, displacement

def plot(x_value, y_value, filename):
    plt.style.use('ggplot')
    fig, ax = plt.subplots(1, 1)

    ax.set_xlabel("X values")
    ax.set_ylabel("Y values")
    ax.set_title("2D RandomWalk")
    ax.tick_params(labelsize=12)

    ax.plot(x_value, y_value, 'r')
    plt.savefig(filename+".pdf")
    plt.show()
