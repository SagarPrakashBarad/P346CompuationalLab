import sys
sys.stdin = open("q5_input.txt", "r")
sys.stdout = open("q5_output.txt", "w")
sys.stderr = open("error.txt", "w")
from LabLibrary import Rk as rk
import math

g =10

def dyx(y, z, x):
    f = z
    return f

def dzx(y, z, x):
    f = - g
    return f

res = rk.ORK4(dyx, dzx,0,10,5,0)
max_hieght = max(res[1])

gamma = 0.02

def dyx(y, z, x):
    f = z
    return f

def dzx(y, z, x):
    f = - g - gamma*z
    return f

z,y,x = rk.ORK4(dyx, dzx, 0,10,max_hieght,0)

plt.scatter(x=y, y=z, label="$V(t) \ vs \ y(t)$")
#plt.scatter(x=t, y=dx, label="$y'(t)$")

plt.xlabel("x", size=18)
plt.ylabel("y", size=18)
plt.legend(prop={"size": 14})
plt.savefig('q5.png')
plt.show()
