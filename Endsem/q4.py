import sys
sys.stdin = open("q4_input.txt", "r")
sys.stdout = open("q4_out.txt", "w")
sys.stderr = open("error.txt", "w")
from LabLibrary import Integration as integrate
import math
L = 1
g = 9.8
Om = math.pi/4
a = math.sin(Om/2)
ll = 0
ul = math.pi/2
#ul = 2*math.asin(a*math.sin(math.pi/2))
N = 10

def Integrad(x):
    return 1/math.sqrt(1 - a**2*(math.sin(x))**2)

res = integrate.SimpsonInt(Integrad,ll,ul,N)
res = 4*math.sqrt(L/g)*res

print('T = {}'.format(res))
