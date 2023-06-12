import sys
from LabLibrary import NonLinearRoot as nonlinear
import math
sys.stdin = open("q3_input.txt", "r")
sys.stdout = open("q3_out.txt", "w")
sys.stderr = open("error.txt", "w")

def f(x):
    return 2.5 - x*math.exp(x)

def g(x):
    return -x*math.exp(x) - math.exp(x)

res = nonlinear.newtonraphson(f,g,3)
print('The string can be streched to x = {}'.format(res))
