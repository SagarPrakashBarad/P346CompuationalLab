import sys
from LabLibrary import NonLinearRoot as nl
from LabLibrary import LeastSquare as ls
from LabLibrary import MatrixAlgebra as ma
import math as m
import pandas as pd

sys.stdin = open("input.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")

a, b = list(map(float,input().strip().split()))[:2]
c, d = list(map(float,input().strip().split()))[:2]
cof = list(map(float,input().strip().split()))[:5]
# functions defined
def f1(x):
    return m.log(5*x/2) - m.sin(5*x/2)

def f2(x):
    return -x - m.cos(x)

def df2(x):
    return -1 + m.sin(x)

# Answer 1
print('\n|-------------------Question 1-------------------------------|\n')
if nl.regularfalsi(f1,a,b) == -1:
    a, b = nl.bracket(f1,a,b)

print("Root of the equation using regular falsi: {}".format(nl.regularfalsi(f1,a,b)))
print('\n****************************************************************\n')
print("Root of the equation using bisection: {}".format(nl.bisect(f1,a,b)))

# Answer 2
print('\n|-------------------Question 2-------------------------------|\n')
print("Root of the equation using newton raphson: {}".format(nl.newtonraphson(f2,df2,c)))
print('\n****************************************************************\n')
if nl.regularfalsi(f2,c,d) == -1:
    c, d = nl.bracket(f2,c,d)

print("Root of the equation using regular falsi: {}".format(nl.regularfalsi(f2,c,d)))
print('\n****************************************************************\n')
print("Root of the equation using bisection: {}".format(nl.bisect(f2,c,d)))

# Answer 3
print('\n|-------------------Question 3-------------------------------|\n')
print(nl.all_roots(cof,0.5))

# Answer 4
print('\n|-------------------Question 4-------------------------------|\n')
df = pd.read_csv("assign4fit.txt", sep=" ")
x = df.iloc[:,0].values
y = df.iloc[:,1].values
# polyfit subroutine
m, j = ls.make_mat(x,y,3)
S = ma.Matrix(dims=[4,4],fill=m)
k, l = ma.gaussjordan(S,j)
print("polyfit: {} + {} x + {} x^2 + {}x^3".format(k[0],k[1],k[2],k[3]))
