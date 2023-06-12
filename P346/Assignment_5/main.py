import sys
from LabLibrary import Integration as I
import math as m


sys.stdin = open("input.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")

a, b = list(map(float,input().strip().split()))[:2]
interval = list(map(int,input().strip().split()))[:3]
c, d = list(map(float,input().strip().split()))[:2]
l, t = list(map(float,input().strip().split()))[:5]
# functions defined
def f1(x):
    return (1+1/x)**(1/2)

def f2(x):
    return m.sin(x)**2

def f3(x):
    return x**3

def f4(x):
    return x**2

# Answer 1
print('\n|-------------------Question 1-------------------------------|\n')
print('-------------------------------------------------------------------------')
print('# of intervals | ', '   Midpoint    | ', '    Trapezoidal    | ', 'Simpson')
print('      10        ', I.MidInt(f1,a,b,interval[0]) , I.trapezInt(f1,a,b,interval[0]), I.SimpsonInt(f1,a,b,interval[0]))
print('      20        ', I.MidInt(f1,a,b,interval[1]) , I.trapezInt(f1,a,b,interval[1]), I.SimpsonInt(f1,a,b,interval[1]))
print('      30        ', I.MidInt(f1,a,b,interval[2]) , I.trapezInt(f1,a,b,interval[2]), I.SimpsonInt(f1,a,b,interval[2]))
print('-------------------------------------------------------------------------')

# Answer 2
print('\n|-------------------Question 2-------------------------------|\n')
val=I.monte_carlo_int(f2,c,d,48)
print(val)

print('\n|-------------------Question 3-------------------------------|\n')
nm = I.SimpsonInt(f3,l,t,30)
dm = I.SimpsonInt(f4,l,t,30)

Xcom = nm/dm
print(Xcom)
