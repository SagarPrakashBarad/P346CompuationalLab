import sys
import LabLibrary
import matplotlib.pyplot as plt

sys.stdin = open("input.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")

N = int(input()) # number of terms
r = float(input()) # common ratio
d = float(input()) # common difference

dim = list(map(int,input().strip().split()))[:2] # dimension of matrix A and B
Aelement = list(map(float,input().strip().split()))[:(dim[0]*dim[1])] # A's elements
Belement = list(map(float,input().strip().split()))[:(dim[0]*dim[1])] # B's elements
di = list(map(int,input().strip().split()))[:2] # dimension of matrix A and B
Celement = list(map(float,input().strip().split()))[:(di[0]*di[1])] # C's elements
Delement = list(map(float,input().strip().split()))[:(di[0]*di[1])] # D's elements

A = LabLibrary.Matrix(dims=(dim[0],dim[1]),fill=Aelement)
B = LabLibrary.Matrix(dims=(dim[0],dim[1]),fill=Belement)
C = LabLibrary.Matrix(dims=(di[0],di[1]),fill=Celement)
D = LabLibrary.Matrix(dims=(di[0],di[1]),fill=Delement)

z1 = list(map(int,input().strip().split()))[:2] # real and imaginary components of complex number
z2 = list(map(int,input().strip().split()))[:2]

# problem 1
print('{}, {}'.format(LabLibrary.factorial(N), LabLibrary.oddsum(N)))

# problem 2
x = LabLibrary.progressionsum(N,r,d)
print('{}, {}, {}'.format(x[0], x[1], x[1]))

# problem 3
n = N*N
y = LabLibrary.summnationseries(n)
x = [x for x in range(1, n+1)]

"""summnationseries vs n plot"""
plt.plot(x, y, 'o', color='blue')
plt.show()

# problem 4
print(A)
print(B)
print(C)
print(D)
print(LabLibrary.matmul(A, B))
print(LabLibrary.dotprod(D, C))
print(LabLibrary.matmul(B,C))

# problem 5

A = LabLibrary.mycomplex(x=z1[0],y=z1[1])
B = LabLibrary.mycomplex(x=z2[0],y=z2[1])
print('{}, {}'.format(A, A.mod))
print('{}, {}'.format(B, B.mod))
C = LabLibrary.mul(A, B)
D = LabLibrary.add(A, B)
print('{}, {}'.format(C, C.mod))
print('{}, {}'.format(D, D.mod))

"""
Expected output:
1. (factorial(N), odd(N))
3628800, 99
2. (APsum(N), GPsum(N), HPsum(N))
63.0, 1.99609375, 1.99609375
3.
A
--------------------- output ---------------------
|   2.000,   -3.000,    1.400|
|   2.500,    1.000,   -2.000|
|  -0.800,    0.000,    3.100|
--------------------------------------------------
B
--------------------- output ---------------------
|   0.000,   -1.000,    1.000|
|   1.500,    0.500,    2.000|
|   3.000,    0.000,   -2.000|
--------------------------------------------------
C
--------------------- output ---------------------
|  -2.000|
|   0.500|
|   1.500|
--------------------------------------------------
D
--------------------- output ---------------------
|   1.000|
|   0.000|
|  -1.000|
--------------------------------------------------
AB
--------------------- output ---------------------
|  -0.300,   -4.500,    9.300|
|  -3.500,   -2.000,    0.800|
|  -6.800,    8.500,   -7.000|
--------------------------------------------------
D.C
-2.0
BC
--------------------- output ---------------------
|   1.000,    0.250,   -9.000|
--------------------------------------------------
5.
z1,|z1|
3 + -2i, 3.605551275463989
z2,|z2|
1 + 2i, 2.23606797749979
z1*z2,|z1*z2|
7 + 4i, 8.06225774829855
z1+z2,|z1+z2|
4 + 0i, 4.0
"""
