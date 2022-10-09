import sys
import LabLibrary
import os

sys.stdin = open("input.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")

if os.path.getsize('input.txt') == 0:
    print("File is empty! \n INPUT format: \n dim 1 dim 2 (dimension of matrix) \n element of matrix in form \n of string seperated\n by spaces 'a b c d e ...'")




dim =list(map(int,input().strip().split()))[:2] # dimensions of matrix for Question 1 and Question 2

"""Answer 1"""
e1 = list(map(float,input().strip().split()))[:(dim[0]*dim[1])]
b1 = list(map(float,input().strip().split()))[:dim[1]]
A = LabLibrary.Matrix(dims=(dim[0],dim[1]),fill=e1)

print('Answer 1: LU decomposition solution, gaussjordan solution.')
x = LabLibrary.solveLU(A,b1)
print(x) # output of LU decomposition
X, TMat = LabLibrary.gaussjordan(A,b1)
print(X) # output of gaussjordan elimination
print(TMat) # transformation matrix (RREF)
print('\n')


"""Answer 2"""
e2 = list(map(float,input().strip().split()))[:(dim[0]*dim[1])]
b2 = list(map(float,input().strip().split()))[:dim[1]]
B = LabLibrary.Matrix(dims=(dim[0],dim[1]),fill=e2)

print('Answer 2: cholesky decomposition solution, gaussSiedel solution')
x = LabLibrary.solvecholemat(B,b2)
print(x) # output of cholesky decomposition
X,c = LabLibrary.GaussSiedel(B,b2)
print(X) # output of gauss GaussSiedel
print(c) # no of interations taken by the algorithm.
print('\n')


di =list(map(int,input().strip().split()))[:2] # dimensions of matrix for Question 3
"""Answer 3"""
e3 = list(map(float,input().strip().split()))[:(di[0]*di[1])]
b3 = list(map(float,input().strip().split()))[:di[1]]
C = LabLibrary.Matrix(dims=(di[0],di[1]),fill=e3)

print('Answer 2: LU decomposition solution, gaussSiedel solution & jacobi method solution.')
x = LabLibrary.solveLU(A,b1)
print(x)
LabLibrary.makeDD(C,b3) # making the matrix daigonally dominant
X, c = LabLibrary.GaussSiedel(C,b3)
print(X) # output of gauss GaussSiedel
print(c) # no of interations taken by the algorithm.
xb, c = LabLibrary.jacobi(C,b3)
print(xb) # output of gauss jacobi method
print(c) # no of interations taken by the algorithm.

"""---------------------------OUTPUT-------------------------------------------

1. LU decomposition solution
   gaussjordan solution.
   RREF matrix
2. cholesky decomposition solution
   gaussSiedel solution
3. LU decomposition solution,
   gaussSiedel solution
   jacobi method solution.

--------------------------------------------------------------------------------"""
