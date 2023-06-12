import sys
sys.stdin = open("input.txt", "r")
sys.stdout = open("q6_out.txt", "w")
sys.stderr = open("error.txt", "w")

from LabLibrary import Eigenvalueandvectors as eigen

dim =list(map(int,input().strip().split()))[:2]
ele = list(map(float,input().strip().split()))[:(dim[0]*dim[1])]

A = eigen.Eigenvalue(dims=(dim[0],dim[1]),fill = ele)
print(A)
k = A.eigenvalue(x0 = [1,1,1,1])
print('Founded in {}th iteration'.format(k))
print('Eigenvalue: {}'.format(A.ev1))
A.eigenvectors()
print('Eigenvectors: {}'.format(A.ev))
