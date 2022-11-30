import sys
sys.stdin = open("q1_input.txt", "r")
sys.stdout = open("q1_out.txt", "w")
sys.stderr = open("error.txt", "w")

from LabLibrary import MatrixAlgebra as mat
dim =list(map(int,input().strip().split()))[:2]
ele = list(map(float,input().strip().split()))[:(dim[0]*dim[1])]

A = mat.Matrix(dims=(dim[0],dim[1]),fill=ele)
A.LUdecomp()
A.Inverse()
Print(A.Inverse())
