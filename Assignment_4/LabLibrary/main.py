import NonLinearRoot as NL
import MatrixAlgebra as AL
import LeastSquare as LS
import pandas as pd
import numpy as np
import math

df = pd.read_csv("assign4fit.txt", sep=" ")
x = df.iloc[:,0].values
y = df.iloc[:,1].values
l = np.log(y)
print(x,l)
coeff = LS.least_squares(x,l)
LS.plot(x,y)


"""A = AL.Matrix(dims = [2,2],fill = [1,1,1,-1])
print(A)
A.gauss_siedel(B=[4,2],init_val=[1,1])"""
