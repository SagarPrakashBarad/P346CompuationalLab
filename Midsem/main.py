from LabLibrary import NonLinearRoot as NL
from LabLibrary import MatrixAlgebra as AL
from LabLibrary import LeastSquare as LS
from LabLibrary import PRNG as PR
import pandas as pd
import numpy as np
import sys

sys.stdin = open("msem_gs.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")

#Question 1
N = int(input())
val = PR.area(N)
print("ellipse area: {}| error in area: {} %| No of random point chosen: {}".format(val[0], val[1], N))# The area of ellipse

#Question 2
h = 6.626*10**(-34)
k = 1.381*10**(-23)
c = 3*10**(8)
x = NL.newtonraphson(6)
b = (h*c)/(k*x)
print("Wien Constant: {0:.4f} x 10^-3".format(b*10**(3)))

print('\n')
#Question 3
dim =list(map(int,input().strip().split()))[:2]
ele = list(map(float,input().strip().split()))[:(dim[0]*dim[1])]
A = AL.Matrix(dims=[dim[0],dim[1]],fill=ele)
b = list(map(float,input().strip().split()))[:(dim[0])]
A.gauss_siedel(b,init_val=[1,1,1,1,1,1])

#Question 4
print('\n')
df = pd.read_csv("msem_fit.txt", sep=" ")
x = df.iloc[:,0].values
y = df.iloc[:,1].values
lx = np.log(x)
ly = np.log(y)
coeff = LS.least_squares(lx,ly)
print("log (a): {} | b: {} | r^2: {}".format(coeff[0], coeff[1], coeff[2]))
LS.plot_2(lx,ly)
loeff = LS.least_squares(x,ly)
print("log (a): {} | -b: {} | r^2: {}".format(loeff[0], -1*loeff[1], loeff[2]))
LS.plot_1(x,ly)
print("The power law model fits better for the data.")
