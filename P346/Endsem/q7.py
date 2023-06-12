import sys
import pandas as pd
from LabLibrary import LeastSquare as ls
from LabLibrary import MatrixAlgebra as ma
import matplotlib.pyplot as plt
sys.stdin = open("q7_input.txt", "r")
sys.stdout = open("q7_out.txt", "w")
sys.stderr = open("error.txt", "w")

df = pd.read_csv("esem4fit.txt", sep=" ")
x = df.iloc[:,0].values
y = df.iloc[:,1].values
# polyfit subroutine
m, j = ls.make_mat(x,y,4)
S = ma.Matrix(dims=[5,5],fill=m)
k, l = ma.gaussjordan(S,j)
print("polyfit: {} + {} x + {} x^2 + {}x^3 + {}x^4 ".format(k[0],k[1],k[2],k[3],k[4]))

def fit(x):
    return k[0] + k[1]*x + k[2]*x**2 + k[3]*x**3 + k[4]*x**4

(max,min) = (max(x),min(x))

yfit = []
xfit = []
xfit.append(min)
xfit.append(max)
yfit.append(fit(min))
yfit.append(fit(max))

plt.scatter(x=x, y=y, label="$data$")
plt.plot(xfit, yfit)
plt.xlabel("x", size=18)
plt.ylabel("y", size=18)
plt.legend(prop={"size": 14})
plt.savefig("q7.png")
plt.show()
