import sys
from LabLibrary import Radioactivity as rad
import matplotlib.pyplot as plt
sys.stdin = open("input.txt", "r")
sys.stdout = open("output.txt", "w")
sys.stderr = open("error.txt", "w")


nl, nr, t = rad.BoxProblem(5000,0.5)
plt.scatter(x=t, y=nl, label="$N_l(t)$")
plt.scatter(x=t, y=nr, label="$N_r(t)$")
plt.xlabel("Time t (in secs)", size=18)
plt.ylabel("# of particles left.", size=18)
plt.legend(prop={"size": 14})
plt.savefig("q2.png")
plt.show()
