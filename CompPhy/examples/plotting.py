# importing relvant libraries
import math
import matplotlib.pyplot as plt
plt.style.use(['science','ieee'])
plt.rcParams.update({'text.usetex': False})


def f_1(x): 
    return math.exp(-x) - x

def g_1(x): 
    return math.exp(-x)

# f and g
f = f_1
g = g_1

# parameters
a = 0
b = 1
iterations = 10
x = [a + (b - a) * i / 99 for i in range(100)]


fig, axs = plt.subplots(1, 2, figsize=((8, 4)))

# f(x)
axs[0].plot(x, [f(xi) for xi in x], label='$f(x) = e^{-x} - x$')
axs[0].plot([a, b], [0, 0], 'k--', label='$y=0$')
axs[0].set_title('$y = f(x)$ and $y=0$')
axs[0].grid(True)
axs[0].legend()

# g(x)
axs[1].plot(x, [g(xi) for xi in x], label='$g(x) = e^{-x}$')
axs[1].plot(x, x, 'r--', label='$y=x$')
axs[1].set_title('$y = g(x)$ and $y=x$')
axs[1].grid(True)
axs[1].legend()

plt.tight_layout()
plt.show()