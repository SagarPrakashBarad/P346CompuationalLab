import numpy as np
import math
import matplotlib.pyplot as plt

def solve_bvp_finite_difference(N, h, A, y_bc_left, y_bc_right, plot_inverse=False, plot_exact_comparison=False):
    """
    Solve a boundary value problem using the finite difference method.

    Args:
    - N: Number of intervals
    - h: Step size
    - A: Matrix A
    - y_bc_left: Boundary condition at left boundary
    - y_bc_right: Boundary condition at right boundary
    - plot_inverse: Boolean flag to indicate whether to plot the inverse of matrix A
    - plot_exact_comparison: Boolean flag to indicate whether to compare with the exact solution

    Returns:
    - y: Solution vector
    """
    y = np.zeros(N+1)
    
    # Set boundary conditions
    y[0] = y_bc_left
    y[N] = y_bc_right

    b = np.zeros(N-1)
    b[0] = -y[0] / (h * h)
    b[N-2] = -y[N] / (h * h)

    invA = np.linalg.inv(A)

    if plot_inverse:
        plt.imshow(invA)
        plt.xlabel('i', fontsize=16)
        plt.ylabel('j', fontsize=16)
        plt.yticks(range(N-1), range(1, N))
        plt.xticks(range(N-1), range(1, N))
        clb = plt.colorbar()
        clb.set_label('Matrix value')
        plt.title(r'Matrix $A^{-1}$', fontsize=20)
        plt.tight_layout()
        plt.show()

    y[1:N] = np.dot(invA, b)

    if plot_exact_comparison:
        x = np.linspace(0, 1, N+1)
        fig = plt.figure(figsize=(8, 4))
        plt.plot(x, y, 'v', label='Finite Difference')
        plt.plot(x, [math.sinh(2*xi+1) for xi in x], 'k:', label='Exact')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend(loc='best')
        plt.show()

    return y


def generate_matrix_A(N, h, v, plot=False):
    """
    Generate the matrix A for a given N and h.

    Args:
    - N: Number of intervals
    - h: Step size
    - plot: Boolean flag to indicate whether to plot the matrix A

    Returns:
    - A: The generated matrix A
    """
    A = [[0] * (N-1) for _ in range(N-1)]
    
    # Diagonal
    for i in range(N-1):
        A[i][i] = -(2 / (h * h) + v)

    # Off-diagonal
    for i in range(N-2):
        A[i+1][i] = 1 / (h * h)
        A[i][i+1] = 1 / (h * h)
    
    if plot:
        plt.figure(figsize=(4,4))
        plt.imshow(A, interpolation='none', cmap='plasma')
        plt.xlabel('i', fontsize=8)
        plt.ylabel('j', fontsize=8)
        #plt.yticks(range(1, N), range(1, N))
        #plt.xticks(range(1, N), range(1, N))
        clb = plt.colorbar()
        clb.set_label('Matrix value')
        plt.title('Matrix A', fontsize=20)
        plt.tight_layout()
        plt.show()
    
    return A


def discretize_axis(N, a, b, plot=False):
    """
    Discretize the axis.

    Args:
    - N: Number of intervals
    - a: Start of the axis
    - b: End of the axis
    - plot: Boolean flag to indicate whether to plot the discretized axis

    Returns:
    - x: Discretized axis
    """
    L = b - a
    h = L / N
    x = [a + i * h for i in range(N+1)]
    
    if plot:
        # Plot the discretized axis
        fig = plt.figure(figsize=(6,4))
        plt.plot(x, [0]*(N+1), 'o:', color='red')
        plt.plot(x[0], 0, 'o:', color='green')
        plt.plot(x[-1], 0, 'o:', color='green')
        plt.xlim((a -1, b + 1))
        plt.xlabel('x', fontsize=8)
        plt.title('Illustration of discrete time points for N={}, L={}'.format(N, L), fontsize=14)
        plt.show()
    
    return x, h

# Test the function