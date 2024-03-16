import math as m
import matplotlib.pyplot as plt
from tabulate import tabulate

global mi, eps, d
mi = 1000
eps = 1/10**6
d = 0.1

def mf(x,cof):
    y = 0
    for i in range(len(cof)):
        y += cof[i]*x**i
    return y

def df(x,cof):
    y = 0
    t = cof.copy()
    n = len(cof)
    for i in range(len(cof)):
        t[i] *= i
        y += t[i]*x**(i-1)
    t = t[1:n]
    return y,t

def ddf(x,cof):
    y = 0
    n = len(cof)
    t = cof.copy()
    for i in range(len(cof)):
        t[i] *= i*(i-1)
        y += t[i]*x**(i-2)
    col = t[2:n]
    return y,col

def regularfalsi(f,a,b):
    if f(a)*f(b) >= 0:
        return -1
    c = a
    for i in range(mi):
        cb = c
        c = (a*f(b) - b*f(a))/(f(b) - f(a))

        if (f(c) * f(a)) < 0:
            b = c
        elif (f(c) * f(b)) < 0:
            a = c
        print("Iteration: {} | Current Solution: {}".format(i, f(c)))
        if abs(cb -c) < eps:
            return c

def newtonraphson(f,df,p):
    for i in range(mi):
        tb = p
        if abs(df(p)) > eps:
            p -= (f(p)/df(p))
            print("Iteration: {} | Current Solution: {}".format(i, f(p)))
            if abs(tb - p) < eps or abs(f(p) - d) < eps:
                return p

def bisect(f, a, b):
    for i in range(mi):
        c = (a + b)/2
        print("Iteration: {} | Current Solution: {}".format(i, f(c)))
        if abs(a - b) < eps:
            return c
        if (f(c)*f(a)) < 0:
            b = c
        elif (f(c)*f(b) < 0):
            a = c


def bracket(f,a,b):
    for i in range(mi):
        if f(a)*f(b) < 0:
            return a, b
        elif (f(a) * f(b)) > 0:
            if abs(f(a)) > abs(f(b)):
                b += d*(b-a)
            elif abs(f(a)) < abs(f(b)):
                a -= d*(b-a)

def deflation(root,cof):
    t = cof.copy()
    n = len(cof)
    for i in range(n-2,-1,-1):
        t[i] += root*t[i+1]
    t = t[1:n]
    return t

def Laguerre(cof,p):

    for i in range(mi):
        f = mf(p,cof)
        t = p
        l1 = df(p,cof)
        l2 = ddf(p,cof)
        G = l1[0]/f
        H = G*G - l2[0]/f
        n = len(cof)
        if G >= 0:
            a = n/(G+m.sqrt((n-1)*abs((n*H - G*G))))
        else:
            a = n/(G-m.sqrt((n-1)*abs((n*H - G*G))))
        p -= a
        if abs(t - p) < eps or abs(f) < eps:
            return p

def all_roots(cof,p):
    roots = []
    a = Laguerre(cof,p)
    t1 = deflation(a,cof)
    roots.append(a)
    b = Laguerre(t1,a)
    t2 = deflation(b,t1)
    roots.append(b)
    c = Laguerre(t2,b)
    t3 = deflation(c,t2)
    roots.append(c)
    d = Laguerre(t3,c+1)
    t4 = deflation(d,t3)
    roots.append(d)

    return roots

# Bisection Method for Root Finding
def bisection(f, a, b, tolerance=1e-4, plot=False):
    """
    Bisection method for finding roots of a function f(x).

    Args:
    f: The function for which roots need to be found.
    a: The lower bound of the interval [a, b].
    b: The upper bound of the interval [a, b].
    tolerance: The desired maximum error allowed.
    plot: Whether to plot the convergence graph.

    Returns:
    root: Approximation of the root found.
    """
    iterations = 0
    root_values = []

    while True:
        c = (a + b) / 2
        error = (b - a) / 2
        root_values.append(c)

        if error < tolerance or abs(f(c)) < tolerance:
            print("Convergence reached.")
            break
        else:
            iterations += 1
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c

    if plot:
        plt.figure(figsize=(8, 6))
        plt.plot(root_values, 'o-', label='Bisection Method')
        plt.xlabel('Iterations')
        plt.ylabel('Root Approximation')
        plt.title('Convergence of Bisection Method')
        plt.legend()
        plt.grid(True)
        plt.show()

    iteration_data = []
    for i, root in enumerate(root_values):
        iteration_data.append([i, root, f(root)])

    print(tabulate(iteration_data, headers=['Iteration', 'Root Approximation', 'f(Root)'], tablefmt='pretty'))

    return c

def regula_falsi(f, a, b, tolerance=1e-4, plot=False):
    """
    False Position (Regula Falsi) method for finding roots of a function f(x).

    Args:
    f: The function for which roots need to be found.
    a: The lower bound of the interval [a, b].
    b: The upper bound of the interval [a, b].
    tolerance: The desired maximum error allowed.
    plot: Whether to plot the convergence graph.

    Returns:
    root: Approximation of the root found.
    """
    iterations = 0
    root_values = []

    while True:
        c = (a * f(b) - b * f(a)) / (f(b) - f(a))
        root_values.append(c)

        if abs(f(c)) < tolerance:
            print("Convergence reached.")
            break
        else:
            iterations += 1
            if f(a) * f(c) < 0:
                b = c
            else:
                a = c

    if plot:
        plt.figure(figsize=(8, 6))
        plt.plot(root_values, 'o-', label='False Position (Regula Falsi) Method')
        plt.xlabel('Iterations')
        plt.ylabel('Root Approximation')
        plt.title('Convergence of False Position (Regula Falsi) Method')
        plt.legend()
        plt.grid(True)
        plt.show()

    iteration_data = []
    for i, root in enumerate(root_values):
        iteration_data.append([i, root, f(root)])

    print(tabulate(iteration_data, headers=['Iteration', 'Root Approximation', 'f(Root)'], tablefmt='pretty'))

    return c

def newton_raphson(f, f_prime, x0, tolerance=1e-4, plot=False):
    """
    Newton-Raphson method for finding roots of a function f(x).

    Args:
    f: The function for which roots need to be found.
    f_prime: The derivative of the function f(x).
    x0: The initial seed value for the method.
    tolerance: The desired maximum error allowed.
    plot: Whether to plot the convergence graph.

    Returns:
    root: Approximation of the root found.
    """
    iterations = 0
    root_values = []

    while True:
        x1 = x0 - f(x0) / f_prime(x0)
        root_values.append(x1)

        if abs(x1 - x0) < tolerance or abs(f(x1)) < tolerance:
            print("Convergence reached.")
            break
        else:
            iterations += 1
            x0 = x1

    if plot:
        plt.figure(figsize=(8, 6))
        plt.plot(root_values, 'o-', label='Newton-Raphson Method')
        plt.xlabel('Iterations')
        plt.ylabel('Root Approximation')
        plt.title('Convergence of Newton-Raphson Method')
        plt.legend()
        plt.grid(True)
        plt.show()

    iteration_data = []
    for i, root in enumerate(root_values):
        iteration_data.append([i, root, f(root)])

    print(tabulate(iteration_data, headers=['Iteration', 'Root Approximation', 'f(Root)'], tablefmt='pretty'))

    return x1

# Example usage:
"""
def f(x):
    return x**3 - 2*x - 5

a = 2
b = 3
tolerance = 1e-4

root_approximation = bisection(f, a, b, tolerance)
"""

# Multivariate Newton-Raphson Method for Systems of Nonlinear Equations
def solve_system(matrix, vector):
    """
    Solve a linear system of equations Ax = b.

    Args:
    matrix: Coefficient matrix (A) as a list of lists.
    vector: Right-hand side vector (b) as a list.

    Returns:
    solution: Solution vector (x) as a list.
    """
    n = len(matrix)
    m = len(matrix[0])
    assert len(vector) == n, "Dimensions of matrix and vector do not match."

    # Gaussian elimination
    augmented_matrix = [row + [vector[i]] for i, row in enumerate(matrix)]
    for i in range(n):
        # Partial pivoting
        max_row = max(range(i, n), key=lambda j: abs(augmented_matrix[j][i]))
        augmented_matrix[i], augmented_matrix[max_row] = augmented_matrix[max_row], augmented_matrix[i]
        # Elimination
        for j in range(i + 1, n):
            factor = augmented_matrix[j][i] / augmented_matrix[i][i]
            for k in range(i, m + 1):
                augmented_matrix[j][k] -= factor * augmented_matrix[i][k]

    # Back substitution
    solution = [0] * n
    for i in range(n - 1, -1, -1):
        solution[i] = augmented_matrix[i][-1] / augmented_matrix[i][i]
        for j in range(i - 1, -1, -1):
            augmented_matrix[j][-1] -= augmented_matrix[j][i] * solution[i]

    return solution

def negate_vector(vector):
    """
    Negate all elements of a vector.

    Args:
    vector: Input vector as a list.

    Returns:
    negated_vector: Negated vector as a list.
    """
    return [-x for x in vector]

def add_vectors(vector1, vector2):
    """
    Add two vectors element-wise.

    Args:
    vector1: First input vector as a list.
    vector2: Second input vector as a list.

    Returns:
    result_vector: Resultant vector as a list.
    """
    return [x + y for x, y in zip(vector1, vector2)]

def subtract_vectors(vector1, vector2):
    """
    Subtract one vector from another element-wise.

    Args:
    vector1: First input vector as a list.
    vector2: Second input vector as a list.

    Returns:
    result_vector: Resultant vector as a list.
    """
    return [x - y for x, y in zip(vector1, vector2)]

def dot_product(vector1, vector2):
    """
    Compute the dot product of two vectors.

    Args:
    vector1: First input vector as a list.
    vector2: Second input vector as a list.

    Returns:
    dot_product_value: Dot product of the two input vectors as a scalar.
    """
    return sum(x * y for x, y in zip(vector1, vector2))

def norm(vector):
    """
    Compute the Euclidean norm of a vector.

    Args:
    vector: Input vector as a list.

    Returns:
    norm_value: Euclidean norm of the input vector as a scalar.
    """
    return sum(x ** 2 for x in vector) ** 0.5

def plot_convergence(convergence_data):
    """
    Plot the convergence path of the Newton-Raphson method.

    Args:
    convergence_data: List of approximation vectors at each iteration.
    """
    x_values = [iteration[0] for iteration in convergence_data]
    y_values = [iteration[1] for iteration in convergence_data]

    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, 'o-', label='Convergence')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Convergence of Newton-Raphson Method')
    plt.legend()
    plt.grid(True)
    plt.show()

def multivariate_newton_raphson(f, jacobian, x0, tolerance=1e-6, max_iterations=1000, plot=False):
    """
    Multivariate Newton-Raphson method for finding roots of a multivariable function f(x).

    Args:
    f: The function for which roots need to be found. This function should return a list or tuple.
    jacobian: The Jacobian matrix of the function f(x).
    x0: The initial guess for the root. Should be a list or tuple.
    tolerance: The desired maximum error allowed.
    max_iterations: Maximum number of iterations allowed.
    plot: Whether to plot the convergence graph.

    Returns:
    root: Approximation of the root found.
    """
    x = x0
    convergence_data = []

    for i in range(max_iterations):
        # Step I: Compute the update direction
        s = solve_system(jacobian(x), negate_vector(f(x)))

        # Step II: Update the approximation
        x_new = add_vectors(x, s)
        convergence_data.append(x_new)

        # Check for convergence
        if norm(subtract_vectors(x_new, x)) < tolerance:
            print("Convergence reached.")
            break

        x = x_new

    else:
        print("Maximum number of iterations reached.")

    if plot:
        print(tabulate(convergence_data, headers=['x', 'y'], tablefmt='pretty'))
        plot_convergence(convergence_data)

    return x


# truncated_number: Truncated number with or without rounding off.
def truncate_number(number, decimal_places, round_off=False):
    """
    Truncate a number with or without rounding off.

    Args:
    number: The number to truncate.
    decimal_places: Number of decimal places to keep.
    round_off: Whether to round off or not. Default is True.

    Returns:
    truncated_number: Truncated number with or without rounding off.
    """
    factor = 10 ** decimal_places
    if round_off:
        truncated_number = round(number * factor) / factor
    else:
        truncated_number = int(number * factor) / factor
    return truncated_number


# fixed-point method for finding roots of the equation x = g(x)

def fixed_point_method(g, a, b, iterations=10, tolerance=1e-6, plot=True):
    """
    Fixed-point method for finding roots of the equation x = g(x).

    Args:
    g: The fixed-point function.
    a: The left initial guess.
    b: The right initial guess.
    iterations: The maximum number of iterations.
    tolerance: The tolerance for convergence.
    plot: Whether to plot the convergence graph.

    Returns:
    root: Approximation of the root found.
    """
    
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))

    x = [a + (b - a) * i / 99 for i in range(100)]
    # Start at left
    print(f"Solving x = g(x) starting to the left, at x_0 = {a}")
    x_k = a
    axs[0].set_title(f"Solving $x = g(x)$ starting to the left, at $x_0$ = {a}")
    axs[0].plot(x, x, 'g', label='$y = x$')
    axs[0].plot(x, [g(xi) for xi in x], 'r', label='$y = g(x)$')
    axs[0].grid(True)
    for k in range(iterations):
        x_k_plus_1 = g(x_k)
        axs[0].annotate(f'${truncate_number(x_k_plus_1, 4)}$', xy=(x_k, x_k), xytext=(x_k + 0.1, x_k - 0.1), arrowprops=dict(facecolor='blue', shrink=0.05))
        axs[0].plot([x_k, x_k], [x_k, x_k_plus_1], 'b')
        axs[0].plot([x_k, x_k_plus_1], [x_k_plus_1, x_k_plus_1], 'b')
        x_k = x_k_plus_1
        print(f"x_{k} = {truncate_number(x_k_plus_1, 4)}")

    # Start at right
    print(f"Solving x = g(x) starting to the right, at x_0 = {b}")
    x_k = b
    axs[1].set_title(f"Solving $x = g(x)$ starting to the right, at $x_0$ = {b}")
    axs[1].plot(x, x, 'g', label='$y = x$')
    axs[1].plot(x, [g(xi) for xi in x], 'r', label='$y = g(x)$')
    axs[1].grid(True)
    for k in range(iterations):
        x_k_plus_1 = g(x_k)
        axs[1].annotate(f'${truncate_number(x_k_plus_1, 4)}$', xy=(x_k, x_k), xytext=(x_k + 0.1, x_k - 0.1), arrowprops=dict(facecolor='blue', shrink=0.05))
        axs[1].plot([x_k, x_k], [x_k, x_k_plus_1], 'b')
        axs[1].plot([x_k, x_k_plus_1], [x_k_plus_1, x_k_plus_1], 'b')
        x_k = x_k_plus_1
        print(f"x_{k} = {truncate_number(x_k_plus_1, 4)}")

    plt.tight_layout()
    if plot:
        plt.show()
