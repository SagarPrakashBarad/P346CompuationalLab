import math 
from tabulate import tabulate

global a, c, m
a = 1103515245
c = 12345
m = 32768

x0 = 0.1

def LCGPRNG():
    global x0
    x0 = (a*x0 + c) % m
    return x0

def MidInt(f, a, b, N):
    h = (b - a) / N
    M = 0
    for i in range(1, N + 1):
        M += f(a + (2*i - 1) / 2 * h)

    M *= h
    return M

def TrapInt(f, a, b, N):
    h = (b - a) / N
    M = (f(a) + f(b)) / 2
    for i in range(1, N):
        M += f(a + (i) * h)
    M *= h
    return M

def SimpsonInt(f, a, b, N):
    h = (b - a) / N
    M = f(a) + f(b)
    for i in range(1, N):
        if i % 2 == 0:
            M += 2 * f(a + (i) * h)
        else:
            M += 4 * f(a + (i) * h)
    M *= h / 3
    return M

def MonteCarloInt(f, a, b, N):
    """
    Monte Carlo integration
    """
    est = 0
    sigma = 0

    for _ in range(N):
        x = a + (b - a) * (LCGPRNG() / m)
        est += f(x)
        sigma += f(x) ** 2

    est *= (b - a) / N
    sigma -= (est / (b - a)) ** 2
    sigma = math.sqrt(sigma / N)

    return est


def gaussian_quadrature(f, a, b, n):
    if n == 2:
        roots = [-1 / math.sqrt(3), 1 / math.sqrt(3)]
        weights = [1, 1]
    elif n == 3:
        roots = [-math.sqrt(0.6), 0, math.sqrt(0.6)]
        weights = [5 / 9, 8 / 9, 5 / 9]
    elif n == 4:
        roots = [-0.861136, -0.339981, 0.339981, 0.861136]
        weights = [0.347855, 0.652145, 0.652145, 0.347855]
    # Add more cases for higher n as needed
    else:
        raise ValueError("Gaussian quadrature not implemented for n > 4")

    x_mapped = [0.5 * (b - a) * root + 0.5 * (a + b) for root in roots]
    integral = 0.5 * (b - a) * sum(weight * f(x) for weight, x in zip(weights, x_mapped))
    return integral

def compare_integration_methods(f, a, b, num_points, method_names, convergence=False):
    methods = {
        "Midpoint": MidInt,
        "Trapezoidal": TrapInt,
        "Simpson's": SimpsonInt,
        "Monte Carlo": MonteCarloInt,
        "Gaussian Quadrature": gaussian_quadrature,
    }

    headers = ["Method", "Integral Value"]
    data = []

    for method_name in method_names:
        if method_name in methods:
            method_func = methods[method_name]
            integral = method_func(f, a, b, num_points) if method_name != "Gaussian Quadrature" else method_func(f, a, b, 3)
            data.append([method_name, integral])
        else:
            print(f"Method '{method_name}' is not available.")

    print(tabulate(data, headers=headers))
    
    if convergence:
        # Compare convergence
        convergence_headers = ["N"] + method_names
        convergence_data = []
        for N in range(20, 121, 20):
            row = [N]
            for method_name in method_names:
                if method_name in methods:
                    method_func = methods[method_name]
                    integral = method_func(f, a, b, N) if method_name != "Gaussian Quadrature" else method_func(f, a, b, 3)
                    row.append(integral)
                else:
                    row.append(None)
            convergence_data.append(row)

        print("\nComparison of Convergence:")
        print(tabulate(convergence_data, headers=convergence_headers))
