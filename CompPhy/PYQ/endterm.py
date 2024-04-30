import itertools
from statistics import median_high
import numpy as np
import plotly.graph_objects as go

from numpy.typing import ArrayLike
from scipy.optimize import fsolve

### ---- USED IN ENDTERM ----
## Polyfit & Condition Number
# Gassian Elimination / Gauss-Jordan
def gauss_elim(A, b):
    """
    Solves Ax = b linear systems using Gassian Elimination

    Parameters
    ----------
    A: numpy.ndarray
        Contains coefficients of the variables (Matrix, A)
    b: numpy.ndarray
        Contains the constants on RHS (Vector, b)
    
    Returns
    -------
    bool
        Whether the process converged within ``iter_lim``
    numpy.ndarray
        Obtained solution
    float
        Error in the obtained solution
    
    """
    # Prepping the Augmented Matrix
    aug_mat = np.concatenate((A, np.reshape(b, (-1, 1))), axis=1)
    # Convergence Flag
    CONV_FLAG = True
    # Position of leading nonzero, nondiagonal-element in a row / pivot
    lead = 0
    
    # aug_mat.shape[0] == No. of rows
    # aug_mat[0].shape or aug_mat.shape[1] == Number of Columns
    rowCount = aug_mat.shape[0]
    columnCount = aug_mat.shape[1]

    for r in range(rowCount):
        if lead >= columnCount:
            CONV_FLAG = False
            break
        i = r

        # Finding the pivot in a column
        while aug_mat[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    CONV_FLAG = False
                    break

        aug_mat[i], aug_mat[r] = aug_mat[r], aug_mat[i] # Swapping rows
        lv = aug_mat[r][lead]
        aug_mat[r] = [mrx / float(lv) for mrx in aug_mat[r]]
        for i in range(rowCount):
            if i != r:
                lv = aug_mat[i][lead]         
                aug_mat[i] = [iv - lv*rv for rv,iv in zip(aug_mat[r], aug_mat[i])]
        lead += 1
    
    if not CONV_FLAG:
        raise Exception(f"Solution did not converge.")
    
    # Returning convergence flag, solution and associated error
    return CONV_FLAG, aug_mat[:, -1], A @ aug_mat[:, -1] - b

def legendre(n, x):
    legendre = {
        0: 1,
        1: x,
        2: (1 / 2) * (3 * x ** 2 - 1),
        3: (1 / 2) * (5 * x ** 3 - 3 * x),
        4: (1 / 8) * (35 * x ** 4 - 30 * x ** 2 + 3),
        5: (1 / 8) * (63 * x ** 5 - 70 * x ** 3 + 15 * x),
        6: (1 / 16) * (231 * x ** 6 - 315 * x ** 4 + 105 * x ** 2 - 5)
    }

    return legendre[n]

def func_to_fit(x, coeffs:tuple, basis:callable):
    f = 0
    for i, coeff in enumerate(coeffs):
        f += coeffs[i] * basis(i, x)

    return f

def polyfit(x, y, yerr, n=4, basis="poly"):
    """
    Polynomial Least-Squares Fitting

    """
    bases = {
        "poly": poly,
        "cheby": cheby,
        "modcheby": mod_cheby,
        "legendre": legendre,
    }

    basis = bases[basis]

    N = n + 1
    A = np.empty((N, N))
    b = np.empty(N)

    for i in range(N):
        for j in range(N):
            A[i, j] = np.sum(basis(i, x) * basis(j, x) / yerr ** 2)
            b[j] = np.sum((basis(j, x) * y) / yerr ** 2)

    # Using Gauss-Elimination to calculate x
    return A, gauss_elim(A, b)

def chi_sq_by_v(x, y, yerr, coeffs):
    v = y.shape[0] - len(coeffs)
    return np.sum(((y - func_to_fit(x, coeffs, basis=legendre)) / yerr) ** 2) / v

def jacobi_iter(A, b, init_val, iter_lim=100, tol=1e-8, info=False):
    """
    Solves Ax = b linear systems using Jacobi Iteration

    Parameters
    ----------
    A: numpy.ndarray
        Contains coefficients of the variables (Matrix, A)
    b: numpy.ndarray
        Contains the constants on RHS (Vector, b)
    init_val: numpy.ndarray
        Contains an initial guess for x
    iter_lim: int
        Maximum number of iterations
        Defaults to ``100``
    tol: float
        Tolerance value
        Defaults to ``1e-8``
    info: bool
        Whether to store residue & iteration steps
        Defaults to ``False``
    
    Returns
    -------
    bool
        Whether the process converged within ``iter_lim``
    numpy.ndarray
        Obtained solution
    float
        Error in the obtained solution
        
    Optionally Returns
    ------------------
    meta: list
        List containing iteration steps & residue per step
    
    """
    CONV_FLAG = False # Convergence Flag
    var = init_val # Vector, X
    
    # To store residue & iteration steps
    meta = []
    
    for i in range(iter_lim):
        var_new = np.zeros_like(var) # stores updated values of all variables (Vector, X)

        for j in range(A.shape[0]):
            # Matrix Multiplying all elements, before A's diagonal (in a row) with all corresponding vars (in Vector, X)
            d = np.dot(A[j, :j], var[:j])
            # Matrix Multiplying all elements, after A's diagonal (in a row) with all corresponding vars (in Vector, X)
            r = np.dot(A[j, j + 1:], var[j + 1:])
            # Updating values of vars
            var_new[j] = (b[j] - d - r) / A[j, j]

        meta.append([i, np.linalg.norm(var - var_new)]) # Storing iteration step and residue

        # Checks, how close the updated values are, to the previous iteration's values and breaks the loop, if close enough (defined by "tol")
        if np.allclose(var, var_new, atol=tol, rtol=0.):
            CONV_FLAG = True
            break

        var = var_new # Storing the new solution

    # If solution is not obtained (no convergence), after iter_lim iterations | Note that, this "else" block belongs to the previous "for" statement and not any "if" statement
    else:
        CONV_FLAG = False
        
    if not CONV_FLAG:
        raise Exception(f"Solution did not converge, after the specified limit of {iter_lim} iterations.")
    
    # Returning convergence flag, metadata, solution and associated error
    if info:
        return CONV_FLAG, meta, var, A @ var - b
    
    # Returning convergence flag, solution and associated error
    return CONV_FLAG, var, A @ var - b

def inverse(A, iter_lim=100, tol=1e-8):
    """
    Finds the inverse of the input matrix using Jac ???

    Parameters
    ----------
    A: numpy.ndarray
        Input matrix
    method: str
        Method to use to calculate inverse
    iter_lim: int
        Maximum number of iterations
        Defaults to ``100``
    tol: float
        Tolerance value
        Defaults to ``1e-8``
    
    Returns
    -------
    numpy.ndarray
        Inverse of the input matrix

    """
    # if method not in METHODS:
    #     raise NotImplementedError(f"{method} not found or not implemented!")

    B = np.identity(A.shape[0])
    init_val = np.zeros(A.shape[0])
    cols = []

    meta_overall = []
    for b in B:
        _, meta, var, _ = jacobi_iter(A, b, init_val, iter_lim, tol, info=True)
        meta_overall.append(meta)
        cols.append(var)

    return meta_overall, np.c_[cols].T


## Random Numbers
def mlcg(seed, a=1664525, m=2**32):
    """
    Multiplicative / Lehmer Linear Congruential Generator
    Uses values suggested by Numerical Recipes as default

    Parameters
    ----------
    seed: float
        Seed
    a: float
        Multiplier
        Defaults to ``1664525``
    m: float
        Modulus
        Defaults to ``2**32``

    Returns
    -------
    float:
        Sequence of Pseudo-Random Numbers, having a period, that is sensitive to the choice of a and m.

    """
    while True:
        seed = (a * seed) % m
        yield seed

def mlcgList(N:int, range:tuple, seed:float=42, a:float=1664525, m:float=2**32):
    """
    Continuous
    Returns normalized list of MLCG-generated random values
    Uses values suggested by Numerical Recipes as default

    Parameters
    ----------
    N: int
        Number of random numbers to be generated
    range: tuple
        Range from where the numbers will be sampled
    seed: float
        Seed
        Defaults to ``42``
    a: float
        Multiplier
        Defaults to ``1664525``
    m: float
        Modulus
        Defaults to ``2**32``

    Returns
    -------
    numpy.ndarray
        Normalized list of MLCG-generated random values

    """    
    start, end = range

    rnList = np.array(list(itertools.islice(mlcg(seed, a, m), 0, N))) / m

    return end * rnList + start

# General function to get random samples
def sampler(N:int, range:tuple=(0, 1), sampling:str="uniform", seed:int=42, a=1664525, m=2**32, **samp_dict):
    # Uniform
    samples = mlcgList(N, range, seed=seed, a=a, m=m)

    if sampling == "gaussian":
        samples, _ = uniform_to_gauss(N, seed=seed)
    elif sampling == "exponential":
        samples = uniform_to_exp(samples, **samp_dict)
    elif sampling == "powerlaw":
        samples = uniform_to_power(samples, **samp_dict)
    else:
        samples = samples

    return samples

## Random Walk Simulator up to N-D
def random_walk(ndim:int=1, steps:int=1000, movement:str="grid", sampling:str="uniform", seed:int=42, **sampling_dict):
    """
    Simulates random walk
    Upto 3D
    
    NOTE: In case of power law and other distributions with extra parameters, also pass those paramaters (as kwargs).

    """
    origin = np.zeros((1, ndim))

    try:
        samples = sampler(N=ndim * steps, sampling=sampling, seed=seed, **sampling_dict)
    except:
        raise Exception("Sampling method not implemented.")

    # Sample points between 0 & 1 for making the choice
    moves = samples.reshape(steps, ndim)
    
    if movement == "grid": # DEFAULT PARAM
        # 2 choices - along gridlines only
        moves[moves >= 0.5] = 1.
        moves[moves < 0.5] = -1.
    elif movement == "diag":
        # 3 choices - diagonal movement allowed or stopping allowed in case of 1D
        moves[moves <= 0.34] = -1.
        moves[(moves > 0.34) & (moves < 0.69)] = 0
        moves[moves >= 0.69] = 1.

    path = np.concatenate([origin, moves]).cumsum(0)
    
    return path

def rw_plot(path:ArrayLike):
    """
    Random Walk Plotter
    
    """
    ndims = path.shape[1]
    start = path[:1]
    end = path[-1:]

    def plot_1d_rw(path:ArrayLike):
        x = np.arange(path.shape[0])
        path = path.squeeze()

        # Plotting the path
        fig = go.Figure(data=[
            go.Scatter(x=x, y=path,
                mode='lines+markers',
                marker=dict(
                    size=6,
                    opacity=0.8
                ),
                name="Path"
            )
        ])

        # Plotting starting and ending points
        fig.add_trace(go.Scatter(x=x[:1], y=start[0], name="Start", mode='markers', marker_symbol='diamond', marker=dict(size=8, color="green")))
        fig.add_trace(go.Scatter(x=x[-1:], y=end[0], name="End", mode='markers', marker_symbol='x', marker=dict(size=8, color="red")))

        return fig
    
    def plot_2d_rw(path:ArrayLike):
        # Plotting the path
        fig = go.Figure(data=[
            go.Scatter(x=path[:, 0], y=path[:, 1],
                mode='lines+markers',
                marker=dict(
                    size=1,
                    opacity=0.8
                ),
                name="Path"
            )
        ])
        # Plotting starting and ending points
        fig.add_trace(go.Scatter(x=start[:, 0], y=start[:, 1], name="Start", mode='markers', marker_symbol='diamond', marker=dict(size=8, color="green")))
        fig.add_trace(go.Scatter(x=end[:, 0], y=end[:, 1], name="End", mode='markers', marker_symbol='x', marker=dict(size=8, color="red")))

        fig.update_layout(
            width=600,
            height=600,
            scene = dict(
                xaxis = dict(title='X',),
                yaxis = dict(title='Y',),
            ),
        )

        return fig

    def plot_3d_rw(path:ArrayLike):
        # Plotting the path
        fig = go.Figure(data=[
            go.Scatter3d(x=path[:, 0], y=path[:, 1], z=path[:, 2],
                mode='lines+markers',
                marker=dict(
                    size=1,
                    opacity=0.8
                ),
                name="Path"
            )
        ])
        # Plotting starting and ending points
        fig.add_trace(go.Scatter3d(x=start[:, 0], y=start[:, 1], z=start[:, 2], name="Start", mode='markers', marker_symbol='diamond', marker=dict(size=4, color="green")))
        fig.add_trace(go.Scatter3d(x=end[:, 0], y=end[:, 1], z=end[:, 2], name="End", mode='markers', marker_symbol='x', marker=dict(size=4, color="red")))

        fig.update_layout(
            width=600,
            height=600,
            margin=dict(l=0, r=0, b=0, t=0,),
            scene = dict(
                xaxis = dict(title='X',),
                yaxis = dict(title='Y',),
                zaxis = dict(title='Z',),
            ),
        )

        return fig

    if ndims == 1:
        fig = plot_1d_rw(path)
    elif ndims == 2:
        fig = plot_2d_rw(path)
    else:
        fig = plot_3d_rw(path)
        
    return fig

def k_rand_walker(k:int=10, ndim:int=1, steps:int=100, **kwargs):
    """
    Simulates k random walks

    """
    paths = []
    plots = []
    # Simulating k random walks
    for i in range(k):
        path = random_walk(ndim=ndim, steps=steps, seed=(i + 1) * 100, **kwargs)
        paths.append(path)
        plot = rw_plot(path)
        plots.append(plot)

    # Combining all plots
    fig = go.Figure(data=[i.data[0] for i in plots])
    fig.update_layout(
            width=600,
            height=600,
            scene = dict(
                xaxis = dict(title='X',),
                yaxis = dict(title='Y',),
            ),
        )

    return paths, fig


### ---- UNUSED IN ENDTERM ----
## LU & Cholesky for Linear Systems
# Forward & Backward Substitution
def forward_sub(A, b, init_val, iter_lim=100, tol=1e-8):
    """
    Solves Ax = b linear systems using Forward Substitutions

    Parameters
    ----------
    A: numpy.ndarray
        Contains coefficients of the variables (Matrix, A)
    b: numpy.ndarray
        Contains the constants on RHS (Vector, b)
    init_val: numpy.ndarray
        Contains an initial guess for x
    iter_lim: int
        Maximum number of iterations
        Defaults to ``100``
    tol: float
        Tolerance value
        Defaults to ``1e-8``
    
    Returns
    -------
    bool
        Whether the process converged within ``iter_lim``
    numpy.ndarray
        Obtained solution
    float
        Error in the obtained solution
    
    """
    CONV_FLAG = False # Convergence Flag
    var = init_val # Vector, X
    
    for i in range(iter_lim):
        var_new = np.zeros_like(var) # stores updated values of all variables (Vector, X)

        for j in range(A.shape[0]):
            # Matrix Multiplying all elements, after A's diagonal (in a row) with all corresponding vars (in Vector, X)
            l = np.dot(A[j, :j], var[:j])
            # Updating values of vars
            var_new[j] = (b[j] - l) / A[j, j]

        # Checks, how close the updated values are, to the previous iteration's values and breaks the loop, if close enough (defined by "tol")
        if np.allclose(var, var_new, atol=tol, rtol=0.):
            CONV_FLAG = True
            break

        var = var_new # Storing the new solution
    # If solution is not obtained (no convergence), after iter_lim iterations | Note that, this "else" block belongs to the previous "for" statement and not any "if" statement
    else:
        CONV_FLAG = False
        
    if not CONV_FLAG:
        raise Exception(f"Solution did not converge, after the specified limit of {iter_lim} iterations.")

    # Returning convergence flag, solution and associated error
    return CONV_FLAG, var, A @ var - b

def back_sub(A, b, init_val, iter_lim=100, tol=1e-8):
    """
    Solves Ax = b linear systems using Backward Substitutions

    Parameters
    ----------
    A: numpy.ndarray
        Contains coefficients of the variables (Matrix, A)
    b: numpy.ndarray
        Contains the constants on RHS (Vector, b)
    init_val: numpy.ndarray
        Contains an initial guess for x
    iter_lim: int
        Maximum number of iterations
        Defaults to ``100``
    tol: float
        Tolerance value
        Defaults to ``1e-8``
    
    Returns
    -------
    bool
        Whether the process converged within ``iter_lim``
    numpy.ndarray
        Obtained solution
    float
        Error in the obtained solution
    
    """
    CONV_FLAG = False # Convergence Flag
    var = init_val # Vector, X
    
    for i in range(iter_lim):
        var_new = np.zeros_like(var) # stores updated values of all variables (Vector, X)

        for j in range(A.shape[0]):
            # Matrix Multiplying all elements, after A's diagonal (in a row) with all corresponding vars (in Vector, X)
            u = np.dot(A[j, j + 1:], var[j + 1:])
            # Updating values of vars
            var_new[j] = (b[j] - u) / A[j, j]
        
        # Checks, how close the updated values are, to the previous iteration's values and breaks the loop, if close enough (defined by "tol")
        if np.allclose(var, var_new, atol=tol, rtol=0.):
            CONV_FLAG = True
            break

        var = var_new # Storing the new solution
    # If solution is not obtained (no convergence), after iter_lim iterations | Note that, this "else" block belongs to the previous "for" statement and not any "if" statement
    else:
        CONV_FLAG = False
        
    if not CONV_FLAG:
        raise Exception(f"Solution did not converge, after the specified limit of {iter_lim} iterations.")

    # Returning convergence flag, solution and associated error
    return CONV_FLAG, var, A @ var - b

# LU Decomposition
def LU(A):
    """
    Returns LU Decomposition of A

    Parameters
    ----------
    A: numpy.ndarray
        Matrix, to be decomposed

    Returns
    -------
    L: numpy.ndarray
        Lower Triangular Factor of A
    U: numpy.ndarray
        Upper Triangular Factor of A

    """
    U = A.copy()
    L = np.identity(A.shape[0], dtype=float)

    for i in range(A.shape[0]):
        factor = U[i+1:, i] / U[i, i]
        L[i+1:, i] = factor
        U[i+1:] -= factor[:, np.newaxis] * U[i] # :, newaxis helps to turn factor into a 2D array of shape (N, 1) / column vector

    return L, U

# Cholesky Decomposition
def cholesky(A):
    """
    Returns the Cholesky Factor of A
    NOTE: DOES NOT CHECK FOR A's ELIGIBILITY TO BE CHOLESKY-DECOMPOSED

    Parameters
    ----------
    A: numpy.ndarray
        Matrix, to be decomposed

    Returns
    -------
    L: numpy.ndarray
        Cholesky Factor of A

    """
    N = A.shape[0]

    # Create zero matrix for L
    L = np.zeros_like(A)

    for i in range(N):
        for k in range(i + 1):
            tmp = sum(L[i, j] * L[k, j] for j in range(k))
            
            if (i == k): # Diagonal elements
                L[i, k] = np.sqrt(A[i, i] - tmp)
            else:
                L[i, k] = (1.0 / L[k, k] * (A[i, k] - tmp))

    return L

def LUorCholesSolver(A, b, init_val, iter_lim=100, tol=1e-8, method="LU"):
    """
    Solves Ax = b linear systems using LU or Cholesky Decomposition
    This in turn uses Forward and Backward Substitutions

    Parameters
    ----------
    A: numpy.ndarray
        Contains coefficients of the variables (Matrix, A)
    b: numpy.ndarray
        Contains the constants on RHS (Vector, b)
    init_val: numpy.ndarray
        Contains an initial guess for x
    iter_lim: int
        Maximum number of iterations
        Defaults to ``100``
    tol: float
        Tolerance value
        Defaults to ``1e-8``
    method: str
        "LU" or "cholesky"
        Defaults to "LU"

    Returns
    -------
    bool
        Whether the process converged within ``iter_lim``
    numpy.ndarray
        Obtained solution
    float
        Error in the obtained solution

    """
    if method not in ["LU", "cholesky"]:
        raise Exception("Method not implemented yet! Supported methods are 'LU' and 'choles'")

    if method == "LU":
        L, U = LU(A)
    elif method == "cholesky":
        L = cholesky(A)

    CONV_FLAG_FS, y, ERR_FS = forward_sub(L, b, init_val, iter_lim, tol)
    if CONV_FLAG_FS:
        CONV_FLAG_BS, x, ERR_BS = back_sub(U if method == "LU" else L.T, y, init_val, iter_lim, tol)
    else:
        raise Exception(f"Solution did not converge, after the specified limit of {iter_lim} iterations.")

    if not CONV_FLAG_BS:
        raise Exception(f"Solution did not converge, after the specified limit of {iter_lim} iterations.")

    return CONV_FLAG_BS, x, ERR_BS


## Polyfit
def poly(n, x):
    return x ** n

def mod_cheby(n, x):
    mod_chebyshev = {
        0: 1,
        1: 2 * x - 1,
        2: 8 * x ** 2 - 8 * x + 1,
        3: 32 * x ** 3 - 48 * x ** 2 + 18 * x - 1
    }

    return mod_chebyshev[n]

def cheby(n, x):
    chebyshev = {
        0: 1,
        1: x,
        2: 2 * x ** 2 - 1,
        3: 4 * x ** 3 - 3 * x,
        4: 8 * x ** 4 - 8 * x ** 2 + 1,
        5: 16 * x ** 5 - 20 * x ** 3 + 5 * x,
        6: 32 * x ** 6 - 48 * x ** 4 + 18 * x ** 2 - 1,
        7: 64 * x ** 7 - 112 * x ** 5 + 56 * x ** 3 - 7 * x,
        8: 128 * x ** 8 - 256 * x ** 6 + 160 * x ** 4 - 32 * x ** 2 + 1,
        9: 256 * x ** 9 - 576 * x ** 7 + 432 * x ** 5 - 120 * x ** 3 + 9 * x
    }

    return chebyshev[n]

def cond(x:float, coeffs:ArrayLike, basis:callable):
    B = 0
    for i, coeff in enumerate(coeffs):
        B += np.abs(coeffs[i]) * np.abs(basis(i, x))
    
    return B


## Random Numbers
# Gaussian RNDataset Generator
def GRD(U0, U1):
    """
    Box-Muller Transform

    """
    Z0 = np.sqrt(-2 * np.log(U0)) * np.cos(2 * np.pi * U1)
    Z1 = np.sqrt(-2 * np.log(U0)) * np.sin(2 * np.pi * U1)

    return Z0, Z1

def uniform_to_gauss(N:int, seed:int):
    U0 = mlcgList(N, (0, 1), seed=seed)
    U1 = mlcgList(N, (0, 1), seed=seed+20)

    Z0, Z1 = GRD(U0, U1)

    return Z0, Z1

# Exponential - alpha is a parameter here
def p(x, alpha:float):
    return alpha * np.exp(-alpha * x)

def uniform_to_exp(y, alpha:float=1):
    return -(1 / alpha) * np.log(1 - y)

# Power Law
def power_law(x, a:float, n:float):
    return a * x ** n

def uniform_to_power(y, lims:tuple, a:float, n:float):
    """
    Source: https://mathworld.wolfram.com/RandomNumber.html
    Here, y is uniformly sampled from [0, 1].

    """
    x0,  x1 = lims

    return ((x1 ** (n + 1) - x0 ** (n + 1)) * y + x0 ** (n + 1)) ** (1 / (n + 1))

def uniform_to_power_neg(y, lims:tuple, a:float, n:float):
    """
    Source: https://mathworld.wolfram.com/RandomNumber.html, but with p(x) = a * (x + 1) ** (-n)
    Here, y is uniformly sampled from [0, 1].
    And, n is implicitly negative.

    """
    x0,  x1 = lims

    return (((x1 + 1) ** (1 - n) - (x0 + 1) ** (1 - n)) * y + (x0 + 1) ** (1 - n)) ** (1 / (1 - n)) - 1


## Backward Euler
def backward_euler_res(y_next, f, t_curr, y_curr, t_next):
    return y_next - y_curr - (t_curr - t_next) * f (t_next, y_next)

def backward_euler(f, y_init, t_lims, N):
    M = 1 if not np.ndim(y_init) else len(y_init)
    t_start, t_end = t_lims
    dt = (t_end - t_start) / float(N)

    t = np.zeros(N + 1)
    t[0] = t_start
    y = np.zeros((N + 1, M))
    y[0, :] = y_init

    for i in range (0, N):
        t_curr = t[i]
        y_curr = y[i, :]
        t_next = t[i] + dt
        y_next = y_curr + dt * f(t_curr, y_curr)

        y_next = fsolve(backward_euler_res, y_next, args=(f, t_curr, y_curr, t_next))

        t[i + 1]   = t_next
        y[i + 1,:] = y_next[:]

    return t, y
