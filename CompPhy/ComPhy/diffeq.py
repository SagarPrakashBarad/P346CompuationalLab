import matplotlib.pyplot as plt

def forward_euler(f, y0, total_time, dt_values, plot_label='', plot=True):
    """
    Forward (explicit) Euler's method for solving first-order ODE.

    Args:
    f: The function representing the first-order ODE dy/dt = f(t, y).
    y0: Initial value of y at t=0.
    total_time: Total time for which to compute the solution.
    dt_values: List of time step sizes.
    plot_label: Label for the plot in LaTeX format.
    plot: Boolean flag indicating whether to plot the solutions or not.

    Returns:
    solutions: A dictionary containing the computed solutions for each time step size.
    """
    solutions = {}

    for dt in dt_values:
        num_steps = int(total_time / dt) + 1
        t = [i * dt for i in range(num_steps)]
        y = [y0]
        for i in range(1, num_steps):
            y.append(y[-1] + f(t[i-1], y[-1]) * dt)
        solutions[dt] = (t, y)

        if plot:
            plt.plot(t, y, label=f'DT = {dt}')

    if plot:
        plt.xlabel('Time (t)')
        plt.ylabel('y')
        plt.title('Forward Euler Method')
        plt.legend()
        plt.grid(True)
        plt.show()

    return solutions

def backward_euler(f, y0, total_time, dt_values, plot_label='', plot=True):
    """
    Backward (implicit) Euler's method for solving first-order ODE.

    Args:
    f: The function representing the first-order ODE dy/dt = f(t, y).
    y0: Initial value of y at t=0.
    total_time: Total time for which to compute the solution.
    dt_values: List of time step sizes.
    plot_label: Label for the plot in LaTeX format.
    plot: Boolean flag indicating whether to plot the solutions or not.

    Returns:
    solutions: A dictionary containing the computed solutions for each time step size.
    """
    solutions = {}

    for dt in dt_values:
        num_steps = int(total_time / dt) + 1
        t = [i * dt for i in range(num_steps)]
        y = [y0]
        for i in range(1, num_steps):
            y.append(y[-1] / (1 - f(t[i], 1) * dt))
        solutions[dt] = (t, y)

        if plot:
            plt.plot(t, y, label=f'DT = {dt}')

    if plot:
        plt.xlabel('Time (t)')
        plt.ylabel('y')
        plt.title('Backward Euler Method')
        plt.legend()
        plt.grid(True)
        plt.show()

    return solutions

def predictor_corrector(f, y0, total_time, dt_values, plot_label='', plot=True):
    """
    Predictor-Corrector method for solving first-order ODE.

    Args:
    f: The function representing the first-order ODE dy/dt = f(t, y).
    y0: Initial value of y at t=0.
    total_time: Total time for which to compute the solution.
    dt_values: List of time step sizes.
    plot_label: Label for the plot in LaTeX format.
    plot: Boolean flag indicating whether to plot the solutions or not.

    Returns:
    solutions: A dictionary containing the computed solutions for each time step size.
    """
    solutions = {}

    for dt in dt_values:
        num_steps = int(total_time / dt) + 1
        t = [i * dt for i in range(num_steps)]
        y = [y0]
        for i in range(1, num_steps):
            y_pred = y[-1] + f(t[i-1], y[-1]) * dt
            y.append(y[-1] + 0.5 * (f(t[i-1], y[-1]) + f(t[i], y_pred)) * dt)
        solutions[dt] = (t, y)

        if plot:
            plt.plot(t, y, label=f'DT = {dt}')

    if plot:
        plt.xlabel('Time (t)')
        plt.ylabel('y')
        plt.title('Predictor-Corrector Method')
        plt.legend()
        plt.grid(True)
        plt.show()

    return solutions

def runge_kutta_4(f, y0, total_time, dt_values, plot_label=None, plot=True):
    solutions = {}
    plt.figure(figsize=(6, 4))
    for dt in dt_values:
        num_steps = int(total_time / dt)
        t_values = [i * dt for i in range(num_steps + 1)]
        y_values = [y0]

        for i in range(num_steps):
            t = t_values[i]
            y = y_values[i]

            k1 = f(t, y)
            k2 = f(t + 0.5 * dt, [y[j] + 0.5 * dt * k1[j] for j in range(len(y))])
            k3 = f(t + 0.5 * dt, [y[j] + 0.5 * dt * k2[j] for j in range(len(y))])
            k4 = f(t + dt, [y[j] + dt * k3[j] for j in range(len(y))])

            y_next = [y[j] + (dt / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(len(y))]
            y_values.append(y_next)

        solutions[dt] = (t_values, y_values)

        if plot:
            plt.plot(t_values, [y[0] for y in y_values], label=f"dt = {dt}")

    if plot:
        plt.xlabel("Time")
        plt.ylabel("y")
        if plot_label:
            plt.title(plot_label)
        else:
            plt.title("Runge-Kutta 4th Order Method")
        plt.legend()
        plt.grid(True)
        plt.show()

    return solutions


def coupled_runge_kutta_4(f, y0, total_time, dt_values, plot_label=None, plot=True):
    solutions = {}
    for dt in dt_values:
        num_steps = int(total_time / dt)
        t_values = [i * dt for i in range(num_steps + 1)]
        y_values = [y0]

        for i in range(num_steps):
            t = t_values[i]
            y = y_values[i]

            k1 = f(t, y)
            k2 = f(t + 0.5 * dt, [y[j] + 0.5 * dt * k1[j] for j in range(len(y))])
            k3 = f(t + 0.5 * dt, [y[j] + 0.5 * dt * k2[j] for j in range(len(y))])
            k4 = f(t + dt, [y[j] + dt * k3[j] for j in range(len(y))])

            y_next = [y[j] + (dt / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(len(y))]
            y_values.append(y_next)

        solutions[dt] = (t_values, y_values)

        if plot:
            plt.figure(figsize=(10, 6))
            for i in range(len(y0)):
                plt.plot(t_values, [y[i] for y in y_values], label=f"y{i+1}")
            plt.xlabel("Time")
            plt.ylabel("y")
            plt.title(plot_label if plot_label else "Coupled Runge-Kutta 4th Order Method")
            plt.legend()
            plt.grid(True)
            plt.show()

    return solutions


def plot_projections(solutions, plot_label=None):
    plt.figure(figsize=(15, 5))
    for dt, (t_values, y_values) in solutions.items():
        plt.subplot(1, 3, 1)
        plt.plot([y[0] for y in y_values], [y[1] for y in y_values], label=f"dt = {dt}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.title("xy Projection")

        plt.subplot(1, 3, 2)
        plt.plot([y[1] for y in y_values], [y[2] for y in y_values], label=f"dt = {dt}")
        plt.xlabel("y")
        plt.ylabel("z")
        plt.grid(True)
        plt.title("yz Projection")

        plt.subplot(1, 3, 3)
        plt.plot([y[2] for y in y_values], [y[0] for y in y_values], label=f"dt = {dt}")
        plt.xlabel("z")
        plt.ylabel("x")
        plt.grid(True)
        plt.title("zx Projection")
    
    plt.suptitle(plot_label if plot_label else "Projections of Lorenz Attractor")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_3d(solutions, plot_label=None):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for dt, (t_values, y_values) in solutions.items():
        ax.plot([y[0] for y in y_values], [y[1] for y in y_values], [y[2] for y in y_values], label=f"dt = {dt}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(plot_label if plot_label else "3D Plot of Lorenz Attractor")
    ax.legend()
    plt.grid(True)
    plt.show()
    
    
def semi_implicit_euler(f, y0, total_time, dt):
    num_steps = int(total_time / dt)
    t_values = [i * dt for i in range(num_steps + 1)]
    y_values = [y0]
    y = y0[:]
    for i in range(num_steps):
        y[0] += dt * y[1]
        y[1] += dt * (f(t_values[i], y)[1])
        y[2] += dt * y[3]
        y[3] += dt * (f(t_values[i], y)[3])
        y_values.append(y[:])
    return t_values, y_values


def verlet(f, y0, total_time, dt):
    num_steps = int(total_time / dt)
    t_values = [i * dt for i in range(num_steps + 1)]
    y_values = [y0]
    y = y0[:]
    for i in range(num_steps):
        y[0] += dt * y[1] + 0.5 * dt**2 * f(t_values[i], y)[1]
        y[2] += dt * y[3] + 0.5 * dt**2 * f(t_values[i], y)[3]
        y[1] += 0.5 * dt * (f(t_values[i], y)[1] + f(t_values[i+1], y)[1])
        y[3] += 0.5 * dt * (f(t_values[i], y)[3] + f(t_values[i+1], y)[3])
        y_values.append(y[:])
    return t_values, y_values


def leapfrog(f, y0, total_time, dt):
    num_steps = int(total_time / dt)
    t_values = [i * dt for i in range(num_steps + 1)]
    y_values = [y0]
    y = y0[:]
    y[0] += 0.5 * dt * y[1]
    y[2] += 0.5 * dt * y[3]
    for i in range(num_steps):
        y[0] += dt * y[1]
        y[2] += dt * y[3]
        y[1] += dt * f(t_values[i], y)[1]
        y[3] += dt * f(t_values[i], y)[3]
        y_values.append(y[:])
    return t_values, y_values


def shooting_method_bvp(f, a, b, alpha, beta, N, boundary_type, plot_label=None, plot=False):
    """
    Solve a boundary value problem using the shooting method.

    Args:
    - f: The function representing the system of first-order ordinary differential equations.
    - a, b: The endpoints of the interval [a, b].
    - alpha, beta: The boundary conditions at x = a and x = b, respectively.
    - N: The number of steps for the shooting method.
    - boundary_type: The type of boundary conditions (Dirichlet, Neumann, or mixed).
    - plot_label: The label for the plot with LaTeX formatting.
    - plot: Boolean indicating whether to plot the solution.

    Returns:
    - x: The array of x values.
    - y: The array of solution values corresponding to y.
    """

    def rk4(f, a, b, y0, N):
        h = (b - a) / N
        x = [a + i * h for i in range(N+1)]
        y = [[0] * len(y0) for _ in range(N+1)]
        y[0] = y0

        for i in range(N):
            k1 = [h * val for val in f(x[i], *y[i])]
            k2 = [h * val for val in f(x[i] + h/2, *[y[i][j] + k1[j]/2 for j in range(len(y0))])]
            k3 = [h * val for val in f(x[i] + h/2, *[y[i][j] + k2[j]/2 for j in range(len(y0))])]
            k4 = [h * val for val in f(x[i] + h, *[y[i][j] + k3[j] for j in range(len(y0))])]
            y[i+1] = [y[i][j] + (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6 for j in range(len(y0))]

        return x, y

    def shoot_dirichlet(s):
        x, y = rk4(f, a, b, [alpha, s], N)
        return y[-1][0] - beta

    def shoot_neumann(s):
        x, y = rk4(f, a, b, [alpha, s], N)
        return y[-1][1]

    def shoot_mixed(s):
        x, y = rk4(f, a, b, [alpha, s], N)
        return y[-1][0] - beta[0], y[-1][1] - beta[1]

    if boundary_type == 'Dirichlet':
        # Newton-Raphson method to find the root
        initial_slope = 0  # Initial guess for the slope
        for _ in range(10):  # Limit iterations for robustness
            x, y = rk4(f, a, b, [alpha, initial_slope], N)
            residue = y[-1][0] - beta
            if abs(residue) < 1e-6:
                break
            slope_derivative = (shoot_dirichlet(initial_slope + 1e-6) - residue) / 1e-6
            initial_slope -= residue / slope_derivative

    elif boundary_type == 'Neumann':
        # Newton-Raphson method to find the root
        initial_slope = 0  # Initial guess for the slope
        for _ in range(10):  # Limit iterations for robustness
            x, y = rk4(f, a, b, [alpha, initial_slope], N)
            residue = y[-1][1] - beta
            if abs(residue) < 1e-6:
                break
            slope_derivative = (shoot_neumann(initial_slope + 1e-6) - residue) / 1e-6
            initial_slope -= residue / slope_derivative

    elif boundary_type == 'Mixed':
        # Newton-Raphson method to find the root
        initial_slope = 0  # Initial guess for the slope
        for _ in range(10):  # Limit iterations for robustness
            x, y = rk4(f, a, b, [alpha, initial_slope], N)
            residue1 = y[-1][0] - beta[0]
            residue2 = y[-1][1] - beta[1]
            if abs(residue1) < 1e-6 and abs(residue2) < 1e-6:
                break
            slope_derivative1 = (shoot_mixed(initial_slope + 1e-6)[0] - residue1) / 1e-6
            slope_derivative2 = (shoot_mixed(initial_slope + 1e-6)[1] - residue2) / 1e-6
            determinant = slope_derivative1 * slope_derivative2
            if determinant == 0:
                break
            slope_correction1 = (residue1 * slope_derivative2 - residue2 * slope_derivative1) / determinant
            slope_correction2 = (-residue1 * slope_derivative2 + residue2 * slope_derivative1) / determinant
            initial_slope -= slope_correction1

    else:
        raise ValueError("Invalid boundary_type. Supported types are 'Dirichlet', 'Neumann', and 'Mixed'.")

    if plot:
        plt.figure(figsize=(6, 4))
        plt.plot(x, [val[0] for val in y], label='$T(x)$')
        plt.plot(3.78, 100, 'o:', color='red', label='$T\'(3.78) = 100 ^o C$')
        plt.xlabel('x')
        plt.ylabel('$T in (^o C)$')
        plt.title(plot_label)
        plt.legend()
        plt.grid(True)
        plt.show()

    return x, [val[0] for val in y]  