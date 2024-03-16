import math 
import matplotlib.pyplot as plt


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

def plot_f_function(f, a, b, f_label):
    """
    Plot the function f(x) within the interval [a, b].

    Args:
    f: The function f(x) to be plotted.
    a: The lower bound of the interval [a, b].
    b: The upper bound of the interval [a, b].
    f_label: Label for the function f(x).
    """
    # x values
    x = [a + (b - a) * i / 99 for i in range(100)]

    plt.figure(figsize=(6, 4))
    plt.plot(x, [f(xi) for xi in x], label='$' + f_label + '$')
    plt.plot([a, b], [0, 0], 'k--', label='$y=0$')
    plt.title('$y = ' + f_label + '$ and $y=0$')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()

def plot_functions(f, g, a, b, f_label, g_label, compare_with_yx=False):
    """
    Plot two functions f(x) and g(x) within the interval [a, b].

    Args:
    f: The first function f(x).
    g: The second function g(x).
    a: The lower bound of the interval [a, b].
    b: The upper bound of the interval [a, b].
    f_label: Label for the first function f(x).
    g_label: Label for the second function g(x).
    """
    # x values
    x = [a + (b - a) * i / 99 for i in range(100)]

    fig, axs = plt.subplots(1, 2, figsize=(8, 4))

    axs[0].plot(x, [f(xi) for xi in x], label='$' + f_label + '$')
    axs[0].plot([a, b], [0, 0], 'k--', label='$y=0$')
    axs[0].set_title('$y = ' + f_label + '$ and $y=0$')
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(x, [g(xi) for xi in x], label='$' + g_label + '$')
    axs[1].plot([a, b], [0, 0], 'k--', label='$y=0$')
    if compare_with_yx:
        axs[1].plot(x, x, 'r--', label='$y=x$')
        axs[1].set_title('$y = ' + g_label + '$ and $y=x$')
    else:
        axs[0].plot([a, b], [0, 0], 'k--', label='$y=0$')
    axs[1].grid(True)
    axs[1].legend()

    plt.tight_layout()
    plt.show()
