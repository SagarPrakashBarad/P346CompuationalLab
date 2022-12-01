import math
global a, c, m, pi, e
a = 1103515245
c = 12345
m = 32768
pi = 3.141592
e = 2.71828
x0 = 0.1

def LCGPRNG():
    global x0
    x0 = (a*x0 + c) % m
    return x0/m

def get_rand_number(min_value, max_value,i):
    """
    This functions gets a random number from a uniform distribution between
    the two input values [min_value, max_value] inclusively
    Args:
    - min_value (float)
    - max_value (float)
    Return:
    - Random number between this range (float)
    """
    range = max_value - min_value
    choice = LCGPRNG()
    res = min_value + range*choice
    if i == 0:
        return math.floor(res)
    elif i == 1:
        return res
