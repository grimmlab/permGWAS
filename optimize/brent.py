# Brent's method

def brent_search(f, a: float, b: float, x: float = None, fx: float = None, rel_tol: float = 1.48e-08,
                 abs_tol: float = 1.48e-08, max_iter: int = 500, **kwargs) -> tuple:
    """
    Find minimum of a function using Brent's method (see Numerical Recipes 3rd Edition: The Art of Scientific Computing)
    Given a function f with minimum in interval [a,b], find local minimum.

    :param f: function to be minimized
    :param a: lower bound of interval
    :param b: upper bound of interval
    :param x: starting point (initial guess of minimum)
    :param fx: function value of f
    :param rel_tol: relative tolerance, default=1.48e-08
    :param abs_tol: absolute tolerance, default=1.48e-08
    :param max_iter: maximal number of iterations, default=500
    :param kwargs: additional arguments of f

    :return: minimum x, function value of minimum f(x) and number of iterations
    """

    golden = 0.381966011250105097
    if a > b:
        raise ValueError('Interval boundaries do not fit. a must be smaller or equal to b.')
    if x is None:
        x = a + golden * (b-a)
    if fx is None:
        fx = f(x, **kwargs)
    if not (a <= x <= b):
        raise ValueError('Starting value x needs to be within interval boundaries.')

    # initialize values
    x_sec, fx_sec = x, fx  # second best value and function value
    x_trd, fx_trd = x, fx  # third best value and function value
    d, e = 0.0, 0.0  # step size and direction of last two iterations
    i = -1

    for i in range(max_iter):
        mid = 0.5 * (a + b)
        tol1 = rel_tol * abs(x) + abs_tol
        tol2 = 2.0 * tol1

        # check stopping crit
        if abs(x - mid) <= tol2 - 0.5 * (b - a):
            break

        # compute Lagrange polynomial through (x, f(x)), (x_sec, f(x_sec)) and (x_trd, f(x_trd))
        if abs(e) > tol1:
            tmp1 = (x - x_sec) * (fx - fx_trd)
            denominator = (x - x_trd) * (fx - fx_sec)
            numerator = (x - x_trd) * denominator - (x - x_sec) * tmp1
            denominator = 2.0 * (denominator - tmp1)
            if denominator > 0.0:
                numerator = -numerator
            denominator = abs(denominator)
            tmp1 = e
            e = d

            if (abs(numerator) >= abs(0.5 * denominator * tmp1)) or (numerator <= denominator * (a-x)) or \
                    (numerator >= denominator * (b-x)):
                # golden section step
                e = b-x if x < mid else a-x
                d = golden * e
            else:
                # polynomial interpolation step
                d = numerator / denominator
                x_new = x + d
                if (x_new - a < tol2) or (b - x_new < tol2):
                    d = tol1 if x < mid else -tol1
        else:
            # golden section step
            e = b - x if x < mid else a - x
            d = golden * e

        # function must not be evaluated too close to x
        if tol1 <= abs(d):
            x_new = x + d
        elif 0.0 < d:
            x_new = x + tol1
        else:
            x_new = x - tol1
        fx_new = f(x_new, **kwargs)

        # check if x_new is better than previous x
        if fx_new <= fx:
            # decrease interval size
            if x_new >= x:
                a = x
            else:
                b = x
            # replace previous best 3 with current best 3
            x_trd, fx_trd = x_sec, fx_sec
            x_sec, fx_sec = x, fx
            x, fx = x_new, fx_new
        else:
            # decrease interval size
            if x_new < x:
                a = x_new
            else:
                b = x_new
            # check if x_new better than second or third and replace accordingly
            if fx_new <= fx_sec or x_sec == x:
                x_trd, fx_trd = x_sec, fx_sec
                x_sec, fx_sec = x_new, fx_new
            elif fx_new <= fx_trd or x_trd == x or x_trd == x_sec:
                x_trd, fx_trd = x_new, fx_new

    return x, fx, i+1
