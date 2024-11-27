from enum import Enum


class InterpolationSchemes(Enum):
    LAGRANGE = "Lagrange"  # Lagrange interpolation, degree equals number of points
    HERMITE = "Hermite"  # Hermite interpolation, third degree
    BSPLINE = "BSpline"


def do_linear_interpolation(time, ty_start, ty_end):
    '''
    Perform linear interpolation to obtain y value at time from given tuples (t,y) at start and end of the interval
    '''
    t_start, y_start = ty_start
    t_end, y_end = ty_end
    dt = t_end - t_start
    eps = 10**-6 * dt
    assert(t_start - eps <= time <= t_end + eps),(t_start, t_end, time)
    t = time - t_start  # time relative to beginning
    return (t/dt) * y_end + (dt-t)/dt * y_start


def do_lagrange_interpolation(time, ts, ys):
    '''
    Perform Lagrange interpolation to obtain y value at time from given tuples ts and ys
    '''
    assert(len(ts) == len(set(ts)))  # make sure that all time stamps are unique.
    assert(len(ts) == len(ys))
    n = len(ts)

    dt = ts[-1] - ts[0]
    eps = 10**-6 * dt
    # assert(min(ts) - eps <= time <= max(ts) + eps)

    return_val = 0
    for i in range(n):
        lagrange_basis_fun_i = 1
        for j in range(n):
            if i != j:
                lagrange_basis_fun_i *= (time - ts[j]) / (ts[i]-ts[j])
        return_val += ys[i] * lagrange_basis_fun_i

    return return_val


def make_cubic_spline(p0, m0, p1, m1, dt):
    """
    see https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Unit_interval_(0,_1)
    """
    #return lambda t : (2*t**3 - 3*t**2 + 1) * p0 + (t**3 - 2*t**2 + t) * dt * m0 + (-2*t**3+3*t**2) * p1 + (t**3 - t**2) * dt * m1
    return lambda t : (((2*p0 - 2*p1 + dt*m0 + dt*m1) * t + (- 3*p0 - 2*dt*m0 - dt*m1 + 3*p1)) * t + m0*dt)*t + p0


def do_cubic_interpolation(time, ty_start, ty_end):
    '''
    Perform cubic interpolation to obtain y value at time from given tuples (t,y,dydt) at start and end of the interval
    '''
    t0, p0, m0 = ty_start
    t1, p1, m1 = ty_end
    dt = t1 - t0
    interpolant = make_cubic_spline(p0, m0, p1, m1, dt)
    return interpolant((time-t0)/dt)