from typing import Tuple, List
import numpy as np
import numbers
import scipy as sp

class GeneralizedAlpha():
    alpha_f = None
    alpha_m = None
    gamma = None
    beta = None
    mass = None
    stiffness = None

    def __init__(self, stiffness, mass, alpha_f = 0.4, alpha_m = 0.2) -> None:
        self.alpha_f = alpha_f
        self.alpha_m = alpha_m

        self.gamma = 0.5 - self.alpha_m + self.alpha_f
        self.beta = 0.25 * (self.gamma + 0.5)

        self.stiffness = stiffness
        self.mass = mass

    def rhs_eval_points(self, dt) -> List[float]:
        return np.array([(1-self.alpha_f) * dt])

    def do_step(self, u, v, a, f, dt) -> Tuple[float, float, float]:
        if type(f) is list or type(f) is np.ndarray:  # if f is list, turn it into a number
            assert(len(f) == 1)
            f = f[0]

        m = 3*[None]
        m[0] = (1-self.alpha_m)/(self.beta*dt**2)
        m[1] = (1-self.alpha_m)/(self.beta*dt)
        m[2] = (1-self.alpha_m-2*self.beta)/(2*self.beta)

        k_bar = self.stiffness * (1 - self.alpha_f) + m[0] * self.mass

        # do generalized alpha step
        if (type(self.stiffness)) is np.ndarray:
            u_new = np.linalg.solve(
                k_bar,
                (f - self.alpha_f * self.stiffness.dot(u) + self.mass.dot((m[0]*u + m[1]*v + m[2]*a)))
            )
        else:
            u_new = (f - self.alpha_f * self.stiffness * u + self.mass * (m[0]*u + m[1]*v + m[2]*a))/k_bar

        a_new = 1.0 / (self.beta * dt**2) * (u_new - u - dt * v) - (1-2*self.beta) / (2*self.beta) * a
        v_new = v + dt * ((1-self.gamma)*a+self.gamma*a_new)

        return u_new, v_new, a_new

class RungeKutta4():
    a = np.array([[0,   0,   0,   0],
                  [0.5, 0,   0,   0],
                  [0,   0.5, 0,   0],
                  [0,   0,   1.0, 0]])
    b = np.array([1/6, 1/3, 1/3, 1/6])
    c = np.array([0, 0.5, 0.5, 1])

    def __init__(self, ode_system) -> None:
        self.ode_system = ode_system
        pass

    def rhs_eval_points(self, dt) -> List[float]:
        return np.array([self.c[0] * dt, self.c[1]*dt, self.c[2]*dt, self.c[3]*dt])

    def do_step(self, u, v, a, f, dt) -> Tuple[float, float, float]:
        assert(type(u) == type(v))

        n_stages = 4

        if isinstance(f, numbers.Number):  # if f is number, assume constant f
            f = n_stages*[f]

        if type(u) is np.ndarray:
            x = np.concatenate([u, v])
            rhs = [np.concatenate([np.array([0, 0]), f[i]]) for i in range(n_stages)]
        elif isinstance(u, numbers.Number):
            x = np.array([u, v])
            rhs = [np.array([0, f[i]]) for i in range(n_stages)]
        else:
            raise Exception(f"Cannot handle input type {type(u)} of u and v")


        s = n_stages * [None]
        s[0] = self.ode_system.dot(x)                       + rhs[0]
        s[1] = self.ode_system.dot(x+self.a[1,0] * s[0]*dt) + rhs[1]
        s[2] = self.ode_system.dot(x+self.a[2,1] * s[1]*dt) + rhs[2]
        s[3] = self.ode_system.dot(x+self.a[3,2] * s[2]*dt) + rhs[3]

        x_new = x

        for i in range(n_stages):
            x_new += dt * self.b[i] * s[i]

        if type(u) is np.ndarray:
            u_new = x_new[0:2]
            v_new = x_new[2:4]
        elif isinstance(u, numbers.Number):
            u_new = x_new[0]
            v_new = x_new[1]

        a_new = None

        return u_new, v_new, a_new


class RadauIIA():
    def __init__(self, ode_system) -> None:
        self.ode_system = ode_system
        pass

    def rhs_eval_points(self, dt) -> List[float]:
        return np.linspace(0, dt, 5)  # will create an interpolant from this later

    def do_step(self, u, v, a, f, dt) -> Tuple[float, float, float]:
        from brot.interpolation import do_lagrange_interpolation

        ts = self.rhs_eval_points(dt)

        if isinstance(f, numbers.Number):  # if f is number, assume constant f
            f = len(ts)*[f]

        t0 = 0

        assert(type(u) == type(v))

        if type(u) is np.ndarray:
            x0 = np.concatenate([u, v])
            f = np.array(f)
            assert(u.shape[0] == f.shape[1])
            rhs_fun = lambda t,x: np.concatenate([np.array([np.zeros_like(t), np.zeros_like(t)]), [do_lagrange_interpolation(t, ts, f[:,i]) for i in range(u.shape[0])]])
        elif isinstance(u, numbers.Number):
            x0 = np.array([u, v])
            rhs_fun = lambda t,x: np.array([np.zeros_like(t), do_lagrange_interpolation(t, ts, f)])
        else:
            raise Exception(f"Cannot handle input type {type(u)} of u and v")

        def fun(t,x):
            return self.ode_system.dot(x) + rhs_fun(t,x)

        ret = sp.integrate.solve_ivp(fun, [t0, t0+dt], x0, method="Radau", first_step=dt, max_step=dt, rtol=10e10, atol=10e10)  # use large rtol and atol to circumvent error control.

        a_new = None
        if type(u) is np.ndarray:
            u_new, v_new = ret.y[0:2,-1], ret.y[2:4,-1]
        elif isinstance(u, numbers.Number):
            u_new, v_new = ret.y[:,-1]

        return u_new, v_new, a_new