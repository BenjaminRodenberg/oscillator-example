from typing import Tuple, Callable
import numpy as np
from numpy.typing import NDArray
import scipy as sp
from scipy.integrate import OdeSolution
from enum import Enum


def eye_like(array):
    """usage similar to ones_like and zeros_like to create identity matrix"""
    if (len(array.shape) != 2):
        raise ValueError(f"Input array has shape {array.shape}, but 2D array is required")

    return np.eye(array.shape[0], array.shape[1])


class TimeSteppingSchemes(Enum):
    NEWMARK_BETA = "Newmark_beta"
    GENERALIZED_ALPHA = "generalized_alpha"
    RUNGE_KUTTA_4 = "runge_kutta_4"
    Radau_IIA = "radauIIA"  # 5th order implicit


class TimeStepper:
    """Time stepper for a mass spring system with given mass and spring stiffness
    """

    def __init__(self,
                 stiffness: NDArray,
                 mass: NDArray) -> None:
        """Initializes time stepper with parameters describing mass spring system

        Args:
            stiffness (NDArray): the stiffness of the spring connected to the mass
            mass (NDArray): the mass
        """
        raise NotImplementedError()

    def do_step(self,
                u0: NDArray,
                v0: NDArray,
                a0: NDArray,
                rhs: Callable[[float], NDArray],
                dt: float
                ) -> Tuple[NDArray, NDArray, NDArray]:
        """Perform a time step of size dt with given time stepper

        Args:
            u0 (NDArray): displacement at time t0
            v0 (NDArray): velocity at time t0
            a0 (NDArray): acceleration at time t0
            rhs (Callable[[float], NDArray]): time dependent right-hand side
            dt (float): time step size

        Returns:
            Tuple[NDArray, NDArray, NDArray]: returns computed displacement, velocity, and acceleration at time t1 = t0 + dt.
        """
        raise NotImplementedError()


class GeneralizedAlpha(TimeStepper):
    """TimeStepper implementing generalized Alpha or Newmark Beta scheme (depends on parameters alpha_f and alpha_m set in constructor)
    """
    alpha_f: float
    alpha_m: float
    gamma: float
    beta: float
    mass: NDArray
    stiffness: NDArray

    def __init__(self, stiffness: NDArray, mass: NDArray, alpha_f: float = 0.4, alpha_m: float = 0.2) -> None:
        self.alpha_f = alpha_f
        self.alpha_m = alpha_m

        self.gamma = 0.5 - self.alpha_m + self.alpha_f
        self.beta = 0.25 * (self.gamma + 0.5)

        self.stiffness = stiffness
        self.mass = mass

    def do_step(self,
                u0: NDArray,
                v0: NDArray,
                a0: NDArray,
                rhs: Callable[[float], NDArray],
                dt: float
                ) -> Tuple[NDArray, NDArray, NDArray]:
        f = rhs((1.0 - self.alpha_f) * dt)

        m = [(1 - self.alpha_m) / (self.beta * dt**2),
             (1 - self.alpha_m) / (self.beta * dt),
             (1 - self.alpha_m - 2 * self.beta) / (2 * self.beta)]

        k_bar = self.stiffness * (1 - self.alpha_f) + m[0] * self.mass

        # do generalized alpha step
        u1 = np.linalg.solve(
            k_bar,
            (f - self.alpha_f * self.stiffness.dot(u0) +
             self.mass.dot((m[0] * u0 + m[1] * v0 + m[2] * a0)))
        )
        a1 = 1.0 / (self.beta * dt**2) * (u1 - u0 - dt * v0) - \
            (1 - 2 * self.beta) / (2 * self.beta) * a0
        v1 = v0 + dt * ((1 - self.gamma) * a0 + self.gamma * a1)

        return u1, v1, a1


class RungeKutta4(TimeStepper):
    """TimeStepper implementing classic Runge Kutta scheme (4 stages)
    """
    # parameters from Butcher tableau of classic Runge Kutta scheme (RK4)
    a = np.array([[0, 0, 0, 0],
                  [0.5, 0, 0, 0],
                  [0, 0.5, 0, 0],
                  [0, 0, 1.0, 0]])
    b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6])
    c = np.array([0, 0.5, 0.5, 1])

    def __init__(self, stiffness: NDArray, mass: NDArray) -> None:
        mass_inverse = np.linalg.inv(mass)
        self.ode_system = np.block([
            [np.zeros_like(stiffness), eye_like(stiffness)],
            [-stiffness.dot(mass_inverse), np.zeros_like(stiffness)]
        ])

    def do_step(self,
                u0: NDArray,
                v0: NDArray,
                a0: NDArray,
                rhs: Callable[[float], NDArray],
                dt: float
                ) -> Tuple[NDArray, NDArray, NDArray]:
        k = 4 * [None]  # store stages in list

        x0 = np.concatenate([u0, v0])
        def f(t, x): return self.ode_system.dot(x) + \
            np.concatenate([np.zeros_like(u0), rhs(t)])

        k[0] = f(self.c[0] * dt, x0)
        k[1] = f(self.c[1] * dt, x0 + self.a[1, 0] * k[0] * dt)
        k[2] = f(self.c[2] * dt, x0 + self.a[2, 1] * k[1] * dt)
        k[3] = f(self.c[3] * dt, x0 + self.a[3, 2] * k[2] * dt)

        x1 = x0 + dt * sum(b_i * k_i for k_i, b_i in zip(k, self.b))

        return x1[:len(u0)], x1[len(u0):], f(dt, x1)[len(u0):]


class RadauIIA():
    """Perform a step with the adaptive RadauIIA time stepper of scipy.integrate
    """
    _dense_output: OdeSolution

    def __init__(self, stiffness: NDArray, mass: NDArray) -> None:
        mass_inverse = np.linalg.inv(mass)
        self.ode_system = np.block([
            [np.zeros_like(stiffness), eye_like(stiffness)],
            [-stiffness.dot(mass_inverse), np.zeros_like(stiffness)]
        ])

    def do_step(self,
                u0: NDArray,
                v0: NDArray,
                a0: NDArray,
                rhs: Callable[[float], NDArray],
                dt: float
                ) -> Tuple[NDArray, NDArray, NDArray]:
        t0, t1 = 0, dt

        x0 = np.concatenate([u0, v0])
        def f(t, x): return self.ode_system.dot(x) + \
            np.concatenate([np.zeros(len(u0)), rhs(t)])

        # use adaptive time stepping; dense_output=True allows us to sample from continuous function later
        # use large rtol and atol to circumvent error control.
        ret = sp.integrate.solve_ivp(
            f, [t0, t1], x0, method="Radau", first_step=dt, max_step=dt, rtol=10e10, atol=10e10)

        self._dense_output = ret.sol  # store dense output in class

        return ret.y[:len(u0), -1], ret.y[len(u0):, -1], f(dt, ret.y[:, -1])[len(u0):]

    @property
    def dense_output(self) -> OdeSolution:
        """Returns dense output created for the last call of do_step

        Returns:
            OdeSolution: dense output
        """
        return self._dense_output
