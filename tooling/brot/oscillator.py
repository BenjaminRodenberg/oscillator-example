import numpy as np

class SpringLeft:
    k = 4*np.pi**2


class SpringMiddle:
    k = 16*(np.pi**2)


class SpringRight:
    k = 4*np.pi**2


class MassLeft:
    # mass
    m = 1

    # initial conditions
    u0 = 1.0
    v0 = 0.0

    # analytical solution
    u_analytical = lambda t: .5 * (np.cos(2 * np.pi * t) + np.cos(6 * np.pi * t))
    v_analytical = lambda t: .5 * (-2 * np.pi * np.sin(2 * np.pi * t) - 6 * np.pi * np.sin(6 * np.pi * t))


class MassRight:
    # mass
    m = 1

    # initial conditions
    u0 = 0.0
    v0 = 0.0

    # analytical solution
    u_analytical = lambda t: .5 * (np.cos(2 * np.pi * t) - np.cos(6 * np.pi * t))
    v_analytical = lambda t: .5 * (-2 * np.pi * np.sin(2 * np.pi * t) + 6 * np.pi * np.sin(6 * np.pi * t))



# system:
# M ddu + K u = f

# Mass matrix
M = np.array([[MassLeft.m, 0          ],
              [0,          MassRight.m]])

M_inv = np.array([[1 / MassLeft.m, 0              ],
                  [0,              1 / MassRight.m]])

# Stiffness matrix
K = np.array([[SpringLeft.k + SpringMiddle.k, -SpringMiddle.k               ],
              [-SpringMiddle.k,                SpringLeft.k + SpringMiddle.k]])