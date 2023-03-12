from __future__ import division

import argparse
import numpy as np
from numpy.linalg import eig
import pandas as pd
import pathlib

from enums import Cases, TimeSteppingSchemes
from output import add_metainfo

this_file = pathlib.Path(__file__)

# mass_1 = mass_2 = 1  # hard-coded for simplification of equations below
k_1, k_2, k_12 = 4*np.pi**2, 4*np.pi**2, 16*np.pi**2

# system:
# m ddu + k u = 0
#
# formulated as first order system
# dv = - k/m u
# du = v

ode_system = np.array([
    [0, 0, 1, 0],              # du_1
    [0, 0, 0, 1],              # du_2
    [-k_1 - k_12, k_12, 0, 0], # dv_1
    [k_12, -k_2 - k_12, 0, 0]  # dv_2
])

K = np.array([[k_1 + k_12, -k_12], [-k_12, k_2 + k_12]])

eigenvalues, eigenvectors = eig(K)  # should be K/M?
f = np.sqrt(eigenvalues)
A, B = eigenvectors

u0_1 = 1
u0_2 = 0
v0_1 = 0
v0_2 = 0

c = np.linalg.solve(eigenvectors, [u0_1, u0_2])

u0, v0 = np.array([u0_1, u0_2]), np.array([v0_1, v0_2])

analytical_1 = lambda t: c[0]*A[0] * np.cos(f[0] * t) + c[1]*A[1] * np.cos(f[1] * t)
analytical_2 = lambda t: c[0]*B[0] * np.cos(f[0] * t) + c[1]*B[1] * np.cos(f[1] * t)

dts = [0.04, 0.02, 0.01, 0.005, 0.0025]
T = 1

print(f"time stepping scheme being used: RK4")
print()
print("dt, error1, error2")

errors1 = []
errors2 = []

for dt in dts:
    u = u0
    v = v0
    t = 0

    positions_1 = [u[0]]
    positions_2 = [u[1]]
    velocities_1 = [v[0]]
    velocities_2 = [v[1]]
    times = [t]

    e = u[0]**2 + u[1]**2 + (u[1]-u[0])**2 + v[0]**2 + v[1]**2
    energies = [e]

    while t < T:

        # do rk4 step
        x = np.array([u[0], u[1], v[0], v[1]])

        s = 4*[None]  # stages
        s[0] = ode_system.dot(x)
        s[1] = ode_system.dot(x+s[0]*dt/2)
        s[2] = ode_system.dot(x+s[1]*dt/2)
        s[3] = ode_system.dot(x+s[2]*dt)

        x_new = x + dt/6 * (s[0] + 2*s[1] + 2*s[2] + s[3])

        u = (x_new[0], x_new[1])
        v = (x_new[2], x_new[3])
        t += dt

        e = u[0]**2 + u[1]**2 + (u[1]-u[0])**2 + v[0]**2 + v[1]**2

        positions_1.append(u[0])
        positions_2.append(u[1])
        velocities_1.append(v[0])
        velocities_2.append(v[1])
        times.append(t)
        energies.append(e)

    error1 = np.max(abs(analytical_1(np.array(times))-np.array(positions_1)))
    error2 = np.max(abs(analytical_2(np.array(times))-np.array(positions_2)))
    errors1.append(error1)
    errors2.append(error2)

    # print errors
    # analytic solution is only valid for this setup!
    print(f"{dt},{error1},{error2}")

df = pd.DataFrame(index=dts)
df.index.name = "dt"
df["error1"] = errors1
df["error2"] = errors2

time_stepping_scheme = TimeSteppingSchemes.RUNGE_KUTTA_4.value
filepath = this_file.parent / f"{Cases.MONOLITHIC.value}_{time_stepping_scheme}.csv"
df.to_csv(filepath)

add_metainfo(this_file, filepath, time_stepping_scheme)