from __future__ import division

import argparse
import numpy as np
from numpy.linalg import eig
import pandas as pd
import pathlib

from enums import Cases, TimeSteppingSchemes
from output import add_metainfo

this_file = pathlib.Path(__file__)

parser = argparse.ArgumentParser()
parser.add_argument(
    "-ts",
    "--time-stepping",
    help="Time stepping scheme being used.",
    type=str,
    default=TimeSteppingSchemes.NEWMARK_BETA.value,
)
args = parser.parse_args()

m_1, m_2 = 1, 1
k_1, k_2, k_12 = 4*np.pi**2, 4*np.pi**2, 16*np.pi**2

# system:
# m ddu + k u = 0

M = np.array([[m_1, 0], [0, m_2]])
M_inv = np.array([[1 / m_1, 0], [0, 1 / m_2]])
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

a0 = -K.dot(M_inv.dot(u0))

# Generalized Alpha Parameters
if args.time_stepping == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    alpha_f = 0.4
    alpha_m = 0.2
elif args.time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    alpha_f = 0.0
    alpha_m = 0.0
elif args.time_stepping == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    raise Exception(f"Please use monolithic_rk4.py for using --time-stepping=\"{args.time_stepping}\"")
else:
    raise Exception(f"Invalid time stepping scheme {args.time_stepping}. Please use one of {[ts.value for ts in TimeSteppingSchemes]}")

print(f"time stepping scheme being used: {args.time_stepping}")
print()
print("dt, error1, error2")

errors1 = []
errors2 = []

for dt in dts:
    u = u0
    v = v0
    a = a0
    t = 0

    gamma = 0.5 - alpha_m + alpha_f
    beta = 0.25 * (gamma + 0.5)

    m = 3 * [None]
    m[0] = (1 - alpha_m) / (beta * dt ** 2)
    m[1] = (1 - alpha_m) / (beta * dt)
    m[2] = (1 - alpha_m - 2 * beta) / (2 * beta)
    k_bar = K * (1 - alpha_f) + m[0] * M

    positions_1 = [u[0]]
    positions_2 = [u[1]]
    velocities_1 = [v[0]]
    times = [t]

    e = u[0]**2 + u[1]**2 + (u[1]-u[0])**2 + v[0]**2 + v[1]**2
    energies = [e]

    while t < T:

        # do generalized alpha step
        u_new = np.linalg.solve(
            k_bar, (-alpha_f * K.dot(u) + M.dot((m[0] * u + m[1] * v + m[2] * a)))
        )
        a_new = (
            1.0 / (beta * dt ** 2) * (u_new - u - dt * v) - (1 - 2 * beta) / (2 * beta) * a
        )
        v_new = v + dt * ((1 - gamma) * a + gamma * a_new)

        u = u_new
        v = v_new
        a = a_new
        t += dt

        e = u[0]**2 + u[1]**2 + (u[1]-u[0])**2 + v[0]**2 + v[1]**2

        positions_1.append(u[0])
        positions_2.append(u[1])
        velocities_1.append(v[0])
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

filepath = this_file.parent / f"{Cases.MONOLITHIC.value}_{args.time_stepping}.csv"
df.to_csv(filepath)

add_metainfo(this_file, filepath, args.time_stepping)