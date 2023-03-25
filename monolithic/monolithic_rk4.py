from __future__ import division

import numpy as np
import pandas as pd
import pathlib

from brot.enums import Cases, TimeSteppingSchemes
from brot.output import add_metainfo
from brot.timesteppers import RungeKutta4
import brot.oscillator as oscillator

this_file = pathlib.Path(__file__)

# system:
# m ddu + k u = 0
#
# formulated as first order system
# dv = - k/m u
# du = v

ode_system = np.array([
    [0, 0, oscillator.MassLeft.m, 0                     ],              # du_1
    [0, 0, 0,                     oscillator.MassRight.m],              # du_2
    [-oscillator.SpringLeft.k - oscillator.SpringMiddle.k, oscillator.SpringMiddle.k, 0, 0], # dv_1
    [oscillator.SpringMiddle.k, -oscillator.SpringRight.k - oscillator.SpringMiddle.k, 0, 0]  # dv_2
])

u0 = np.array([oscillator.MassLeft.u0, oscillator.MassRight.u0])
v0 = np.array([oscillator.MassLeft.v0, oscillator.MassRight.v0])

analytical_1 = oscillator.MassLeft.u_analytical
analytical_2 = oscillator.MassRight.u_analytical

dts = [0.04, 0.02, 0.01, 0.005, 0.0025]
T = 1

time_stepper = RungeKutta4(ode_system=ode_system)

print(f"time stepping scheme being used: {TimeSteppingSchemes.RUNGE_KUTTA_4.value}")
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

    while t < T:

        # do rk4 step
        f = [np.array([0, 0]) for _ in range(4)]  # no external forces for monolithic system
        u_new, v_new, a_new = time_stepper.do_step(u, v, None, f, dt)
        u = u_new
        v = v_new
        t += dt

        positions_1.append(u[0])
        positions_2.append(u[1])
        velocities_1.append(v[0])
        velocities_2.append(v[1])
        times.append(t)

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