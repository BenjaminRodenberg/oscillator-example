from __future__ import division

import argparse
import numpy as np
from numpy.linalg import eig
import pandas as pd
import pathlib

from brot.enums import Cases, TimeSteppingSchemes
from brot.output import add_metainfo
from brot.timesteppers import GeneralizedAlpha
import brot.oscillator as oscillator

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

M = oscillator.M
M_inv = oscillator.M_inv
K = oscillator.K

u0 = np.array([oscillator.MassLeft.u0, oscillator.MassRight.u0])
v0 = np.array([oscillator.MassLeft.v0, oscillator.MassRight.v0])
a0 = -K.dot(M_inv.dot(u0))

analytical_1 = oscillator.MassLeft.u_analytical
analytical_2 = oscillator.MassRight.u_analytical

dts = [0.04, 0.02, 0.01, 0.005, 0.0025]
T = 1

# Generalized Alpha Parameters
if args.time_stepping == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    time_stepper = GeneralizedAlpha(stiffness=K, mass=M, alpha_f=0.4, alpha_m=0.2)
elif args.time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    time_stepper = GeneralizedAlpha(stiffness=K, mass=M, alpha_f=0.0, alpha_m=0.0)
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

    positions_1 = [u[0]]
    positions_2 = [u[1]]
    times = [t]

    while t < T:

        # do generalized alpha step
        u_new, v_new, a_new = time_stepper.do_step(u, v, a, 0, dt)

        u = u_new
        v = v_new
        a = a_new
        t += dt

        e = u[0]**2 + u[1]**2 + (u[1]-u[0])**2 + v[0]**2 + v[1]**2

        positions_1.append(u[0])
        positions_2.append(u[1])
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

filepath = this_file.parent / f"{Cases.MONOLITHIC.value}_{args.time_stepping}.csv"
df.to_csv(filepath)

add_metainfo(this_file, filepath, args.time_stepping)