from __future__ import division

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from brot.enums import TimeSteppingSchemes
from brot.timeSteppersMonolithic import GeneralizedAlpha, RungeKutta4, RadauIIA
import brot.oscillator as oscillator

from io import TextIOWrapper
from numpy.typing import ArrayLike

parser = argparse.ArgumentParser()
parser.add_argument("-tss",
                    "--time-stepping-scheme",
                    help=f"Time stepping scheme being used. Please use one of {[ts.value for ts in TimeSteppingSchemes]}",
                    type=str,
                    default=TimeSteppingSchemes.NEWMARK_BETA.value)
parser.add_argument("-dt", "--time-step-size", help=f"Time step size being used", type=float, default=0.04)
args = parser.parse_args()

M = oscillator.M
K = oscillator.K

u0 = np.array([oscillator.MassLeft.u0, oscillator.MassRight.u0])
v0 = np.array([oscillator.MassLeft.v0, oscillator.MassRight.v0])
a0 = -K.dot(np.linalg.inv(M).dot(u0))

analytical_1 = oscillator.MassLeft.u_analytical
analytical_2 = oscillator.MassRight.u_analytical

T = 1

# Generalized Alpha Parameters
if args.time_stepping_scheme == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    time_stepper = GeneralizedAlpha(stiffness=K, mass=M, alpha_f=0.4, alpha_m=0.2)
elif args.time_stepping_scheme == TimeSteppingSchemes.NEWMARK_BETA.value:
    time_stepper = GeneralizedAlpha(stiffness=K, mass=M, alpha_f=0.0, alpha_m=0.0)
elif args.time_stepping_scheme == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    time_stepper = RungeKutta4(stiffness=K, mass=M)
elif args.time_stepping_scheme == TimeSteppingSchemes.Radau_IIA.value:
    time_stepper = RadauIIA(stiffness=K, mass=M)
else:
    raise Exception(f"Invalid time stepping scheme {args.time_stepping_scheme}. Please use one of {
                    [ts.value for ts in TimeSteppingSchemes]}")

u = u0
v = v0
a = a0
t = 0

positions_1 = [u[0]]
positions_2 = [u[1]]
times = [t]

dt = args.time_step_size

while t < T:

    # do generalized alpha step
    def f(t: float) -> ArrayLike: return np.zeros(len(u0))  # no external forces for monolithic system
    u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)

    t += dt

    u = u_new
    v = v_new
    a = a_new

    e = u[0]**2 + u[1]**2 + (u[1] - u[0])**2 + v[0]**2 + v[1]**2

    positions_1.append(u[0])
    positions_2.append(u[1])
    times.append(t)

df = pd.DataFrame()
df["times"] = times
df["errors1"] = abs(analytical_1(np.array(times)) - np.array(positions_1))
df["errors2"] = abs(analytical_2(np.array(times)) - np.array(positions_2))
df = df.set_index('times')
metadata = f'''# time_step_size: {dt}
# time stepping scheme: {args.time_stepping_scheme}
'''

participant_name = "Oscillator"

errors_csv = Path(f"errors-{participant_name}.csv")
errors_csv.unlink(missing_ok=True)

print("Error w.r.t analytical solution:")
print(f"{dt},{df['errors1'].max()},,{df['errors2'].max()}")

file: TextIOWrapper
with open(errors_csv, 'a') as file:
    file.write(f"{metadata}")
    df.to_csv(file)
