import argparse
import numpy as np
import pandas as pd
from pathlib import Path

from brot.timeSteppersMonolithic import TimeStepper, TimeSteppingSchemes, GeneralizedAlpha, RungeKutta4, RadauIIA
import brot.oscillator as problemDefinition

from io import TextIOWrapper
from numpy.typing import ArrayLike

parser = argparse.ArgumentParser()
parser.add_argument(
    "-tss",
    "--time-stepping-scheme",
    help=f"Time stepping scheme being used. Please use one of {[ts.value for ts in TimeSteppingSchemes]}",
    type=str,
    default=TimeSteppingSchemes.NEWMARK_BETA.value)
parser.add_argument("-dt", "--time-step-size", help=f"Time step size being used", type=float, default=0.04)
args = parser.parse_args()

M = problemDefinition.M
K = problemDefinition.K

mass_left = problemDefinition.MassLeft
mass_right = problemDefinition.MassRight

u0 = np.array([mass_left.u0, mass_right.u0])
v0 = np.array([mass_left.v0, mass_right.v0])
a0 = -K.dot(np.linalg.inv(M).dot(u0))

# Initial Conditions
u = u0
v = v0
a = a0
t = 0

T = 1

time_stepper: TimeStepper

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

positions_left = [u[0]]
positions_right = [u[1]]
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

    positions_left.append(u[0])
    positions_right.append(u[1])
    times.append(t)

df = pd.DataFrame()
df["times"] = times
df["error Mass-Left"] = abs(mass_left.u_analytical(np.array(times)) - np.array(positions_left))
df["error Mass-Right"] = abs(mass_right.u_analytical(np.array(times)) - np.array(positions_right))
df = df.set_index('times')
metadata = f'''# time_step_size: {dt}
# time stepping scheme: {args.time_stepping_scheme}
'''

participant_name = "Monolithic"

output_csv = Path(f"output-{participant_name}.csv")
output_csv.unlink(missing_ok=True)

print("Error w.r.t analytical solution:")
print(f"{dt},{df['error Mass-Left'].max()},,{df['error Mass-Right'].max()}")

file: TextIOWrapper
with open(output_csv, 'a') as file:
    file.write(f"{metadata}")
    df.to_csv(file)
