import argparse
import numpy as np
import pandas as pd
from pathlib import Path

import precice
from enum import Enum
from typing import Type

from brot.timeSteppers import TimeStepper, TimeSteppingSchemes, GeneralizedAlpha, RungeKutta4, RadauIIA
import brot.oscillator as problemDefinition
from brot.interpolation import InterpolationSchemes, do_lagrange_interpolation

from io import TextIOWrapper

class Participant(Enum):
    MASS_LEFT = "Mass-Left"
    MASS_RIGHT = "Mass-Right"


n_substeps_default = 1


parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str, choices=[p.value for p in Participant])
parser.add_argument("-ts", "--time-stepping", help="Time stepping scheme being used.", type=str, choices=[s.value for s in TimeSteppingSchemes], default=TimeSteppingSchemes.NEWMARK_BETA.value)
parser.add_argument("-nl", "--n-substeps-left", help="Number of substeps in one window for Mass-Left", type=int, default=n_substeps_default)
parser.add_argument("-nr", "--n-substeps-right", help="Number of substeps in one window for Mass-Right", type=int, default=n_substeps_default)
parser.add_argument("-is", "--interpolation-scheme", help=f"Interpolation scheme being used.", type=str, choices=[InterpolationSchemes.LAGRANGE.value, InterpolationSchemes.BSPLINE.value], default=InterpolationSchemes.LAGRANGE.value)
parser.add_argument("-p", "--interpolation-degree", help="Desired degree of interpolation scheme (Only allowed, if using BSpline interpolation).", type=int, default=1)
args = parser.parse_args()

participant_name = args.participantName

this_mass: Type[problemDefinition.Mass]
other_mass: Type[problemDefinition.Mass]
this_spring: Type[problemDefinition.Spring]
connecting_spring = problemDefinition.SpringMiddle

write_data_names = []
read_data_names = []
if participant_name == Participant.MASS_LEFT.value:
    # this participant must know substeps of other participant to use correct time grid in interpolation (breaks black-box!)
    n_substeps_this = args.n_substeps_left
    n_substeps_other = args.n_substeps_right

    for i in range(n_substeps_this):
        write_data_names.append(f"Displacement-Left-{i+1}")

    for i in range(n_substeps_other):
        read_data_names.append(f"Displacement-Right-{i+1}")

    mesh_name = 'Mass-Left-Mesh'

    this_mass = problemDefinition.MassLeft
    this_spring = problemDefinition.SpringLeft
    other_mass = problemDefinition.MassRight

elif participant_name == Participant.MASS_RIGHT.value:
    # this participant must know substeps of other participant to use correct time grid in interpolation (breaks black-box!)
    n_substeps_this = args.n_substeps_right
    n_substeps_other = args.n_substeps_left

    for i in range(n_substeps_this):
        write_data_names.append(f"Displacement-Right-{i+1}")

    for i in range(n_substeps_other):
        read_data_names.append(f"Displacement-Left-{i+1}")

    mesh_name = 'Mass-Right-Mesh'

    this_mass = problemDefinition.MassRight
    this_spring = problemDefinition.SpringRight
    other_mass = problemDefinition.MassLeft

else:
    raise Exception(f"wrong participant name: {participant_name}")

mass = this_mass.m
stiffness = this_spring.k + connecting_spring.k
u0, v0, f0 = this_mass.u0, this_mass.v0, connecting_spring.k * other_mass.u0

solver_process_index = 0
solver_process_size = 1

configuration_file_name = "precice-config.xml"

participant = precice.Interface(participant_name, configuration_file_name, solver_process_index, solver_process_size)

mesh_id = participant.get_mesh_id(mesh_name)
dimensions = participant.get_dimensions()

vertex = np.zeros(dimensions)
vertex_id = participant.set_mesh_vertex(mesh_id, vertex)
read_data_ids = [participant.get_data_id(read_data_name, mesh_id) for read_data_name in read_data_names]
write_data_ids = [participant.get_data_id(write_data_name, mesh_id) for write_data_name in write_data_names]

precice_dt = participant.initialize()
my_dt = precice_dt / n_substeps_this
other_dt = precice_dt / n_substeps_other
dt = np.min([precice_dt, my_dt])

if participant.is_action_required(precice.action_write_initial_data()):
    for write_data_id in write_data_ids:
        participant.write_scalar_data(write_data_id, vertex_id, u0)
    participant.mark_action_fulfilled(precice.action_write_initial_data())

participant.initialize_data()

# Initial Conditions
a0 = (f0 - stiffness * u0) / mass
u = u0
v = v0
a = a0
t = 0

time_stepper: TimeStepper

if args.time_stepping == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.4, alpha_m=0.2)
elif args.time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.0, alpha_m=0.0)
elif args.time_stepping == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    time_stepper = RungeKutta4(stiffness=stiffness, mass=mass)
elif args.time_stepping == TimeSteppingSchemes.Radau_IIA.value:
    time_stepper = RadauIIA(stiffness=stiffness, mass=mass)
else:
    raise Exception(f"Invalid time stepping scheme {args.time_stepping}. Please use one of {[ts.value for ts in TimeSteppingSchemes]}")


positions = []
velocities = []
times = []

u_write = [u]
v_write = [v]
t_write = [t]

substep = 1
write_buffer = {}  # internal buffer for write data until final advance of window is called
u_read = (n_substeps_other + 1) * [other_mass.u0]

while participant.is_coupling_ongoing():
    if participant.is_action_required(precice.action_write_iteration_checkpoint()):
        u_cp = u
        v_cp = v
        a_cp = a
        t_cp = t
        u_read[0] = u_read[-1]  # force at the beginning of the new window is force at end of last window
        t_read = [t + i*other_dt for i in range(n_substeps_other+1)]
        substep = 0

        # store data for plotting and postprocessing
        positions += u_write
        velocities += v_write
        times += t_write

        participant.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

    # read n_substeps_other samples that will be associated with t_read
    for i in range(n_substeps_other):
        u_read[i+1] = participant.read_scalar_data(read_data_ids[i], vertex_id)

    # implementation of waveform iteration in adapter
    if args.interpolation_scheme == InterpolationSchemes.LAGRANGE.value:
        interpolant = lambda t: do_lagrange_interpolation(t, t_read, u_read)
    elif args.interpolation_scheme == InterpolationSchemes.BSPLINE.value:
        from scipy.interpolate import splrep, splev
        b_spline_degree = args.interpolation_degree
        tck = splrep(t_read, u_read, k=b_spline_degree)
        interpolant = lambda t: splev(t, tck)

    f = lambda dt: connecting_spring.k * interpolant(t + dt)

    # do time stepping
    u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)
    t_new = t + dt

    # store result of this substep to buffer
    write_buffer[write_data_ids[substep]] = u_new
    substep += 1

    # write n_substeps_this samples to other participant
    if substep == n_substeps_this:  # this time step concludes window. Write data
        for id, data in write_buffer.items():
            participant.write_scalar_data(id, vertex_id, data)

    precice_dt = participant.advance(dt)
    dt = np.min([precice_dt, my_dt])

    if participant.is_action_required(precice.action_read_iteration_checkpoint()):
        u = u_cp
        v = v_cp
        a = a_cp
        t = t_cp
        substep = 0
        # empty buffers for next window
        u_write = []
        v_write = []
        t_write = []

        participant.mark_action_fulfilled(precice.action_read_iteration_checkpoint())

    else:
        u = u_new
        v = v_new
        a = a_new
        t = t_new
        # write data to buffers
        u_write.append(u)
        v_write.append(v)
        t_write.append(t)

# store final result
positions += u_write
velocities += v_write
times += t_write

participant.finalize()

df = pd.DataFrame()
df["times"] = times
df["position"] = positions
df["velocity"] = velocities
df["errors"] = abs(this_mass.u_analytical(np.array(times)) - np.array(positions))
df = df.set_index('times')
metadata = f'''# time_window_size: {precice_dt}
# time_step_size: {my_dt}
# time stepping scheme: {args.time_stepping}
'''

errors_csv = Path(f"output-{participant_name}.csv")
errors_csv.unlink(missing_ok=True)

print("Error w.r.t analytical solution:")
print(f"{my_dt},{df['errors'].max()}")

file: TextIOWrapper
with open(errors_csv, 'a') as file:
    file.write(f"{metadata}")
    df.to_csv(file)