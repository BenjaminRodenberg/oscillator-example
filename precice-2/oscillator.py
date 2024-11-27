import argparse
import numpy as np
import pandas as pd
from pathlib import Path

import precice
from enum import Enum
from typing import Type

from brot.timeSteppers import TimeStepper, TimeSteppingSchemes, GeneralizedAlpha, RungeKutta4, RadauIIA
import brot.oscillator as problemDefinition

from io import TextIOWrapper

class Participant(Enum):
    MASS_LEFT = "Mass-Left"
    MASS_RIGHT = "Mass-Right"


parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str, choices=[p.value for p in Participant])
parser.add_argument(
    "-ts",
    "--time-stepping",
    help="Time stepping scheme being used.",
    type=str,
    choices=[
        s.value for s in TimeSteppingSchemes],
    default=TimeSteppingSchemes.NEWMARK_BETA.value)
parser.add_argument(
    "-s",
    "--n-substeps",
    help="Number of substeps in one window for this participant",
    type=int,
    default=1)
args = parser.parse_args()

participant_name = args.participantName

this_mass: Type[problemDefinition.Mass]
other_mass: Type[problemDefinition.Mass]
this_spring: Type[problemDefinition.Spring]
connecting_spring = problemDefinition.SpringMiddle

if participant_name == Participant.MASS_LEFT.value:
    write_data_name = 'Displacement-Left'
    read_data_name = 'Displacement-Right'
    mesh_name = 'Mass-Left-Mesh'

    this_mass = problemDefinition.MassLeft
    this_spring = problemDefinition.SpringLeft
    other_mass = problemDefinition.MassRight

elif participant_name == Participant.MASS_RIGHT.value:
    read_data_name = 'Displacement-Left'
    write_data_name = 'Displacement-Right'
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
read_data_id = participant.get_data_id(read_data_name, mesh_id)
write_data_id = participant.get_data_id(write_data_name, mesh_id)

precice_dt = participant.initialize()
my_dt = precice_dt / args.n_substeps  # use my_dt < precice_dt for subcycling
dt = np.min([precice_dt, my_dt])

if participant.is_action_required(precice.action_write_initial_data()):
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

while participant.is_coupling_ongoing():
    if participant.is_action_required(precice.action_write_iteration_checkpoint()):
        u_cp = u
        v_cp = v
        a_cp = a
        t_cp = t
        # store data for plotting and postprocessing
        positions += u_write
        velocities += v_write
        times += t_write

        participant.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

    f = connecting_spring.k * participant.read_scalar_data(read_data_id, vertex_id)

    # do time stepping
    u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)
    t_new = t + dt

    participant.write_scalar_data(write_data_id, vertex_id, u_new)

    precice_dt = participant.advance(dt)
    dt = np.min([precice_dt, my_dt])

    if participant.is_action_required(precice.action_read_iteration_checkpoint()):
        u = u_cp
        v = v_cp
        a = a_cp
        t = t_cp
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