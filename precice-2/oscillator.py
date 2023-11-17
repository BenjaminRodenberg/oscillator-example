from __future__ import division

import argparse
import numpy as np
from numpy.linalg import eig
import precice
from enum import Enum
import pandas as pd
from pathlib import Path

import problemDefinition

from brot.enums import TimeSteppingSchemes
import brot.timesteppers as timeSteppers


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

if participant_name == Participant.MASS_LEFT.value:
    write_data_name = 'Force-Left'
    read_data_name = 'Force-Right'
    mesh_name = 'Mass-Left-Mesh'

    this_mass = problemDefinition.MassLeft
    this_spring = problemDefinition.SpringLeft
    connecting_spring = problemDefinition.SpringMiddle
    other_mass = problemDefinition.MassRight

elif participant_name == Participant.MASS_RIGHT.value:
    read_data_name = 'Force-Left'
    write_data_name = 'Force-Right'
    mesh_name = 'Mass-Right-Mesh'

    this_mass = problemDefinition.MassRight
    this_spring = problemDefinition.SpringRight
    connecting_spring = problemDefinition.SpringMiddle
    other_mass = problemDefinition.MassLeft

else:
    raise Exception(f"wrong participant name: {participant_name}")

mass = this_mass.m
stiffness = this_spring.k + connecting_spring.k
u0, v0, f0, d_dt_f0 = this_mass.u0, this_mass.v0, connecting_spring.k * other_mass.u0, connecting_spring.k * other_mass.v0

num_vertices = 1  # Number of vertices

solver_process_index = 0
solver_process_size = 1

configuration_file_name = "precice-config.xml"

participant = precice.Interface(participant_name, configuration_file_name, solver_process_index, solver_process_size)

mesh_id = participant.get_mesh_id(mesh_name)
dimensions = participant.get_dimensions()

vertex = np.zeros(dimensions)
read_data = np.zeros(num_vertices)
write_data = connecting_spring.k * u0 * np.ones(num_vertices)

vertex_id = participant.set_mesh_vertex(mesh_id, vertex)
read_data_id = participant.get_data_id(read_data_name, mesh_id)
write_data_id = participant.get_data_id(write_data_name, mesh_id)

precice_dt = participant.initialize()
my_dt = precice_dt / args.n_substeps  # use my_dt < precice_dt for subcycling
dt = np.min([precice_dt, my_dt])

if participant.is_action_required(
        precice.action_write_initial_data()):
    participant.write_scalar_data(
        write_data_id, vertex_id, write_data)
    participant.mark_action_fulfilled(precice.action_write_initial_data())

participant.initialize_data()

# Initial Conditions
a0 = (f0 - stiffness * u0) / mass
u = u0
v = v0
a = a0
t = 0

if args.time_stepping == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    time_stepper = timeSteppers.GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.4, alpha_m=0.2)
elif args.time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    time_stepper = timeSteppers.GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.0, alpha_m=0.0)
elif args.time_stepping == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    ode_system = np.array([
        [0, mass],  # du
        [-stiffness, 0],  # dv
    ])
    time_stepper = timeSteppers.RungeKutta4(ode_system=ode_system)
elif args.time_stepping == TimeSteppingSchemes.Radau_IIA.value:
    ode_system = np.array([
        [0, mass],  # du
        [-stiffness, 0],  # dv
    ])
    time_stepper = timeSteppers.RadauIIA(ode_system=ode_system)
else:
    raise Exception(
        f"Invalid time stepping scheme {args.time_stepping}. Please use one of {[ts.value for ts in TimeSteppingSchemes]}")


positions = []
velocities = []
times = []

u_write = [u]
v_write = [v]
t_write = [t]

while participant.is_coupling_ongoing():
    if participant.is_action_required(
            precice.action_write_iteration_checkpoint()):
        u_cp = u
        v_cp = v
        a_cp = a
        t_cp = t
        # store data for plotting and postprocessing
        positions += u_write
        velocities += v_write
        times += t_write

        participant.mark_action_fulfilled(
            precice.action_write_iteration_checkpoint())

    read_data = participant.read_scalar_data(read_data_id, vertex_id)
    f = read_data

    # do time stepping
    u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)
    t_new = t + dt

    write_data = connecting_spring.k * u_new

    participant.write_scalar_data(
        write_data_id, vertex_id, write_data)

    precice_dt = participant.advance(dt)
    dt = np.min([precice_dt, my_dt])

    if participant.is_action_required(
            precice.action_read_iteration_checkpoint()):
        u = u_cp
        v = v_cp
        a = a_cp
        t = t_cp
        # empty buffers for next window
        u_write = []
        v_write = []
        t_write = []

        participant.mark_action_fulfilled(
            precice.action_read_iteration_checkpoint())

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
df["errors"] = abs(this_mass.u_analytical(np.array(times)) - np.array(positions))
df = df.set_index('times')
metadata = f'''# time_window_size: {precice_dt}
# time_step_size: {my_dt}
'''

errors_csv = Path(f"errors-{participant_name}.csv")
errors_csv.unlink(missing_ok=True)

print("Error w.r.t analytical solution:")
print(f"{my_dt},{df['errors'].max()}")

with open(errors_csv, 'a') as f:
    f.write(f"{metadata}")
    df.to_csv(f)

# output trajectory
trajectory_df = pd.DataFrame()
trajectory_df["time"] = times
trajectory_df["position"] = positions
trajectory_df["velocity"] = velocities

trajectory_csv = Path(f"trajectory-{participant_name}.csv")
trajectory_csv.unlink(missing_ok=True)

with open(trajectory_csv, 'a') as f:
    trajectory_df.to_csv(f)