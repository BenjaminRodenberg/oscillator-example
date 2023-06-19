from __future__ import division

import argparse
import numpy as np
from numpy.linalg import eig
import precice
import pathlib

from brot.enums import Cases, TimeSteppingSchemes, ParticipantNames, DataNames, MeshNames
from brot.output import add_metainfo
from brot.interpolation import do_lagrange_interpolation
from brot.timesteppers import GeneralizedAlpha
import brot.oscillator as oscillator

this_file = pathlib.Path(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str)
parser.add_argument("preciceConfig", nargs="?", help="precice-config.xml to be used.", type=pathlib.Path, default=this_file.parent / "configs" / "precice-config-0.04-QN.xml")
parser.add_argument("-ts", "--time-stepping", help="Time stepping scheme being used.", type=str, default=TimeSteppingSchemes.NEWMARK_BETA.value)


try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python ./solverdummy participant-name")
    quit()

participant_name = args.participantName

window_dt = 0.04

if participant_name == ParticipantNames.MASS_LEFT.value:  # left mass uses large time step size == window size
    my_dt = window_dt

    write_data_names = [DataNames.FORCE_LEFT.value]  # write data at end of single time step
    read_data_names = [DataNames.FORCE_RIGHT_1.value, DataNames.FORCE_RIGHT_2.value]  # receives results from two time steps of right mass
    mesh_name = MeshNames.MASS_LEFT_MESH.value

    mass = oscillator.MassLeft.m
    stiffness = oscillator.SpringLeft.k + oscillator.SpringMiddle.k
    u0, v0, f0, d_dt_f0 = oscillator.MassLeft.u0, oscillator.MassLeft.v0, oscillator.SpringMiddle.k * oscillator.MassRight.u0, oscillator.SpringMiddle.k * oscillator.MassRight.v0
    u_analytical = oscillator.MassLeft.u_analytical
    v_analytical = oscillator.MassLeft.v_analytical

elif participant_name == ParticipantNames.MASS_RIGHT.value:  # right mass uses small time step size == 0.5 * window size
    my_dt = window_dt * 0.5

    read_data_names = [DataNames.FORCE_LEFT.value]  # reads one piece of data corresponding to end of window
    write_data_names = [DataNames.FORCE_RIGHT_1.value, DataNames.FORCE_RIGHT_2.value]  # sends results for each of the two substeps
    mesh_name = MeshNames.MASS_RIGHT_MESH.value

    mass = oscillator.MassRight.m
    stiffness = oscillator.SpringLeft.k + oscillator.SpringMiddle.k
    u0, v0, f0, d_dt_f0 = oscillator.MassRight.u0, oscillator.MassRight.v0, oscillator.SpringMiddle.k * oscillator.MassLeft.u0, oscillator.SpringMiddle.k * oscillator.MassLeft.v0
    u_analytical = oscillator.MassRight.u_analytical
    v_analytical = oscillator.MassRight.v_analytical

else:
    raise Exception(f"wrong participant name: {participant_name}. Please use one of {[p.value for p in ParticipantNames]}.")

num_vertices = 1  # Number of vertices

solver_process_index = 0
solver_process_size = 1

configuration_file_name = str(args.preciceConfig)

interface = precice.Interface(participant_name, configuration_file_name,
                            solver_process_index, solver_process_size)

mesh_id = interface.get_mesh_id(mesh_name)
dimensions = interface.get_dimensions()

vertex = np.zeros(dimensions)
read_data = np.zeros(num_vertices)
write_data = oscillator.SpringMiddle.k * u0 * np.ones(num_vertices)

vertex_id = interface.set_mesh_vertex(mesh_id, vertex)
read_data_ids = [interface.get_data_id(read_data_name, mesh_id) for read_data_name in read_data_names]
write_data_ids = [interface.get_data_id(write_data_name, mesh_id) for write_data_name in write_data_names]

precice_dt = interface.initialize()
dt = np.min([precice_dt, my_dt])

if interface.is_action_required(precice.action_write_initial_data()):
    for write_data_id in write_data_ids:
        interface.write_scalar_data(write_data_id, vertex_id, write_data)
    interface.mark_action_fulfilled(precice.action_write_initial_data())

interface.initialize_data()

# Initial Conditions

a0 = (f0 - stiffness * u0) / mass
u = u0
v = v0
a = a0
f_start = f_end = f0
t = 0

# Generalized Alpha Parameters
if args.time_stepping == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.4, alpha_m=0.2)
elif args.time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.0, alpha_m=0.0)
elif args.time_stepping == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    raise Exception(f"Please use monolithic_rk4.py for using --time-stepping=\"{args.time_stepping}\"")
else:
    raise Exception(f"Invalid time stepping scheme {args.time_stepping}. Please use one of {[ts.value for ts in TimeSteppingSchemes]}")

positions = []
velocities = []
times = []

u_write = [u]
v_write = [v]
t_write = [t]

substep = 0

while interface.is_coupling_ongoing():
    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
        u_cp = u
        v_cp = v
        a_cp = a
        t_cp = t
        f_start = f_end  # force at the beginning of the window
        substep = 0
        interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

        # store data for plotting and postprocessing
        positions += u_write
        velocities += v_write
        times += t_write

    t_start = t_cp  # time at beginning of the window
    t_end = t_start + window_dt  # time at end of the window

    if participant_name == ParticipantNames.MASS_LEFT.value:  # does single large time step per window
        f_mid = interface.read_scalar_data(read_data_ids[0], vertex_id)
        t_mid = (t_start + t_end) * 0.5
        f_end = interface.read_scalar_data(read_data_ids[1], vertex_id)
        ts = [t_start, t_mid, t_end]
        fs = [f_start, f_mid, f_end]

    elif participant_name == ParticipantNames.MASS_RIGHT.value:  # does two small time steps per window
        f_end = interface.read_scalar_data(read_data_ids[0], vertex_id)  # read data always corresponds to end of window
        ts = [t_start, t_end]
        fs = [f_start, f_end]

    # do time stepping
    t_f = time_stepper.rhs_eval_points(dt)
    f = do_lagrange_interpolation(t + t_f, ts, fs)
    u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)

    write_data = oscillator.SpringMiddle.k * u_new

    if participant_name == ParticipantNames.MASS_LEFT.value:  # does two substeps per window
        interface.write_scalar_data(write_data_ids[0], vertex_id, write_data)
    elif participant_name == ParticipantNames.MASS_RIGHT.value:  # does one step per window
        # TODO: This is wrong! We should only write ALL data before the last advance of the window, otherwise data will be reset for subsequent advance calls.
        interface.write_scalar_data(write_data_ids[substep], vertex_id, write_data)

    precice_dt = interface.advance(dt)
    dt = np.min([precice_dt, my_dt])

    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
        u = u_cp
        v = v_cp
        a = a_cp
        t = t_cp
        substep = 0

        interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())

        # empty buffers for next window
        u_write = []
        v_write = []
        t_write = []

    else:
        u = u_new
        v = v_new
        a = a_new
        t += dt
        substep += 1

        # write data to buffers
        u_write.append(u)
        v_write.append(v)
        t_write.append(t)

positions += u_write
velocities += v_write
times += t_write

interface.finalize()

# print errors
# analytic solution is only valid for this setup!
error = np.max(abs(u_analytical(np.array(times))-np.array(positions)))
print(f"{dt},{error}")
# todo: write into file!
# todo: name file in a structured way!
