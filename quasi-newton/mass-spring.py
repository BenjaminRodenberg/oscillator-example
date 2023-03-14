from __future__ import division

import argparse
import numpy as np
from numpy.linalg import eig
import precice
import pathlib

from helpers.enums import Cases, TimeSteppingSchemes, ParticipantNames, DataNames, MeshNames
from helpers.output import add_metainfo
from helpers.interpolation import do_linear_interpolation, do_lagrange_interpolation

this_file = pathlib.Path(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str)
parser.add_argument("-ts", "--time-stepping", help="Time stepping scheme being used.", type=str, default=TimeSteppingSchemes.NEWMARK_BETA.value)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python ./solverdummy participant-name")
    quit()

participant_name = args.participantName

m_1, m_2 = 1, 1
k_1, k_2, k_12 = 4*np.pi**2, 4*np.pi**2, 16*(np.pi**2)

# system:
# m ddu + k u = f

M = np.array([[m_1, 0], [0, m_2]])
M_inv = np.array([[1 / m_1, 0], [0, 1 / m_2]])
K = np.array([[k_1 + k_12, -k_12], [-k_12, k_2 + k_12]])

eigenvalues, eigenvectors = eig(K)  # should be K/M?
omega = np.sqrt(eigenvalues)
A, B = eigenvectors

u0_1 = 1
u0_2 = 0
v0_1 = 0
v0_2 = 0

window_dt = 0.04

if participant_name == ParticipantNames.MASS_LEFT.value:  # left mass uses large time step size == window size
    my_dt = window_dt

    write_data_names = [DataNames.FORCE_LEFT.value]  # write data at end of single time step
    read_data_names = [DataNames.FORCE_RIGHT_1.value, DataNames.FORCE_RIGHT_2.value]  # receives results from two time steps of right mass
    mesh_name = MeshNames.MASS_LEFT_MESH.value

    mass = m_1
    stiffness = k_1 + k_12
    u0, v0, f0, d_dt_f0 = u0_1, v0_1, k_12 * u0_2, k_12 * v0_2
    u_analytical = lambda t: .5 * (np.cos(2 * np.pi * t) + np.cos(6 * np.pi * t))
    v_analytical = lambda t: .5 * (-2 * np.pi * np.sin(2 * np.pi * t) - 6 * np.pi * np.sin(6 * np.pi * t))

elif participant_name == ParticipantNames.MASS_RIGHT.value:  # right mass uses small time step size == 0.5 * window size
    my_dt = window_dt * 0.5

    read_data_names = [DataNames.FORCE_LEFT.value]  # reads one piece of data corresponding to end of window
    write_data_names = [DataNames.FORCE_RIGHT_1.value, DataNames.FORCE_RIGHT_2.value]  # sends results for each of the two substeps
    mesh_name = MeshNames.MASS_RIGHT_MESH.value

    mass = m_2
    stiffness = k_2 + k_12
    u0, v0, f0, d_dt_f0 = u0_2, v0_2, k_12 * u0_1, k_12 * v0_1
    u_analytical = lambda t: .5 * (np.cos(2 * np.pi * t) - np.cos(6 * np.pi * t))
    v_analytical = lambda t: .5 * (-2 * np.pi * np.sin(2 * np.pi * t) + 6 * np.pi * np.sin(6 * np.pi * t))

else:
    raise Exception(f"wrong participant name: {participant_name}. Please use one of {[p.value for p in ParticipantNames]}.")

num_vertices = 1  # Number of vertices

solver_process_index = 0
solver_process_size = 1

configuration_file_name = f"precice-config-{window_dt}.xml"

interface = precice.Interface(participant_name, configuration_file_name,
                            solver_process_index, solver_process_size)

mesh_id = interface.get_mesh_id(mesh_name)
dimensions = interface.get_dimensions()

vertex = np.zeros(dimensions)
read_data = np.zeros(num_vertices)
write_data = k_12 * u0 * np.ones(num_vertices)

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
    alpha_f = 0.4
    alpha_m = 0.2
elif args.time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    alpha_f = 0.0
    alpha_m = 0.0
elif args.time_stepping == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    raise Exception(f"Please use monolithic_rk4.py for using --time-stepping=\"{args.time_stepping}\"")
else:
    raise Exception(f"Invalid time stepping scheme {args.time_stepping}. Please use one of {[ts.value for ts in TimeSteppingSchemes]}")

gamma = 0.5 - alpha_m + alpha_f
beta = 0.25 * (gamma + 0.5)

m = 3*[None]
m[0] = (1-alpha_m)/(beta*dt**2)
m[1] = (1-alpha_m)/(beta*dt)
m[2] = (1-alpha_m-2*beta)/(2*beta)
k_bar = stiffness * (1 - alpha_f) + m[0] * mass

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

    t_f = (1-alpha_f) * dt

    t_start = t_cp  # time at beginning of the window
    t_end = t_start + window_dt  # time at end of the window

    if participant_name == ParticipantNames.MASS_LEFT.value:  # does single large time step per window
        f_mid = interface.read_scalar_data(read_data_ids[0], vertex_id)
        t_mid = (t_start + t_end) * 0.5
        f_end = interface.read_scalar_data(read_data_ids[1], vertex_id)
        f = do_lagrange_interpolation(t + t_f, [t_start, t_mid, t_end], [f_start, f_mid, f_end])

    elif participant_name == ParticipantNames.MASS_RIGHT.value:  # does two small time steps per window
        f_end = interface.read_scalar_data(read_data_ids[0], vertex_id)  # read data always corresponds to end of window
        print(f"evaluation time for substep {substep} is {t + t_f}")
        f = do_linear_interpolation(t + t_f, (t_start, f_start), (t_end, f_end))
        print(f"f={f}")

    # do generalized alpha step
    u_new = (f - alpha_f * stiffness * u + mass*(m[0]*u + m[1]*v + m[2]*a)) / k_bar
    a_new = 1.0 / (beta * dt**2) * (u_new - u - dt * v) - (1-2*beta) / (2*beta) * a
    v_new = v + dt * ((1-gamma)*a+gamma*a_new)

    write_data = k_12 * u_new

    if participant_name == ParticipantNames.MASS_LEFT.value:  # does two substeps per window
        interface.write_scalar_data(write_data_ids[0], vertex_id, write_data)
    elif participant_name == ParticipantNames.MASS_RIGHT.value:  # does one step per window
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
