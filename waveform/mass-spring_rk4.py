from __future__ import division

import argparse
import numpy as np
import pandas as pd
import precice
import pathlib

from brot.enums import Cases, TimeSteppingSchemes, ReadWaveformSchemes, ParticipantNames, DataNames, MeshNames
from brot.output import add_metainfo
from brot.interpolation import do_linear_interpolation
from brot.timesteppers import RungeKutta4
import brot.oscillator as oscillator

this_file = pathlib.Path(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python ./solverdummy participant-name")
    quit()

participant_name = args.participantName

if participant_name == ParticipantNames.MASS_LEFT.value:
    write_data_name = DataNames.FORCE_LEFT.value
    read_data_name = DataNames.FORCE_RIGHT.value
    mesh_name = MeshNames.MASS_LEFT_MESH.value

    mass = oscillator.MassLeft.m
    stiffness = oscillator.SpringLeft.k + oscillator.SpringMiddle.k
    u0, v0, f0, d_dt_f0 = oscillator.MassLeft.u0, oscillator.MassLeft.v0, oscillator.SpringMiddle.k * oscillator.MassRight.u0, oscillator.SpringMiddle.k * oscillator.MassRight.v0
    u_analytical = oscillator.MassLeft.u_analytical
    v_analytical = oscillator.MassLeft.v_analytical

elif participant_name == ParticipantNames.MASS_RIGHT.value:
    read_data_name = DataNames.FORCE_LEFT.value
    write_data_name = DataNames.FORCE_RIGHT.value
    mesh_name = MeshNames.MASS_RIGHT_MESH.value

    mass = oscillator.MassRight.m
    stiffness = oscillator.SpringLeft.k + oscillator.SpringMiddle.k
    u0, v0, f0, d_dt_f0 = oscillator.MassRight.u0, oscillator.MassRight.v0, oscillator.SpringMiddle.k * oscillator.MassLeft.u0, oscillator.SpringMiddle.k * oscillator.MassLeft.v0
    u_analytical = oscillator.MassRight.u_analytical
    v_analytical = oscillator.MassRight.v_analytical

else:
    raise Exception(f"wrong participant name: {participant_name}. Please use one of {[p.value for p in ParticipantNames]}.")

# system:
# m ddu + k u = 0
#
# formulated as first order system
# dv = - k/m u
# du = v

ode_system = np.array([
    [0,          mass], # du
    [-stiffness, 0   ], # dv
])

num_vertices = 1  # Number of vertices

solver_process_index = 0
solver_process_size = 1

time_stepper = RungeKutta4(ode_system=ode_system)

print(f"time stepping scheme being used: {TimeSteppingSchemes.RUNGE_KUTTA_4.value}")
print(f"participant: {participant_name}")
print()
print("dt, error")

dts = [0.04, 0.02, 0.01, 0.005, 0.0025]

errors = []

for dt in dts:
    # use same dt for both solvers and preCICE
    configuration_file_name = f"configs/precice-config-{dt}.xml"
    my_dt = dt

    interface = precice.Interface(participant_name, configuration_file_name,
                                  solver_process_index, solver_process_size)

    mesh_id = interface.get_mesh_id(mesh_name)
    dimensions = interface.get_dimensions()

    vertex = np.zeros(dimensions)
    read_data = np.zeros(num_vertices)
    write_data = oscillator.SpringMiddle.k * u0 * np.ones(num_vertices)

    vertex_id = interface.set_mesh_vertex(mesh_id, vertex)
    read_data_id = interface.get_data_id(read_data_name, mesh_id)
    write_data_id = interface.get_data_id(write_data_name, mesh_id)

    precice_dt = interface.initialize()
    dt = np.min([precice_dt, my_dt])

    if interface.is_action_required(
            precice.action_write_initial_data()):
        interface.write_scalar_data(
            write_data_id, vertex_id, write_data)
        interface.mark_action_fulfilled(precice.action_write_initial_data())

    interface.initialize_data()

    # Initial Conditions

    u = u0
    v = v0
    f_start = f_end = f0
    t = 0

    positions = []
    velocities = []
    times = []

    u_write = [u]
    v_write = [v]
    t_write = [t]

    while interface.is_coupling_ongoing():
        if interface.is_action_required(
                precice.action_write_iteration_checkpoint()):
            u_cp = u
            v_cp = v
            t_cp = t
            f_start = f_end  # force at the beginning of the window
            interface.mark_action_fulfilled(
                precice.action_write_iteration_checkpoint())

            # store data for plotting and postprocessing
            positions += u_write
            velocities += v_write
            times += t_write

        read_data = interface.read_scalar_data(read_data_id, vertex_id)

        # do rk4 step
        t_stage = time_stepper.rhs_eval_points(dt)

        # implementation of waveform iteration in adapter
        f_end = read_data  # preCICE v2 returns value at end of window by default
        t_start = t_cp  # time at beginning of the window
        t_end = t_start + dt  # time at end of the window

        f = 4*[None]
        for i in range(4):
            f[i] = do_linear_interpolation(t + t_stage[i], (t_start, f_start), (t_end, f_end))

        u_new, v_new, a_new = time_stepper.do_step(u, v, None, f, dt)

        write_data = oscillator.SpringMiddle.k * u_new

        interface.write_scalar_data(
            write_data_id, vertex_id, write_data)

        precice_dt = interface.advance(dt)
        dt = np.min([precice_dt, my_dt])

        if interface.is_action_required(
                precice.action_read_iteration_checkpoint()):
            u = u_cp
            v = v_cp
            t = t_cp
            interface.mark_action_fulfilled(
                precice.action_read_iteration_checkpoint())

            # empty buffers for next window
            u_write = []
            v_write = []
            t_write = []

        else:
            u = u_new
            v = v_new
            t += dt

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
    errors.append(error)
    print(f"{dt},{error}")

df = pd.DataFrame(index=dts)
df.index.name = "dt"
df["error"] = errors

time_stepping_scheme = TimeSteppingSchemes.RUNGE_KUTTA_4.value
filepath = this_file.parent / f"{Cases.WAVEFORM.value}_{time_stepping_scheme}.csv"
df.to_csv(filepath)

add_metainfo(this_file, filepath, time_stepping_scheme, precice.__version__, ReadWaveformSchemes.LAGRANGE_LINEAR.value)
