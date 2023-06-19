from __future__ import division

import argparse
import numpy as np
import pandas as pd
import precice
import pathlib

from brot.enums import Cases, TimeSteppingSchemes, ReadWaveformSchemes, MultirateMode, ParticipantNames, DataNames, MeshNames
from brot.output import add_metainfo
from brot.interpolation import do_linear_interpolation, do_lagrange_interpolation
from brot.timesteppers import GeneralizedAlpha, RungeKutta4
import brot.oscillator as oscillator

this_file = pathlib.Path(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str)
parser.add_argument("-ts", "--time-stepping", help="Time stepping scheme being used.", type=str, default=TimeSteppingSchemes.NEWMARK_BETA.value)
parser.add_argument("-is", "--interpolation-scheme", help="Interpolation scheme being used.", type=str, default=ReadWaveformSchemes.LAGRANGE_QUADRATIC.value)
parser.add_argument("-mr", "--multirate", help="Pick one of 3 modes for multirate", type=str, default=MultirateMode.TWO_SUBSTEPS.value)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python ./solverdummy participant-name")
    quit()

participant_name = args.participantName

if participant_name == ParticipantNames.MASS_LEFT.value:
    write_data_names = [DataNames.FORCE_LEFT_1.value, DataNames.FORCE_LEFT_2.value]  # sends results for each of the two substeps
    read_data_names = [DataNames.FORCE_RIGHT_1.value, DataNames.FORCE_RIGHT_2.value]  # receives results from two time steps of right mass
    mesh_name = MeshNames.MASS_LEFT_MESH.value

    mass = oscillator.MassLeft.m
    stiffness = oscillator.SpringLeft.k + oscillator.SpringMiddle.k
    u0, v0, f0, d_dt_f0 = oscillator.MassLeft.u0, oscillator.MassLeft.v0, oscillator.SpringMiddle.k * oscillator.MassRight.u0, oscillator.SpringMiddle.k * oscillator.MassRight.v0
    u_analytical = oscillator.MassLeft.u_analytical
    v_analytical = oscillator.MassLeft.v_analytical

elif participant_name == ParticipantNames.MASS_RIGHT.value:
    read_data_names = [DataNames.FORCE_LEFT_1.value, DataNames.FORCE_LEFT_2.value]  # receives results from two time steps of left mass
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

if args.time_stepping == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.4, alpha_m=0.2)
elif args.time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.0, alpha_m=0.0)
elif args.time_stepping == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    ode_system = np.array([
        [0,          mass], # du
        [-stiffness, 0   ], # dv
    ])
    time_stepper = RungeKutta4(ode_system=ode_system)
else:
    raise Exception(f"Invalid time stepping scheme {args.time_stepping}. Please use one of {[ts.value for ts in TimeSteppingSchemes]}")

print(f"time stepping scheme being used: {args.time_stepping}")
print(f"participant: {participant_name}")
print()
print("configured_precice_dt, my_dt, error")

dts = [0.04, 0.02, 0.01, 0.005, 0.0025]

if args.interpolation_scheme == ReadWaveformSchemes.LAGRANGE_LINEAR.value:
    interpolation_scheme = ReadWaveformSchemes.LAGRANGE_LINEAR.value
elif args.interpolation_scheme == ReadWaveformSchemes.LAGRANGE_QUADRATIC.value:
    interpolation_scheme = ReadWaveformSchemes.LAGRANGE_QUADRATIC.value
else:
    raise Exception(f"wrong interpolation scheme name: {args.interpolation_scheme}. Please use one of {[p.value for p in ReadWaveformSchemes]}.")

errors = []
my_dts = []

for dt in dts:

    if args.multirate == MultirateMode.TWO_SUBSTEPS.value:
        configured_precice_dt = dt
        my_dt = dt/2
    else:
        raise Exception(f"wrong multirate mode: {args.multirate}. Please use one of {[m.value for m in MultirateMode]}.")

    if interpolation_scheme == ReadWaveformSchemes.LAGRANGE_QUADRATIC.value:
      configuration_file_name = f"lagrange_2_configs/precice-config-{configured_precice_dt}.xml"


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

    positions = []
    velocities = []
    times = []

    u_write = [u]
    v_write = [v]
    t_write = [t]

    substep = 0

    write_buffer = {}  # use this dict to buffer write data until final advance of window is called

    while interface.is_coupling_ongoing():
        if interface.is_action_required(
                precice.action_write_iteration_checkpoint()):
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
        t_end = t_start + configured_precice_dt  # time at end of the window

        f_mid = interface.read_scalar_data(read_data_ids[0], vertex_id)
        t_mid = (t_start + t_end) * 0.5

        f_end = interface.read_scalar_data(read_data_ids[1], vertex_id)

        ts = [t_start, t_mid, t_end]
        fs = [f_start, f_mid, f_end]

        # do time stepping
        t_f = time_stepper.rhs_eval_points(dt)

        f = len(t_f)*[None]
        for i in range(len(f)):
            '''Lagrange Interpolation
            f[i] = do_lagrange_interpolation(t + t_f[i], ts, fs)
            '''

            '''scipy Bspline
            '''
            from scipy.interpolate import splrep, splev
            tck = splrep(ts, fs, k=2)
            interpolant = lambda t: splev(t, tck)
            f[i] = interpolant(t + t_f[i])

        u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)

        write_data = oscillator.SpringMiddle.k * u_new

        write_buffer[write_data_ids[substep]] = write_data

        if dt == precice_dt:  # this time step concludes window. Write data
            for id, data in write_buffer.items():
                interface.write_scalar_data(id, vertex_id, data)


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
    errors.append(error)
    my_dts.append(my_dt)
    print(f"{configured_precice_dt}, {my_dt}, {error}")

df = pd.DataFrame(index=dts)
df.index.name = "dt"
df["my_dt"] = my_dts
df["error"] = errors

time_stepping_scheme = args.time_stepping
filepath = this_file.parent / f"{Cases.MULTIRATE.value}_{participant_name}_{time_stepping_scheme}_{interpolation_scheme}_{args.multirate}.csv"
df.to_csv(filepath)

add_metainfo(this_file, filepath, time_stepping_scheme, precice.__version__, interpolation_scheme)
