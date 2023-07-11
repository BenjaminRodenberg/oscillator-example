from __future__ import division

import argparse
import numpy as np
import pandas as pd
import precice
import pathlib
from configs.create_config import render

from brot.enums import Cases, TimeSteppingSchemes, ReadWaveformSchemes, ParticipantNames, DataNames, MeshNames, AccelerationSchemes
from brot.output import add_metainfo
from brot.timesteppers import GeneralizedAlpha, RungeKutta4
import brot.oscillator as oscillator

this_file = pathlib.Path(__file__)

degree_default = 3
n_substeps_default = 4
substeps_default = True

parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str)
parser.add_argument("-tsl", "--time-stepping-left", help="Time stepping scheme being used for Mass-Left.", type=str, default=TimeSteppingSchemes.RUNGE_KUTTA_4.value)
parser.add_argument("-tsr", "--time-stepping-right", help="Time stepping scheme being used for Mass-Right.", type=str, default=TimeSteppingSchemes.RUNGE_KUTTA_4.value)
parser.add_argument("-acc", "--acceleration-scheme", help="Acceleration scheme being used", type=str, default=AccelerationSchemes.NONE.value)
parser.add_argument("-is", "--interpolation-scheme", help="Interpolation scheme being used.", type=str, default=ReadWaveformSchemes.ZERO.value)
parser.add_argument("-p", "--interpolation-degree", help="Desired degree of interpolation scheme.", type=int, default=degree_default)
parser.add_argument("-s", "--no-substeps", help="Set substeps=\"false\" precice-config.", action='store_const', dest="substeps", default=substeps_default, const=not(substeps_default))
parser.add_argument("-nl", "--n-substeps-left", help="Number of substeps in one window for Mass-Left", type=int, default=n_substeps_default)
parser.add_argument("-nr", "--n-substeps-right", help="Number of substeps in one window for Mass-Right", type=int, default=n_substeps_default)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python ./solverdummy participant-name")
    quit()

participant_name = args.participantName

if participant_name == ParticipantNames.MASS_LEFT.value:
    n_substeps_this = args.n_substeps_left
    n_substeps_other = args.n_substeps_right

    time_stepping = args.time_stepping_left

    write_data_name = DataNames.FORCE_LEFT.value
    read_data_name = DataNames.FORCE_RIGHT.value
    mesh_name = MeshNames.MASS_LEFT_MESH.value

    mass = oscillator.MassLeft.m
    stiffness = oscillator.SpringLeft.k + oscillator.SpringMiddle.k
    u0, v0, f0, d_dt_f0 = oscillator.MassLeft.u0, oscillator.MassLeft.v0, oscillator.SpringMiddle.k * oscillator.MassRight.u0, oscillator.SpringMiddle.k * oscillator.MassRight.v0
    u_analytical = oscillator.MassLeft.u_analytical
    v_analytical = oscillator.MassLeft.v_analytical

elif participant_name == ParticipantNames.MASS_RIGHT.value:
    n_substeps_this = args.n_substeps_right
    n_substeps_other = args.n_substeps_left

    time_stepping = args.time_stepping_right

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

num_vertices = 1  # Number of vertices

solver_process_index = 0
solver_process_size = 1

if time_stepping == TimeSteppingSchemes.GENERALIZED_ALPHA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.4, alpha_m=0.2)
elif time_stepping == TimeSteppingSchemes.NEWMARK_BETA.value:
    time_stepper = GeneralizedAlpha(stiffness=stiffness, mass=mass, alpha_f=0.0, alpha_m=0.0)
elif time_stepping == TimeSteppingSchemes.RUNGE_KUTTA_4.value:
    ode_system = np.array([
        [0,          mass], # du
        [-stiffness, 0   ], # dv
    ])
    time_stepper = RungeKutta4(ode_system=ode_system)
else:
    raise Exception(f"Invalid time stepping scheme {time_stepping}. Please use one of {[ts.value for ts in TimeSteppingSchemes]}")

print(f"time stepping scheme being used: {time_stepping}")
print(f"participant: {participant_name}")
print()
print("configured_precice_dt, my_dt, error, avg it / window, min it / window, max it / window")

dts = [0.2, 0.1, 0.05, 0.025, 0.0125]

errors = []
my_dts = []
avg_iterations = []
min_iterations = []
max_iterations = []

for dt in dts:

    configured_precice_dt = dt
    my_dt = dt/n_substeps_this
    other_dt = dt/n_substeps_other

    waveform_order = args.interpolation_degree
    render(configured_precice_dt, waveform_order, args.substeps, args.acceleration_scheme)

    configuration_file_name = f"configs/precice-config.xml"

    participant = precice.Participant(participant_name, configuration_file_name, solver_process_index, solver_process_size)

    dimensions = participant.get_mesh_dimensions(mesh_name)

    vertex = np.zeros(dimensions)
    read_data = np.zeros(num_vertices)
    write_data = oscillator.SpringMiddle.k * u0 * np.ones(num_vertices)

    vertex_ids = [participant.set_mesh_vertex(mesh_name, vertex)]

    if participant.requires_initial_data():
        participant.write_data(mesh_name, write_data_name, vertex_ids, write_data)

    participant.initialize()
    precice_dt = participant.get_max_time_step_size()
    dt = np.min([precice_dt, my_dt])

    # Initial Conditions

    a0 = (f0 - stiffness * u0) / mass
    u = u0
    v = v0
    a = a0
    t = 0

    positions = []
    velocities = []
    times = []
    n_iterations = []

    u_write = [u]
    v_write = [v]
    t_write = [t]

    iterations = 0

    while participant.is_coupling_ongoing():
        if participant.requires_writing_checkpoint():
            u_cp = u
            v_cp = v
            a_cp = a
            t_cp = t

            # store data for plotting and postprocessing
            positions += u_write
            velocities += v_write
            times += t_write

            if iterations > 0:  # don't write data, when writing checkpoint at beginning of first window
                n_iterations.append(iterations)

            iterations = 1

        t_f = time_stepper.rhs_eval_points(dt)
        f = len(t_f)*[None]

        for i in range(len(t_f)):
            read_data = participant.read_data(mesh_name, read_data_name, vertex_ids, t_f[i])
            f[i] = read_data[0]

        # do time stepping
        u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)

        write_data = [oscillator.SpringMiddle.k * u_new]

        participant.write_data(mesh_name, write_data_name, vertex_ids, write_data)

        participant.advance(dt)
        precice_dt = participant.get_max_time_step_size()
        dt = np.min([precice_dt, my_dt])

        if participant.requires_reading_checkpoint():
            u = u_cp
            v = v_cp
            a = a_cp
            t = t_cp

            # empty buffers for next window
            u_write = []
            v_write = []
            t_write = []
            iterations += 1

        else:
            u = u_new
            v = v_new
            a = a_new
            t += dt

            # write data to buffers
            u_write.append(u)
            v_write.append(v)
            t_write.append(t)

    positions += u_write
    velocities += v_write
    times += t_write

    participant.finalize()

    # print errors
    # analytic solution is only valid for this setup!
    error = np.max(abs(u_analytical(np.array(times))-np.array(positions)))
    errors.append(error)
    my_dts.append(my_dt)
    avg_iterations.append(np.average(n_iterations))
    min_iterations.append(np.min(n_iterations))
    max_iterations.append(np.max(n_iterations))
    print(f"{configured_precice_dt}, {my_dt}, {error}, {np.average(n_iterations)}, {np.min(n_iterations)}, {np.max(n_iterations)}")

df = pd.DataFrame(index=dts)
df.index.name = "dt"
df["my_dt"] = my_dts
df["error"] = errors
df["avg(iterations)"] = avg_iterations
df["min(iterations)"] = min_iterations
df["max(iterations)"] = max_iterations

interpolation_scheme = ReadWaveformSchemes.BSPLINE

# TODO: Also add acceleration scheme in metadata!

filepath = this_file.parent / f"{Cases.PRECICE3.value}_{participant_name}_{args.time_stepping_left}_{args.time_stepping_right}_{interpolation_scheme}_{args.interpolation_degree}_{args.n_substeps_left}_{args.n_substeps_right}_{args.acceleration_scheme}.csv"

df.to_csv(filepath)

add_metainfo(runner_file=this_file,
             csv_file=filepath,
             time_stepping_scheme_left=args.time_stepping_left,
             time_stepping_scheme_right=args.time_stepping_right,
             precice_version=precice.__version__,
             read_waveform_scheme=interpolation_scheme,
             read_waveform_order=args.interpolation_degree,
             acceleration_scheme=args.acceleration_scheme,
             n_substeps_left=args.n_substeps_left,
             n_substeps_right=args.n_substeps_right)