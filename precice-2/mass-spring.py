from __future__ import division

import argparse
import numpy as np
import pandas as pd
import precice
import pathlib

from brot.enums import Cases, TimeSteppingSchemes, MultirateMode, ParticipantNames, DataNames, MeshNames
from brot.output import add_metainfo
from brot.timesteppers import GeneralizedAlpha, RungeKutta4
import brot.oscillator as oscillator

this_file = pathlib.Path(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("participantName", help="Name of the solver.", type=str)
parser.add_argument("-ts", "--time-stepping", help="Time stepping scheme being used.", type=str, default=TimeSteppingSchemes.NEWMARK_BETA.value)
parser.add_argument("-mr", "--multirate", help="Pick one of 3 modes for multirate", type=str, default=MultirateMode.NONE.value)

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

dts = [0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625]

errors = []
my_dts = []

for dt in dts:
    if args.multirate == MultirateMode.NONE.value:
        # use same dt for both solvers and preCICE
        my_dt = dt
        configured_precice_dt = dt
    elif args.multirate == MultirateMode.SUBCYCLING.value:
        # use fixed dt for preCICE
        configured_precice_dt = np.max(dts)
        my_dt = dt / 4
    elif args.multirate == MultirateMode.FOUR_SUBSTEPS.value:
        configured_precice_dt = dt
        # always use four substeps
        my_dt = dt / 4
    elif args.multirate == MultirateMode.MULTIRATE.value:
        # use fixed dt for preCICE
        configured_precice_dt = 0.04
        if participant_name == ParticipantNames.MASS_LEFT.value:
            # use fixed dt for left participant
            my_dt = 0.01
        elif participant_name == ParticipantNames.MASS_RIGHT.value:
            my_dt = dt
    else:
        raise Exception(f"wrong multirate mode: {args.multirate}. Please use one of {[m.value for m in MultirateMode]}.")

    # use same dt for both solvers and preCICE
    configuration_file_name = f"configs/precice-config-{configured_precice_dt}.xml"

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

    while interface.is_coupling_ongoing():
        if interface.is_action_required(
                precice.action_write_iteration_checkpoint()):
            u_cp = u
            v_cp = v
            a_cp = a
            t_cp = t
            f_start = f_end  # force at the beginning of the window
            interface.mark_action_fulfilled(
                precice.action_write_iteration_checkpoint())

            # store data for plotting and postprocessing
            positions += u_write
            velocities += v_write
            times += t_write

        read_data = interface.read_scalar_data(read_data_id, vertex_id)
        f = read_data

        # do time stepping
        u_new, v_new, a_new = time_stepper.do_step(u, v, a, f, dt)

        write_data = oscillator.SpringMiddle.k * u_new

        interface.write_scalar_data(
            write_data_id, vertex_id, write_data)

        precice_dt = interface.advance(dt)
        dt = np.min([precice_dt, my_dt])

        if interface.is_action_required(
                precice.action_read_iteration_checkpoint()):
            u = u_cp
            v = v_cp
            a = a_cp
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
            a = a_new
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
    my_dts.append(my_dt)
    print(f"{configured_precice_dt}, {my_dt}, {error}")

df = pd.DataFrame(index=dts)
df.index.name = "dt"
df["my_dt"] = my_dts
df["error"] = errors

time_stepping_scheme = args.time_stepping
filepath = this_file.parent / f"{Cases.PRECICE2.value}_{participant_name}_{time_stepping_scheme}_{args.multirate}.csv"
df.to_csv(filepath)

add_metainfo(this_file, filepath, time_stepping_scheme, precice.__version__)
