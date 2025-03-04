from pathlib import Path
import uuid
import argparse

from prepesthel.participant import Participant, Participants
from prepesthel.runner import run, postproc
from prepesthel.io import Results, Executors


if __name__ == "__main__":
    n_supported_participants = 2

    parser = argparse.ArgumentParser(description="Solving oscillator example.")
    parser.add_argument(
        "template_path",
        help="template for the preCICE configuration file",
        type=str)
    parser.add_argument(
        "--silent",
        help="Deactivates result output to command line",
        action='store_true')
    parser.add_argument(
        "--executor",
        help="Define type of executor",
        type=str,
        choices=[e.value for e in Executors],
        default=Executors.LOCAL.value)
    parser.add_argument(
        "-T",
        "--max-time",
        help="Max simulation time",
        type=float,
        default=1.0)
    parser.add_argument(
        "-dt",
        "--base-time-window-size",
        help="Base time window size",
        type=float,
        default=0.04)
    parser.add_argument(
        "-w",
        "--time-window-refinements",
        help="Number of refinements by factor 2 for the time window size",
        type=int,
        default=5)
    parser.add_argument(
        "-sb",
        "--base-time-step-refinement",
        help="Base factor for time window size / time step size",
        type=int,
        nargs=n_supported_participants,
        default=n_supported_participants*[1])
    parser.add_argument(
        "-s",
        "--time-step-refinements",
        help="Number of refinements by given factor for the time step size of each participant ( >1 will result in subcycling)",
        type=int,
        default=1)
    parser.add_argument(
        "-sf",
        "--time-step-refinement-factor",
        help="Factor of time step refinements for each participant (use 1, if you want to use a fixed time step / time window relationship for one participant while refining the time steps for the other participant)",
        type=int,
        nargs=n_supported_participants,
        default=n_supported_participants*[2])
    ## add solver specific arguments below, if needed
    parser.add_argument(
        "-tss",
        "--time-stepping-scheme",
        help="Define time stepping scheme used by each solver",
        type=str,
        nargs=n_supported_participants,
        default=n_supported_participants*["Newmark_beta"])
    parser.add_argument(
        "-wd",
        "--waveform-degree",
        help="Waveform degree being used",
        type=int,
        default=1)
    parser.add_argument(
        "-o",
        "--out-filename",
        help="Provide a file name. If no file name is provided a UUID will be generated as name. Abort if file already exists.",
        type=str,
    )
    args = parser.parse_args()

    # Define values that will be inserted into precice-config-template.xml here
    precice_config_params = {
        'time_window_size': None,  # will be defined later
        'max_time': args.max_time,
        'waveform_degree': args.waveform_degree,
        'substeps': True,
    }

    root_folder = Path(__file__).parent.absolute()

    # Define how participants will be executed here
    participants: Participants = {
        "Mass-Left": Participant("Mass-Left", root_folder, [".venv/bin/python3", "oscillator.py"], ["Mass-Left"], {'--time-stepping': args.time_stepping_scheme[0], '--n-substeps': None}),
        "Mass-Right": Participant("Mass-Right", root_folder, [".venv/bin/python3", "oscillator.py"], ["Mass-Right"], {'--time-stepping': args.time_stepping_scheme[1], '--n-substeps': None}),        
    }

    if len(participants) != n_supported_participants:
        raise Exception(f"Currently only supports coupling of {n_supported_participants} participants")

    results_file_path = root_folder
    if args.out_filename:  # use file name given by user
        results_file_path = results_file_path / args.out_filename
    else:  # no file name is given. Create UUID for file name
        results_file_path = results_file_path / "convergence-studies" / f"{uuid.uuid4()}.csv"

    results = Results(results_file_path)

    for dt in [args.base_time_window_size * 0.5**i for i in range(args.time_window_refinements)]:
        for refinement in range(args.time_step_refinements):
            precice_config_params['time_window_size'] = dt
            i = 0
            for p in participants.values():
                p.kwargs['--n-substeps'] = args.base_time_step_refinement[i]*args.time_step_refinement_factor[i]**refinement
                i += 1

            run(participants, args.template_path, precice_config_params)
            summary = postproc(participants, precice_config_params, silent=args.silent)

            results.append(summary)
            results.output_preliminary(silent=args.silent)
    
    results.output_final(participants, args, precice_config_params, silent=args.silent, executor=args.executor)