from pathlib import Path
import uuid
import argparse

from prepesthel.participant import Participant, Participants
from prepesthel.runner import run, postproc
from prepesthel.io import Results, Executors


if __name__ == "__main__":
    n_supported_participants = 1  # convergence study just for a monolithic simulation

    parser = argparse.ArgumentParser(description="Solving oscillator example.")
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
        "-dt",
        "--base-time-step-size",
        help="Base time step size",
        type=float,
        default=0.04)
    parser.add_argument(
        "-w",
        "--time-step-refinements",
        help="Number of refinements by factor 2 for the time step size",
        type=int,
        default=5)
    # add solver specific arguments below, if needed
    parser.add_argument(
        "-tss",
        "--time-stepping-scheme",
        help="Define time stepping scheme used by each solver",
        type=str,
        nargs=n_supported_participants,
        default=n_supported_participants * ["Newmark_beta"])
    parser.add_argument(
        "-o",
        "--out-filename",
        help="Provide a file name. If no file name is provided a UUID will be generated as name. Abort if file already exists.",
        type=str,
    )
    args = parser.parse_args()

    root_folder = Path()

    # Define how participants will be executed here
    participants: Participants = {
        "Monolithic": Participant("Monolithic", root_folder, [".venv/bin/python3", "oscillator.py"], [], {'--time-stepping': args.time_stepping_scheme[0], '--time-step-size': None})
    }

    if len(participants) != n_supported_participants:
        raise Exception(f"Currently only supports coupling of {n_supported_participants} participants")

    results_file_path = root_folder
    if args.out_filename:  # use file name given by user
        results_file_path = results_file_path / args.out_filename
    else:  # no file name is given. Create UUID for file name
        results_file_path = results_file_path / "convergence-studies" / f"{uuid.uuid4()}.csv"

    results = Results(results_file_path)

    for dt in [args.base_time_step_size * 0.5**i for i in range(args.time_step_refinements)]:
        for p in participants.values():
            p.kwargs['--time-step-size'] = dt

        run(participants)
        summary = postproc(participants, silent=args.silent)

        results.append(summary)
        results.output_preliminary(silent=args.silent)

    results.output_final(participants, args, silent=args.silent, executor=args.executor)
