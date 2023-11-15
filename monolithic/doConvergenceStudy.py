from jinja2 import Environment, select_autoescape, FileSystemLoader
import pandas as pd
from pathlib import Path
import subprocess
import datetime
import os
import uuid
import argparse
import sys


def do_run(participants):
    print(f"{datetime.datetime.now()}: Running ...")

    for participant in participants:
        participant['logfile'] = f"stdout-{participant['name']}.log"

    for participant in participants:
        with open(participant["folder"] / participant['logfile'], "w") as outfile:
            cmd = participant["exec"] + participant["params"] + [f"{keyword}={value}" for keyword, value in participant['kwargs'].items()]
            print(cmd)
            p = subprocess.Popen(cmd,
                                 cwd=participant["folder"],
                                 stdout=outfile)
            participant["proc"] = p

    for participant in participants:
        participant["proc"].wait()

    for participant in participants:
        if participant["proc"].returncode != 0:
            raise Exception(f"Experiment failed for participant {participant['name']}. See logs {[p['logfile'] for p in participants]}.")

    print("Done.")
    print("Postprocessing...")
    summary = {}
    for participant in participants:
        df = pd.read_csv(participant["folder"] / f"errors-{participant['name']}.csv", comment="#")
        if abs(df.times.diff().var() / df.times.diff().mean()) > 10e-10:
            term_size = os.get_terminal_size()
            print('-' * term_size.columns)
            print("WARNING: times vary stronger than expected. Note that adaptive time stepping is not supported.")
            print(df)
            print('-' * term_size.columns)
        summary[f"time step size {participant['name']}"] = df.times.diff().mean()
        summary[f"error1 {participant['name']}"] = df.errors1.abs().max()
        summary[f"error2 {participant['name']}"] = df.errors1.abs().max()
    print("Done.")

    return summary


if __name__ == "__main__":
    n_supported_participants = 1  # convergence study just for a monolithic simulation

    parser = argparse.ArgumentParser(description="Solving oscillator example.")
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
    ## add solver specific arguments below, if needed
    parser.add_argument(
        "-tss",
        "--time-stepping-scheme",
        help="Define time stepping scheme used by each solver",
        type=str,
        nargs=n_supported_participants,
        default=n_supported_participants*["Newmark_beta"])
    args = parser.parse_args()

    df = pd.DataFrame()

    root_folder = Path(__file__).parent.absolute()

    # Define how participants will be executed here
    participants = [
        {
            "name": "Oscillator",  # identifier of this participant
            "folder": root_folder,  # root folder of this participant
            "exec": ["python3", "oscillator.py"],  # how to execute the participant, e.g. python3 script.py
            "params": [],  # list of positional arguments that will be used. Results in python3 script.py param1 ...
            "kwargs": {  # dict with keyword arguments that will be used. Results in python3 script.py param1 ... k1=v1 k2=v2 ...
                '--time-stepping': args.time_stepping_scheme[0],
                '--time-step-size': None,  # will be defined later
            },
        }
    ]

    if len(participants) != n_supported_participants:
        raise Exception(f"Currently only supports coupling of {n_supported_participants} participants")

    summary_file = root_folder / "convergence-studies" / f"{uuid.uuid4()}.csv"

    for dt in [args.base_time_step_size * 0.5**i for i in range(args.time_step_refinements)]:
        i = 0
        for p in participants:
            p['kwargs']['--time-step-size'] = dt
            i += 1

        summary = do_run(participants)
        df = pd.concat([df, pd.DataFrame(summary, index=[0])], ignore_index=True)

        print(f"Write preliminary output to {summary_file}")
        df.to_csv(summary_file)

        term_size = os.get_terminal_size()
        print('-' * term_size.columns)
        print(df)
        print('-' * term_size.columns)

    df = df.set_index([f"time step size {p['name']}" for p in participants])
    print(f"Write final output to {summary_file}")

    import git

    repo = git.Repo(__file__, search_parent_directories=True)
    chash = str(repo.head.commit)[:7]
    if repo.is_dirty():
        chash += "-dirty"

    metadata = {
        "git repository": repo.remotes.origin.url,
        "git commit": chash,
        "run cmd": "python3 " + " ".join(sys.argv),
        "args": args,
        "participants": participants,
    }

    summary_file.unlink()

    with open(summary_file, 'a') as f:
        for key, value in metadata.items():
            f.write(f"# {key}:{value}\n")
        df.to_csv(f)

    print('-' * term_size.columns)
    for key, value in metadata.items():
        print(f"{key}:{value}")
    print()
    print(df)
    print('-' * term_size.columns)
