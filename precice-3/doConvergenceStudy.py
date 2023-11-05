from jinja2 import Environment, select_autoescape, FileSystemLoader
import pandas as pd
from pathlib import Path
import subprocess
import datetime
import os
import uuid
import argparse
from enum import Enum


def render(precice_config_params):
    base_path = Path(__file__).parent.absolute()

    env = Environment(
        loader=FileSystemLoader(base_path),
        autoescape=select_autoescape(['xml'])
    )

    precice_config_template = env.get_template('precice-config-template.xml')

    precice_config_name = base_path / "precice-config.xml"

    with open(precice_config_name, "w") as file:
        file.write(precice_config_template.render(precice_config_params))


def do_run(precice_config_params, participants):
    render(precice_config_params)
    print(f"{datetime.datetime.now()}: Start run with parameters {precice_config_params}")
    print("Running...")

    for participant in participants:
        participant['logfile'] = f"stdout-{participant['name']}.log"

    for participant in participants:
        with open(participant["folder"] / participant['logfile'], "w") as outfile:
            cmd = ["python3", participant["folder"] / participant["exec"]] + participant["params"] + [f"{keyword}={value}" for keyword, value in participant['kwargs'].items()]
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
    time_window_size = precice_config_params['time_window_size']
    summary = {"time window size": time_window_size}
    for participant in participants:
        df = pd.read_csv(participant["folder"] / f"errors-{participant['name']}.csv", comment="#")
        summary[f"time step size {participant['name']}"] = time_window_size / participant['kwargs']['--n-substeps']
        summary[f"error {participant['name']}"] = df["errors"].abs().max()
    print("Done.")

    return summary


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Solving heat equation for simple or complex interface case")
    parser.add_argument(
        "-dt",
        "--base-time-window-size",
        help="Base time window / time step size",
        type=float,
        default=0.04)
    parser.add_argument(
        "-w",
        "--time-window-refinements",
        help="Number of refinements by factor 2 for the time window size",
        type=int,
        default=5)
    parser.add_argument(
        "-s",
        "--time-step-refinements",
        help="Number of refinements by factor 2 for the time step size ( >1 will result in subcycling)",
        type=int,
        default=1)
    args = parser.parse_args()

    df = pd.DataFrame()

    # Define values that will be inserted into precice-config-template.xml here
    precice_config_params = {
        'time_window_size': '',  # will be defined later
    }

    # Define how participants will be executed here
    participants = [
        {
            "name": "Mass-Left",  # identifier of this participant
            "folder": Path(__file__).parent.absolute(),  # root folder of this participant
            "exec": "oscillator.py",  # name of the python script, e.g. script.py
            "params": ["Mass-Left"],  # list of positional arguments that will be used. Results in python3 script.py param1 ...
            "kwargs": {  # dict with keyword arguments that will be used. Results in python3 script.py param1 ... k1=v1 k2=v2 ...
                '--time-stepping': 'radauIIA',
                '--n-substeps': '',  # will be defined later
            },
        },
        {
            "name": "Mass-Right",
            "folder": Path(__file__).parent.absolute(),
            "exec": "oscillator.py",
            "params": ["Mass-Right"],
            "kwargs": {
                '--time-stepping': 'radauIIA',
                '--n-substeps': '',  # will be defined later
            },
        },
    ]

    summary_file = Path("convergence-studies") / f"{uuid.uuid4()}.csv"

    for dt in [args.base_time_window_size * 0.5**i for i in range(args.time_window_refinements)]:
        for n in [2**i for i in range(args.time_step_refinements)]:

            precice_config_params['time_window_size'] = dt
            for p in participants:
                p['kwargs']['--n-substeps'] = n

            summary = do_run(precice_config_params, participants)
            df = pd.concat([df, pd.DataFrame(summary, index=[0])], ignore_index=True)

            print(f"Write preliminary output to {summary_file}")
            df.to_csv(summary_file)

            term_size = os.get_terminal_size()
            print('-' * term_size.columns)
            print(df)
            print('-' * term_size.columns)

    df = df.set_index(["time window size"] + [f"time step size {p['name']}" for p in participants])
    print(f"Write final output to {summary_file}")

    import git

    repo = git.Repo(__file__, search_parent_directories=True)
    chash = str(repo.head.commit)[:7]
    if repo.is_dirty():
        chash += "-dirty"

    metadata = {
        "git repository": repo.remotes.origin.url,
        "git commit": chash,
        "args": args,
        "precice_config_params": precice_config_params,
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
