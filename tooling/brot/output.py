import git
import pathlib

from .enums import Cases, TimeSteppingSchemes

def add_metainfo(runner_file, csv_file, time_stepping_scheme_left='not available', time_stepping_scheme_right='not available',  precice_version='not available', read_waveform_scheme='not available', read_waveform_order='not available', acceleration_scheme='not available', n_substeps_left='not available', n_substeps_right='not available'):
    repo_base = runner_file.parent / ".."

    repo = git.Repo(repo_base)
    chash = str(repo.head.commit)[:7]
    if repo.is_dirty():
        chash += "-dirty"

    repourl = repo.remotes.origin.url

    metainfo = (
        f"# git repo: {repourl}\n"
        f"# git commit: {chash}\n"
        f"# precice version: {precice_version}\n"
        f"# executable: {pathlib.Path(runner_file.resolve()).relative_to(repo_base.resolve())}\n"
        f"# time stepping scheme left: {time_stepping_scheme_left}\n"
        f"# time stepping scheme right: {time_stepping_scheme_right}\n"
        f"# read waveform scheme: {read_waveform_scheme}\n"
        f"# read waveform order: {read_waveform_order}\n"
        f"# acceleration scheme: {acceleration_scheme}\n"
        f"# n substeps left: {n_substeps_left}\n"
        f"# n substeps right: {n_substeps_right}\n"
    )

    with open(csv_file, "r") as original:
        data = original.read()
    with open(csv_file, "w") as modified:
        modified.write(metainfo + data)