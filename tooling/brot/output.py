import git
import pathlib

from .enums import Cases, TimeSteppingSchemes

def add_metainfo(runner_file, csv_file, time_stepping_scheme, precice_version='not available'):
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
        f"# experiment: {Cases.MONOLITHIC.value}\n"
        f"# time stepping scheme: {time_stepping_scheme}\n"
    )

    with open(csv_file, "r") as original:
        data = original.read()
    with open(csv_file, "w") as modified:
        modified.write(metainfo + data)