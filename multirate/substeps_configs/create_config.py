from jinja2 import Environment, select_autoescape, FileSystemLoader
import os
from pathlib import Path


def render(dt, n_left, n_right):
    base_path = Path(__file__).parent.absolute()

    env = Environment(
        loader=FileSystemLoader(base_path),
        autoescape=select_autoescape(['xml'])
    )

    precice_config_template = env.get_template('precice-config-template.xml')
    precice_config_name = base_path / "precice-config.xml"

    with open(os.path.join( ".", precice_config_name), "w") as file:
        file.write(precice_config_template.render(time_window_size=dt,
                                                forces_left=[f"Force-Left-{i}" for i in range(1, n_left+1)],
                                                forces_right=[f"Force-Right-{i}" for i in range(1, n_right+1)]))

if __name__ == "__main__":
    render(dt=0.01, n_left=4, n_right=4)
