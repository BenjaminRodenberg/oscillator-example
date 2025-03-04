from jinja2 import Environment, select_autoescape, FileSystemLoader
import os
from pathlib import Path
from enum import Enum


class AccelerationSchemes(Enum):
    NONE = "None"
    CONSTANT = "Constant"
    REDUCED_QUASI_NEWTON = "rQN"
    FULL_QUASI_NEWTON = "QN"


def render(dt, n_left, n_right, acceleration_scheme):
    base_path = Path(__file__).parent.absolute()

    env = Environment(
        loader=FileSystemLoader(base_path),
        autoescape=select_autoescape(['xml'])
    )

    if acceleration_scheme == AccelerationSchemes.NONE.value:
        precice_config_template = env.get_template('precice-config-template-NoAcceleration.xml')
    elif acceleration_scheme == AccelerationSchemes.CONSTANT.value:
        precice_config_template = env.get_template('precice-config-template-UR.xml')
    elif acceleration_scheme == AccelerationSchemes.REDUCED_QUASI_NEWTON.value:
        precice_config_template = env.get_template('precice-config-template-rQN.xml')
    elif acceleration_scheme == AccelerationSchemes.FULL_QUASI_NEWTON.value:
        precice_config_template = env.get_template('precice-config-template-QN.xml')
    else:
        raise Exception(f"Unknown acceleration_scheme {acceleration_scheme}")

    precice_config_name = base_path / "precice-config.xml"

    with open(os.path.join(".", precice_config_name), "w") as file:
        file.write(precice_config_template.render(time_window_size=dt,
                                                  displacements_left=[f"Displacement-Left-{i}" for i in range(1, n_left + 1)],
                                                  displacements_right=[f"Displacement-Right-{i}" for i in range(1, n_right + 1)],
                                                  convergence_limit=10e-6))


if __name__ == "__main__":
    render(dt=0.01, n_left=4, n_right=4, acceleration_scheme=AccelerationSchemes.REDUCED_QUASI_NEWTON.value)
