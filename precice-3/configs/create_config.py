from jinja2 import Environment, select_autoescape, FileSystemLoader
import os
from pathlib import Path


def render(dt, waveform_order):
    base_path = Path(__file__).parent.absolute()

    env = Environment(
        loader=FileSystemLoader(base_path),
        autoescape=select_autoescape(['xml'])
    )

    precice_config_template = env.get_template('precice-config-template.xml')
    precice_config_name = base_path / "precice-config.xml"

    with open(os.path.join( ".", precice_config_name), "w") as file:
        file.write(precice_config_template.render(time_window_size=dt,
                                                  waveform_order=waveform_order))

if __name__ == "__main__":
    render(dt=0.01, waveform_order=1)
