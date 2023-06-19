from jinja2 import Environment, select_autoescape, FileSystemLoader
import os, stat

env = Environment(
    loader=FileSystemLoader('.'),
    autoescape=select_autoescape(['xml', 'json'])
)

precice_config_template = env.get_template('precice-config-template.xml')
precice_config_name = "precice-config.xml"

n_left = n_right = 4

with open(os.path.join( ".", precice_config_name), "w") as file:
    file.write(precice_config_template.render(time_window_size=0.01,
                                              forces_left=[f"Force-Left-{i}" for i in range(1, n_left+1)],
                                              forces_right=[f"Force-Right-{i}" for i in range(1, n_right+1)]))