`mass-spring-substeps.py` solves mass-spring-system in partitioned fashion.

# Multirate case

## Config file

Is created at runtime by executing `./substeps_configs` with the subcycling configuration for the two participants. This fills a jinja2 template and creates an instance of the configuration file.

## Running

```
python3 mass-spring.py Mass-Left
```

and

```
python3 mass-spring.py Mass-Right
```

## Experiments

### Higher-degree Lagrange interpolation

Use the options `-is=Lagrange -tsl=runge_kutta_4 -tsr=runge_kutta_4`. Increase the number of substeps `-nl=1 -nr=1` ... `-nl=8 -nr=8`. With increasing number of substeps the degree of the Lagrange interpolation will increase. This allows 4-th order convergence for `-nl=3 -nr=3` and larger.

### Higher-degree interpolation with B-splines

Use the options `-is=B_spline -tsl=runge_kutta_4 -tsr=runge_kutta_4`. Increase the number of substeps `-nl=1 -nr=1` ... `-nl=8 -nr=8`. Pick B-spline degree with `-p=1`. This allows 4-th order convergence for `-nl=2 -nr=2 -p=2` and larger.

### Use multirate to compensate lower order time stepping scheme