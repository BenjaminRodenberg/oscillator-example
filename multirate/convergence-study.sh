. .venv/bin/activate
# BSpline Interpolation
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -is BSpline -wd 1 -sb 1 1 -o convergence-studies/bspline/multirate_radauIIA_degree1.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -is BSpline -wd 2 -sb 2 2 -o convergence-studies/bspline/multirate_radauIIA_degree2.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -is BSpline -wd 3 -sb 4 4 -o convergence-studies/bspline/multirate_radauIIA_degree3.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -is BSpline -wd 4 -sb 8 8 -o convergence-studies/bspline/multirate_radauIIA_degree4.csv
