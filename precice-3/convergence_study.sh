. .venv/bin/activate
# Higher-order time stepping
python3 doConvergenceStudy.py precice-config-template.xml -tss Newmark_beta Newmark_beta -o convergence-studies/higher-order/precice_3_Newmark_beta.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss generalized_alpha generalized_alpha -o convergence-studies/higher-order/precice_3_generalized_alpha.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss runge_kutta_4 runge_kutta_4 -o convergence-studies/higher-order/precice_3_runge_kutta_4.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -o convergence-studies/higher-order/precice_3_radauIIA.csv
# Subcycling
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss Newmark_beta Newmark_beta -o convergence-studies/subcycling/precice_3_Newmark_beta.csv
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss generalized_alpha generalized_alpha -o convergence-studies/subcycling/precice_3_generalized_alpha.csv
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss runge_kutta_4 runge_kutta_4 -o convergence-studies/subcycling/precice_3_runge_kutta_4.csv
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss radauIIA radauIIA -o convergence-studies/subcycling/precice_3_radauIIA.csv
# Runge Kutta 4
python3 doConvergenceStudy.py precice-config-template.xml -tss runge_kutta_4 runge_kutta_4 -wd 1 -dt 0.2 -sb 5 5 -o convergence-studies/runge-kutta-4/precice_3_runge_kutta_4_degree1.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss runge_kutta_4 runge_kutta_4 -wd 2 -dt 0.2 -sb 5 5 -o convergence-studies/runge-kutta-4/precice_3_runge_kutta_4_degree2.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss runge_kutta_4 runge_kutta_4 -wd 3 -dt 0.2 -sb 5 5 -o convergence-studies/runge-kutta-4/precice_3_runge_kutta_4_degree3.csv
# RadauIIA
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -wd 1 -dt 0.2 -sb 5 5 -o convergence-studies/radauIIA/precice_3_radauIIA_degree1.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -wd 2 -dt 0.2 -sb 5 5 -o convergence-studies/radauIIA/precice_3_radauIIA_degree2.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -wd 3 -dt 0.2 -sb 5 5 -o convergence-studies/radauIIA/precice_3_radauIIA_degree3.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -wd 4 -dt 0.2 -sb 5 5 -o convergence-studies/radauIIA/precice_3_radauIIA_degree4.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -wd 5 -dt 0.2 -sb 5 5 -o convergence-studies/radauIIA/precice_3_radauIIA_degree5.csv
