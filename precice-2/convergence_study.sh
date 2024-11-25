. .venv/bin/activate
# Higher-order time stepping
python3 doConvergenceStudy.py precice-config-template.xml -tss Newmark_beta Newmark_beta -o convergence-studies/higher-order/precice_2_Newmark_beta.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss generalized_alpha generalized_alpha -o convergence-studies/higher-order/precice_2_generalized_alpha.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss runge_kutta_4 runge_kutta_4 -o convergence-studies/higher-order/precice_2_runge_kutta_4.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -o convergence-studies/higher-order/precice_2_radauIIA.csv
# Subcycling
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss Newmark_beta Newmark_beta -o convergence-studies/subcycling/precice_2_Newmark_beta.csv
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss generalized_alpha generalized_alpha -o convergence-studies/subcycling/precice_2_generalized_alpha.csv
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss runge_kutta_4 runge_kutta_4 -o convergence-studies/subcycling/precice_2_runge_kutta_4.csv
python3 doConvergenceStudy.py precice-config-template.xml -s 5 -sb 1 1 -sf 2 1 -w 1 -tss radauIIA radauIIA -o convergence-studies/subcycling/precice_2_radauIIA.csv