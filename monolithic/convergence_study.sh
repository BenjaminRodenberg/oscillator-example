. .venv/bin/activate
python3 doConvergenceStudy.py -tss Newmark_beta -w 10 -o convergence-studies/monolithic_Newmark_beta.csv
python3 doConvergenceStudy.py -tss generalized_alpha -w 10 -o convergence-studies/monolithic_generalized_alpha.csv
python3 doConvergenceStudy.py -tss runge_kutta_4 -w 10 -o convergence-studies/monolithic_runge_kutta_4.csv
python3 doConvergenceStudy.py -tss radauIIA -w 10 -o convergence-studies/monolithic_radauIIA.csv