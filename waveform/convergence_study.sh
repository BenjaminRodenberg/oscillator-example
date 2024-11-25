. .venv/bin/activate
# Lagrange Interpolation
python3 doConvergenceStudy.py precice-config-template.xml -tss Newmark_beta Newmark_beta -o convergence-studies/lagrange/waveform_Newmark_beta.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss generalized_alpha generalized_alpha -o convergence-studies/lagrange/waveform_generalized_alpha.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss runge_kutta_4 runge_kutta_4 -o convergence-studies/lagrange/waveform_runge_kutta_4.csv
python3 doConvergenceStudy.py precice-config-template.xml -tss radauIIA radauIIA -o convergence-studies/lagrange/waveform_radauIIA.csv
# Hermite Interpolation
python3 doConvergenceStudy.py precice-config-hermite-template.xml -is Hermite -tss Newmark_beta Newmark_beta -o convergence-studies/hermite/waveform_Newmark_beta.csv
python3 doConvergenceStudy.py precice-config-hermite-template.xml -is Hermite -tss generalized_alpha generalized_alpha -o convergence-studies/hermite/waveform_generalized_alpha.csv
python3 doConvergenceStudy.py precice-config-hermite-template.xml -is Hermite -tss runge_kutta_4 runge_kutta_4 -o convergence-studies/hermite/waveform_runge_kutta_4.csv
python3 doConvergenceStudy.py precice-config-hermite-template.xml -is Hermite -tss radauIIA radauIIA -o convergence-studies/hermite/waveform_radauIIA.csv