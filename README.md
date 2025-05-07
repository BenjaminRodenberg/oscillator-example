# Oscillator example for higher-order and multirate partitioned time stepping with preCICE

This example was developed in the scope of the dissertation of Benjamin Rodenberg[^1]. A modified version of the example was already published in the conference proceedings of WCCM-APCOM 2022[^2].

The examples provided in this repository illustrate the whole development from a monolithic simulation of the two-mass-spring system (oscillator example) over a partitioned approach using preCICE version 2.5 to an approach using the newly developed higher-order coupling scheme implemented in preCICE 3.

The following examples are included here:

1) `monolithic`: A monolithic solution of the oscillator example.
2) `precice-2`: A partitioned solution of the oscillator using preCICE 2.5.0[^3].
3) `waveform`: A partitioned solution of the oscillator using preCICE 2.5.0[^3]. Implementing waveform iteration ideas from[^1] in the adapter code.
4) `multirate`: A partitioned solution of the oscillator using preCICE 2.5.0[^3]. Implementing waveform iteration ideas from[^1] in the adapter code. Additionally illustrating improvements for subcycling.
5) `acceleration`: A partitioned solution of the oscillator using preCICE 2.5.0[^3]. Implementing waveform ideas from [^1] in the adapter code. Additionally illustrating how to implement underrelaxation and quasi-Newton waveform iteration following the approach from [^5].
6) `precice-3`: A partitioned solution of the oscillator using preCICE 3.1.1[^4]. The native API from preCICE 3.x is used to apply waveform iteration.

## Running the code

Each of the provided cases offers the script `doConvergenceStudy.py` to automatically run convergence studies for each case with a given set of parameters. Refer to `python3 doConvergenceStudy.py --help` for available parameters.

If you want to individually run the case, please run `python3 oscillator.py`. For the monolithic case (1) it is sufficient to run this command in a single terminal. If you want to run one of the partitioned cases (2-6) please open two terminals and run `python3 oscillator.py Mass-Left` and `python3 oscillator.py Mass-Right`. Again, refer to `python3 oscillator.py --help` for a list of available parameters.

### Dependencies

Please make sure to install the package `brot` included in this repository. Refer to `tooling/README.md` for further information.
If you want to perform convergence studies, please install the package `prepesthel`. Refer to https://github.com/BenjaminRodenberg/precice-performance-study-helpers for further information.

## Code Reproducibility

It is generally recommended to use the appropriate version of the preCICE virtual machine[^3][^4] for run all experiments. See https://precice.org/installation-vm.html for details.

You can also refer to the folder `.github/workflows` for CI jobs that run all experiments provided in this repository.

[^1]: Rodenberg, Benjamin. *Flexible and robust time stepping for partitioned multiphysics*. Technical University of Munich, 2025.
[^2]: Schüller, Valentina; Rodenberg, Benjamin; Uekermann, Benjamin; Bungartz, Hans-Joachim. *A Simple Test Case for Convergence Order in Time and Energy Conservation of Black-Box Coupling Schemes*. WCCM-APCOM2022, 2022. Available at: https://www.scipedia.com/public/Rodenberg_2022a.
[^3]: Chourdakis, Gerasimos; Davis, Kyle; Desai, Ishaan; Rodenberg, Benjamin; Schneider, David; Simonis, Frédéric; Uekermann, Benjamin; Ariguib, Boshra; Cardiff, Philip; Jaust, Alexander; Kharitenko, Pavel; Klöfkorn, Robert; Kotarsky, Niklas; Martin, Boris; Scheurer, Erik; Schüller, Valentina; van Zwieten, Gertjan; Yurt, Kursat. *preCICE Distribution Version v2211.0*. DaRUS, 2023, V1. https://doi.org/10.18419/darus-3576.
[^4]: Chen, Jun; Chourdakis, Gerasimos; Desai, Ishaan; Homs-Pons, Carme; Rodenberg, Benjamin; Schneider, David; Simonis, Frédéric; Uekermann, Benjamin; Davis, Kyle; Jaust, Alexander; Kelm, Mathis; Kotarsky, Niklas; Kschidock, Helena; Mishra, Durganshu; Mühlhäußer, Markus; Schrader, Timo Pierre; Schulte, Miriam; Seitz, Valentin; Signorelli, Joseph; van Zwieten, Gertjan; Vinnitchenko, Niklas; Vladimirova, Tina; Willeke, Leonard; Zonta, Elia. *preCICE Distribution Version v2404.0*. DaRUS, 2024, V1. https://doi.org/10.18419/darus-4167.
[^5]: Rüth, Benjamin; Uekermann, Benjamin; Mehl, Mehl; Birken, Philipp; Monge, Azahar; Bungartz, Hans-Joachim *Quasi-Newton waveform iteration for partitioned surface-coupled multiphysics applications*. International Journal for Numerical Methods in Engineering, 2021, 122: 5236–5257. https://doi.org/10.1002/nme.6443.

