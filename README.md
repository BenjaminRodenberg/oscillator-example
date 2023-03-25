# oscillator-example

This example was developed in the scope of the dissertation of Benjamin Rodenberg[^1]. A modified version of the example was already published in the conference proceedings of WCCM-APCOM 2022[^2].

The examples provided in this repository illustrate the whole development from a monolithic simulation of the two-mass-spring system (oscillator example) over a partitioned approach using preCICE version 2.5 to an approach using the newly developed higher-order coupling scheme implemented in preCICE 3.

The following examples are included here:

1) `monolithic`: A monolithic solution of the oscillator example.
2) `precice-2`: A partitioned solution of the oscillator using preCICE 2.5[^3].
3) `waveform`: A partitioned solution of the oscillator using preCICE 2.5[^3]. Implementing waveform iteration ideas from[^1] in the adapter code.
4) `multirate`: A partitioned solution of the oscillator using preCICE 2.5[^3]. Implementing waveform iteration ideas from[^1] in the adapter code. Additionally illustrating improvements for subcycling.
5) `acceleration`: A partitioned solution of the oscillator using preCICE 2.5[^3]. Implementing waveform ideas from [^1] in the adapter code. Additionally illustrating how to implement underrelaxation and quasi-Newton waveform iteration following the approach from [^5]
6) `precice-3`: A partitioned solution of the oscillator using preCICE 3.x[^4]. The native API from preCICE 3.x is used to apply waveform iteration.

## Code Reproducibility

It is generally recommended to use the appropriate version of the preCICE virtual machine[^3,^4] for run all experiments. See https://precice.org/installation-vm.html for details.

[^1]: B. Rodenberg, **in preparation**
[^2]: V. Schüller, B. Rodenberg, B. Uekermann and H. Bungartz, A Simple Test Case for Convergence Order in Time and Energy Conservation of Black-Box Coupling Schemes, in: WCCM-APCOM2022. URL https://www.scipedia.com/public/Rodenberg_2022a
[^3]: preCICE VM with distribution version v202211.0.0 https://github.com/precice/vm/releases/tag/v202211.0.0
[^4]: preCICE VM with distribution version vYYYYMM.0.0 **not yet available**
[^5]: Rüth, B, Uekermann, B, Mehl, M, Birken, P, Monge, A, Bungartz, H-J. Quasi-Newton waveform iteration for partitioned surface-coupled multiphysics applications. Int J Numer Methods Eng. 2021; 122: 5236– 5257. https://doi.org/10.1002/nme.6443
