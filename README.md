# oscillator-example

This example was developed in the scope of the dissertation of Benjamin Rodenberg[^1]. A modified version of the example was already published in the conference proceedings of WCCM-APCOM 2022[^2].

The examples provided in this repository illustrate the whole development from a monolithic simulation of the two-mass-spring system (oscillator example) over a partitioned approach using preCICE 2 to an approach using the newly developed higher-order coupling scheme implemented in preCICE 3.

The following examples are included here:

1) `monolithic`: A monolithic solution of the oscillator example.
2) `precice-2`: A partitioned solution of the oscillator using preCICE 2.
3) `waveform`: A particioned solution of the oscillator using preCICE 2. Implementing waveform iteration ideas from[^1] in the adapter code.
4) `multirate`: A particioned solution of the oscillator using preCICE 2. Implementing waveform iteration ideas from[^1] in the adapter code. Additionally illustrating improvements for subcycling.
5) `precice-3`: A particioned solution of the oscillator using preCICE 3. The native API from preCICE 3 is used to apply waveform iteration.

[^1]: B. Rodenberg, **in preparation**
[^2]: V. Schüller, B. Rodenberg, B. Uekermann and H. Bungartz, A Simple Test Case for Convergence Order in Time and Energy Conservation of Black-Box Coupling Schemes, in: WCCM-APCOM2022. URL https://www.scipedia.com/public/Rodenberg_2022a
