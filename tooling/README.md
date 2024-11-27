# Tooling for oscillator example from Benjamin Rodenberg's PhD thesis

The packaage `brot` (short for Benjamin Rodenberg Oscillator Tooling) provides tooling for the oscillator example presented in Benjamin Rodenberg's PhD thesis this includes:

* Different interpolation schemes in `brot.interpolation`
* Some shared parameters (mass, spring stiffness) and functions (analytical solution) for the oscillator example are provided in `brot.oscillator`
* Time stepping schemes for the partitioned case are provided in `brot.timeSteppers`
* Time stepping schemes for the monolithic case are provided in `brot.timeSteppersMonolithic`

## Installation

Please install this package by running `pip3 install --user .` from this folder.