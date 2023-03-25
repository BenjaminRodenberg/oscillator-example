# Tooling for oscillator example from Benjamin Rodenberg's PhD thesis

The packaage `brot` (short for Benjamin Rodenberg Oscillator Tooling) provides tooling for the oscillator example presented in Benjamin Rodenberg's PhD thesis this includes:

* Enums that are uses consistently across the example for consistent naming in `brot.enums`
* Different interpolation schemes in `brot.interpolation`
* Some shared parameters (mass, spring stiffness) and functions (analytical solution) for the oscillator example are provided in `brot.oscillator`
* Tooling for creation of output files with metadata in `brot.output`
* Time stepping schemes are provided in `brot.timesteppers`

## Installation

Please install this package by running `pip3 install --user .` from this folder.