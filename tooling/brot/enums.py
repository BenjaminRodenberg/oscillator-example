from enum import Enum

class Cases(Enum):
    MONOLITHIC = "monolithic"
    PRECICE2 = "precice2"
    WAVEFORM = "waveform"
    MULTIRATE = "multirate"
    PRECICE3 = "precice3"


class TimeSteppingSchemes(Enum):
    NEWMARK_BETA = "Newmark_beta"
    GENERALIZED_ALPHA = "generalized_alpha"
    RUNGE_KUTTA_4 = "runge_kutta_4"


class ReadWaveformSchemes(Enum):
    ZERO = "Zero"  # No interpolation. Just piecewise constant
    LAGRANGE = "Lagrange"
    LAGRANGE_LINEAR = "Lagrange_linear"  # Lagrange interpolation, first order
    LAGRANGE_QUADRATIC = "Lagrange_quadratic"  # Lagrange interpolation, second order
    HERMITE_CUBIC = "Hermite_cubic"  # Hermite interpolation, third order
    HERMITE = "Hermite"
    BSPLINE = "BSpline"
    BSPLINE_LINEAR = "BSpline_linear"  # Linear BSplines
    BSPLINE_CUBIC = "BSpline_cubic"  # Cubic BSplines
    BSPLINE_TEN = "BSpline_ten"  # tenth order BSplines

class MultirateMode(Enum):
    NONE = "None"  # Don't use any multirate
    SUBCYCLING = "Subcycling"  # Use fixed window size and apply subcycling on both sides
    SUBSTEPS = "substeps"  # Perform substeps
    TWO_SUBSTEPS = "Two_substeps"  # Reduce window size and apply subcycling to always use two substeps
    FOUR_SUBSTEPS = "Four_substeps"  # Reduce window size and apply subcycling to always use four substeps
    EIGHT_SUBSTEPS = "Eight_substeps"  # Reduce window size and apply subcycling to always use eight substeps
    MULTIRATE = "Multirate"  # Use fixed window size, perform 4 substeps on the left side and apply subcycling on the right


class ParticipantNames(Enum):
    MASS_LEFT = "Mass-Left"
    MASS_RIGHT = "Mass-Right"


class DataNames(Enum):
    FORCE_LEFT = "Force-Left"
    FORCE_LEFT_1 = "Force-Left-1"
    FORCE_LEFT_2 = "Force-Left-2"
    FORCE_RIGHT = "Force-Right"
    FORCE_RIGHT_1 = "Force-Right-1"
    FORCE_RIGHT_2 = "Force-Right-2"


class MeshNames(Enum):
    MASS_LEFT_MESH = "Mass-Left-Mesh"
    MASS_RIGHT_MESH = "Mass-Right-Mesh"