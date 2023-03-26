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
    LAGRANGE_LINEAR = "Lagrange interpolation, first order"
    LAGRANGE_QUADRATIC = "Lagrange interpolation, second order"
    HERMITE_CUBIC = "Hermite interpolation, third order"


class ParticipantNames(Enum):
    MASS_LEFT = "Mass-Left"
    MASS_RIGHT = "Mass-Right"


class DataNames(Enum):
    FORCE_LEFT = "Force-Left"
    FORCE_RIGHT = "Force-Right"
    FORCE_RIGHT_1 = "Force-Right-1"
    FORCE_RIGHT_2 = "Force-Right-2"


class MeshNames(Enum):
    MASS_LEFT_MESH = "Mass-Left-Mesh"
    MASS_RIGHT_MESH = "Mass-Right-Mesh"