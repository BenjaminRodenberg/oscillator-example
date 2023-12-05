from enum import Enum

class Cases(Enum):
    MONOLITHIC = "monolithic"
    PRECICE2 = "precice2"
    WAVEFORM = "waveform"
    MULTIRATE = "multirate"
    ACCELERATION = "acceleration"
    PRECICE3 = "precice3"


class TimeSteppingSchemes(Enum):
    NEWMARK_BETA = "Newmark_beta"
    GENERALIZED_ALPHA = "generalized_alpha"
    RUNGE_KUTTA_4 = "runge_kutta_4"
    Radau_IIA = "radauIIA"  # 5th order implicit


class ReadWaveformSchemes(Enum):
    LAGRANGE = "Lagrange"  # Lagrange interpolation, degree equals number of points
    HERMITE = "Hermite"  # Hermite interpolation, third degree
    BSPLINE = "BSpline"


class ParticipantNames(Enum):
    MASS_LEFT = "Mass-Left"
    MASS_RIGHT = "Mass-Right"


class DataNames(Enum):
    FORCE_LEFT = "Force-Left"
    FORCE_RIGHT = "Force-Right"


class MeshNames(Enum):
    MASS_LEFT_MESH = "Mass-Left-Mesh"
    MASS_RIGHT_MESH = "Mass-Right-Mesh"

class AccelerationSchemes(Enum):
    NONE = "None"
    CONSTANT = "Constant"
    REDUCED_QUASI_NEWTON = "rQN"
    FULL_QUASI_NEWTON = "QN"