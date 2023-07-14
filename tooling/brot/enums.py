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
    ZERO = "Zero"  # No interpolation. Just piecewise constant
    LAGRANGE = "Lagrange"
    LAGRANGE_LINEAR = "Lagrange_linear"  # Lagrange interpolation, first degree
    LAGRANGE_QUADRATIC = "Lagrange_quadratic"  # Lagrange interpolation, second degree
    HERMITE_CUBIC = "Hermite_cubic"  # Hermite interpolation, third degree
    HERMITE = "Hermite"
    BSPLINE = "B_spline"
    BSPLINE_LINEAR = "B_spline_linear"  # Linear B-splines
    BSPLINE_CUBIC = "B_spline_cubic"  # Cubic B-splines
    BSPLINE_TEN = "B_spline_ten"  # tenth degree B-splines


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

class AccelerationSchemes(Enum):
    NONE = "None"
    CONSTANT = "Constant"
    REDUCED_QUASI_NEWTON = "rQN"
    FULL_QUASI_NEWTON = "QN"