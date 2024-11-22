from enum import Enum


class ReadWaveformSchemes(Enum):
    LAGRANGE = "Lagrange"  # Lagrange interpolation, degree equals number of points
    HERMITE = "Hermite"  # Hermite interpolation, third degree
    BSPLINE = "BSpline"


class AccelerationSchemes(Enum):
    NONE = "None"
    CONSTANT = "Constant"
    REDUCED_QUASI_NEWTON = "rQN"
    FULL_QUASI_NEWTON = "QN"