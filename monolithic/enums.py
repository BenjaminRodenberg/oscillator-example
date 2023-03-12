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