import numpy as np

# --- Constants ---
KMIN = 1e-4
STEEP = 1 / 4
MLDEPTH = 25

# ------------------------------------------------------------------------
def Gompertz_Eddy_Diffusion_1em1(z, t):
    KMAX = 1e-1
    return KMIN + (KMAX - KMIN) * (1 - np.exp(-np.exp((MLDEPTH - z) * STEEP)))


# ------------------------------------------------------------------------
def Gompertz_Eddy_Diffusion_1em2(z, t):
    KMAX = 1e-2
    return KMIN + (KMAX - KMIN) * (1 - np.exp(-np.exp((MLDEPTH - z) * STEEP)))


# ------------------------------------------------------------------------
def smooth_sawtooth(t):
    """Periodic function with period 2*pi, ranging between 0 and 1."""
    return (
        1
        + (np.sin(t) + 0.3 * np.sin(2 * t) + 0.1 * np.sin(3 * t))
        / 1.1283820906416413
    ) / 2


# ------------------------------------------------------------------------
def Time_Varying_Gompertz_Eddy_Diffusion(z, t):
    PERIOD_DAYS = 7

    KSURF_MAX = 1e-1
    KSURF_MIN = 1e-2
    MLDEPTH_MAX = 40
    MLDEPTH_MIN = 15
    KDEEP = 1e-4

    Ksurf = KSURF_MIN + smooth_sawtooth(2 * np.pi * t / (PERIOD_DAYS * 86400)) * (
        KSURF_MAX - KSURF_MIN
    )

    MLdepth = MLDEPTH_MIN + smooth_sawtooth(2 * np.pi * t / (PERIOD_DAYS * 86400)) * (
        MLDEPTH_MAX - MLDEPTH_MIN
    )

    return KDEEP + (Ksurf - KDEEP) * (
        -np.expm1(-np.exp((MLdepth - z) * STEEP))
    )
