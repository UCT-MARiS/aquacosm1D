from pylab import *

#------------------------------------------------------------------------

def Gompertz_Eddy_Diffusion_1em1(z, t):
    Kmax    = 1.e-1
    Kmin    = 1.e-4
    steep   = 1/4
    MLdepth = 25
    return Kmin + (Kmax-Kmin)*( 1 - exp(-exp((MLdepth-z)*steep)) )

#------------------------------------------------------------------------

def Gompertz_Eddy_Diffusion_1em2(z, t):
    Kmax    = 1.e-2
    Kmin    = 1.e-4
    steep   = 1/4
    MLdepth = 25
    return Kmin + (Kmax-Kmin)*( 1 - exp(-exp((MLdepth-z)*steep)) )

#------------------------------------------------------------------------
def smooth_sawtooth(t):
    """This is periodic of period 2*pi, ranging between 0 and 1.
    """
    return (1+(sin(t)+0.3*sin(2*t)+0.1*sin(3*t))/1.1283820906416413)/2

def Time_Varying_Gompertz_Eddy_Diffusion(z, t):
    Period = 7 #days
    Ksurfmax   = 1.e-1
    Ksurfmin   = 1.e-2
    MLdepthmax = 40
    MLdepthmin = 15
    Kdeep      = 1.e-4
    steep      = 1/4
    Ksurf   = (Ksurfmin +
               smooth_sawtooth(2*pi*t/(Period*86400)) *
               (Ksurfmax-Ksurfmin)
               )
    MLdepth = (MLdepthmin +
               smooth_sawtooth(2*pi*t/(Period*86400)) *
               (MLdepthmax-MLdepthmin)
               )
    return Kdeep + (Ksurf-Kdeep)*( -expm1(-exp((MLdepth-z)*steep)) )

