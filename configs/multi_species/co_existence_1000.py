import sys
sys.path.append("../../aquacosm1D_lib/")
from aquacosm1D.aquacosm1D_lib.aquacosm1D import *
from eddy_diffusivity_models import *
from numpy.random import random, seed
ion()
#-----------------------------------------------------------------------
seed(1234567)

dt        = 5. # time step in seconds
Ndays     = 100 #4*365 #length of the simulation
Nloops    = int(24*3600  *  Ndays  / dt)
Nstore    = int(1*3600 / dt) #store the particles every Nshow time steps
Npts      = 200  #number of particles
Nspecies  = 1000 #number of distinct plankton species
Nscalars  = Nspecies+1   #number of scalars carried by each particle
Pstore    = []   #list that stores the particles every Nstore time steps
P25store  = []   #list that stores the particles every Nstore time steps

days = arange(0, Nloops, Nstore)*dt/86400 #time index of Pstore in days

#------------------------------------------------------------------------
def daily_light_cycle(t):
    s = 2*pi*t/86400
    A = 100 # W/m^2 of equivalent constant light source
    return pi*A*sin(s)*(sin(s)>0.)

#------------------------------------------------------------------------
#-- Set-up the simulation
#--
wc = water_column(
   #kappa                = Time_Varying_Gompertz_Eddy_Diffusion,
    kappa                = Gompertz_Eddy_Diffusion_1em1,
   #kappa                = Gompertz_Eddy_Diffusion_1em2,
   #kappa                = 0.1,
    max_depth            = 50,
    short_wave_radiation = 100, # W/m^2
   #short_wave_radiation = daily_light_cycle, 
)
#--
Particles = create_particles(Npts, Nscalars, wc)
Particles[:,2:-1] = 0.2/Nspecies #initial C content of each particle
Particles[:,  -1] = 1.  #initial R content of each particle

Transport = set_up_transport(wc, dt, SODE='Milstein')

AlphaEpsilon    = (0.1e-5+1.4e-5*random(Nspecies))*(0.4/0.217)
MaxPhotoRate    = 1.+random(Nspecies)*9
BasalMetabolism = 0.01+random(Nspecies)*0.2
HalfSaturation  = 0.01+random(Nspecies)*2
React = set_up_reaction(wc, dt, Chemostat_multi_species,
    AlphaEpsilon    = AlphaEpsilon,
    MaxPhotoRate    = MaxPhotoRate,
    BasalMetabolism = BasalMetabolism, 
    HalfSaturation  = HalfSaturation,
)

Diffuse = set_up_diffusion(Npts, Nscalars,
                           radius=10.,
                           p=1.e-5, ###################
                           dt=dt,
                           wc=wc)

#------------------------------------------------------------------------
#-- Run the simulation
#--
for loop in range(Nloops):
    #-------------------------
    if loop % Nstore == 0:
        Ptmp = Particles.copy()
        sort_by_index(Ptmp)
        Pstore.append(Ptmp)
        print('Days elapsed: ', loop*dt/(24*3600.))
    Diffuse(Particles)
    React(Particles)
    Transport(Particles)
    wc.increment_current_time(dt)


#------------------------------------------------------------------------
#-- Store results
#--
basename = f"Species_{Nspecies}.p{Diffuse.p:1.0e}.Gompertz_Ktop_1.e-1.Kbt_1.e-4"
print("Saving")
Pstore = array(Pstore)
Pstore.tofile(basename+".Pstore.pydat")
AlphaEpsilon.tofile(basename+".AlphaEpsilon.pydat")
MaxPhotoRate.tofile(basename+".MaxPhotoRate.pydat")
BasalMetabolism.tofile(basename+".BasalMetabolism.pydat")
HalfSaturation.tofile(basename+".HalfSaturation.pydat")
print("All done!")

