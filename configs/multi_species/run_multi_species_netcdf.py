## ------------------------------------------------------------------------
## Hardcoded path (to be improved)
## ------------------------------------------------------------------------
import sys
sys.path.append("../../aquacosm1D_lib/")

import numpy as np
from netCDF4 import Dataset

from aquacosm1D import *
from eddy_diffusivity_models import *

## ------------------------------------------------------------------------
## --- Constants ---
## ------------------------------------------------------------------------

SECONDS_PER_DAY   = 86400        # [s day-1]
HOURS_TO_SECONDS  = 3600         # [s hour-1]

dt        = 5.0      # time step [s]
Ndays     = 100       # simulation length [days]
Npts      = 200      # number of particles (aquacosms)
Nspecies  = 300      # number of phytoplankton species

## ------------------------------------------------------------------------
## --- Derived constants ---
## ------------------------------------------------------------------------

Nloops    = int(24 * HOURS_TO_SECONDS * Ndays / dt)   # total time steps
Nstore    = int(6 * HOURS_TO_SECONDS / dt)            # store every X hours
Nscalars  = Nspecies + 1                              # species + nutrient

EXP       = "T02"

## ------------------------------------------------------------------------
## --- RNG (random generator) ---
## ------------------------------------------------------------------------

rng = np.random.default_rng(1234567)

## ------------------------------------------------------------------------
## --- Utility functions ---
## ------------------------------------------------------------------------

def uniform_samples(n, min_val, max_val):
    """Generate uniformly distributed samples."""
    return rng.uniform(min_val, max_val, n)

## ------------------------------------------------------------------------
## --- Set-up the physical environment ---
## ------------------------------------------------------------------------

wc = water_column(
    kappa = Gompertz_Eddy_Diffusion_1em1,
    max_depth = 50,                 # [m]
    short_wave_radiation = 100      # [W m-2]
)

## Transport operator (stochastic vertical motion)
Transport = set_up_transport(wc, dt, SODE='Milstein')


## ------------------------------------------------------------------------
## Initial conditions  and parameters (to be varied in sensitivity analysis)
## ------------------------------------------------------------------------
P0  = 0.2     # initial phytoplankton carbon [mg C m-3]
N0  = 1.0     # initial nutrient concentration [mmol N m-3]
p   = 1e-5    # aquacosm interaction parameter
radius = 10.0 # interaction radius [m]

## Create particles (Lagrangian aquacosms)
Particles = create_particles(Npts, Nscalars, wc)
Particles[:, 2:-1] = P0
Particles[:, -1]   = N0


## ------------------------------------------------------------------------
## --- Species traits (uniform distributions) ---
## ------------------------------------------------------------------------

Alpha            = uniform_samples(Nspecies, 0.1e-5, 2e-5)
MaxPhotoRate     = uniform_samples(Nspecies, 1.0, 8.0)
BasalMetabolism  = uniform_samples(Nspecies, 0.01, 0.3)
HalfSaturation   = uniform_samples(Nspecies, 0.1, 4.0)

## ------------------------------------------------------------------------
## --- Biological reaction model ---
## ------------------------------------------------------------------------

React = set_up_reaction(
    wc,
    dt,
    Chemostat_multi_species,
    Alpha           = Alpha,
    MaxPhotoRate    = MaxPhotoRate,
    BasalMetabolism = BasalMetabolism,
    HalfSaturation  = HalfSaturation,
)

## ------------------------------------------------------------------------
## --- Particle interaction / diffusion ---
## ------------------------------------------------------------------------

Diffuse = set_up_diffusion(
    Npts,
    Nscalars,
    radius = radius,
    p      = p,
    dt     = dt,
    wc     = wc
)

## ------------------------------------------------------------------------
## --- NetCDF file setup (incremental writing) ---
## ------------------------------------------------------------------------

basename  = f"{EXP}_{Nspecies}.p{Diffuse.p:1.0e}.nc"
ds = Dataset(basename, "w", format="NETCDF4")

## Define dimensions
ds.createDimension("time", None)    # unlimited dimension
ds.createDimension("particle", Npts)
ds.createDimension("nutrient", 1) 
ds.createDimension("species", Nspecies)

## Coordinate variables
time_var     = ds.createVariable("time", "f8", ("time",))
particle_var = ds.createVariable("particle", "i4", ("particle",))
species_var   = ds.createVariable("species", "i4", ("species",))
nutrient_var   = ds.createVariable("nutrient", "i4", ("nutrient",))

## ID variable
ID_var = ds.createVariable(
    "ID",
    "i4",
    ("particle")
)
  
## Depth variable
depth_var = ds.createVariable(
    "Depth",
    "f8",
    ("time","particle"),
)

## Main data variables
P_var = ds.createVariable(
    "Phytoplankton",
    "f8",
    ("time", "particle", "species"),
    zlib=True,                 # compression
    complevel=4
)
P_var.longname = "Phytoplankton carbon"
P_var.units = "mg C/m3"

N_var = ds.createVariable(
    "Nutrient",
    "f8",
    ("time", "particle"),
    zlib=True,                 # compression
    complevel=4
)
N_var.longname = "Nitrate concentration"
N_var.units = "mmol N/m3"


## ------------------------------------------------------------------------
## --- Metadata (self-describing output) ---
## ------------------------------------------------------------------------

ds.description = "Aquacosm multi-species Lagrangian simulation"
ds.dt          = dt
ds.Ndays       = Ndays
ds.Nspecies    = Nspecies
ds.Npts        = Npts

## Store species trait parameters
ds.createVariable("Alpha", "f8", ("species",))[:] = Alpha
ds.createVariable("MaxPhotoRate", "f8", ("species",))[:] = MaxPhotoRate
ds.createVariable("BasalMetabolism", "f8", ("species",))[:] = BasalMetabolism
ds.createVariable("HalfSaturation", "f8", ("species",))[:] = HalfSaturation

## Fill coordinate arrays
particle_var[:] = np.arange(Npts)
species_var[:]  = np.arange(Nspecies)
nutrient_var[:] = np.arange(1)

## ------------------------------------------------------------------------
## --- Run the simulation ---
## ------------------------------------------------------------------------

frame = 0                 # index along "time" dimension at storage frequency
ID_var   = Particles[:,0] # store the ID for each particle (invariant)

for loop in range(Nloops):

    ## ------------------------------------------------------------
    ## --- Store results every Nstore steps ---
    ## ------------------------------------------------------------

    if loop % Nstore == 0:

        ## Convert simulation time to days
        time_days = loop * dt / SECONDS_PER_DAY

        ## Save to NetCDF (incremental write)
        depth_var            = Particles[:,1]
        time_var[frame]      = time_days
        P_var[frame, :, :]   = Particles[:,2:-1]
        N_var[frame, :]      = Particles[:,-1]

        ## Flush to disk periodically (safety + performance balance)
        if frame % 10 == 0:
            ds.sync()

        frame += 1
        print('Days elapsed: ', time_days)

    ## ------------------------------------------------------------
    ## --- Model time step ---
    ## ------------------------------------------------------------

    ## Vertical transport
    Transport(Particles)

    ## Biological reactions
    React(Particles)

    ## Particle interaction / diffusion
    Diffuse(Particles)

## ------------------------------------------------------------------------
## --- Finalize simulation ---
## ------------------------------------------------------------------------

ds.close()

print("Simulation complete.")
print("Output saved to:", basename)

