### ------------------------------------------------------------------------
### Hardcoded path (to be improved)
### ------------------------------------------------------------------------

import sys
sys.path.append("../../aquacosm1D_lib/")
import numpy as np
from netCDF4 import Dataset
from aquacosm1D import *
from eddy_diffusivity_models import *
from multiprocessing import Pool
from tqdm import tqdm

### ------------------------------------------------------------------------
### --- Constants ---
### ------------------------------------------------------------------------

SECONDS_PER_DAY   = 86400        # [s day-1]
HOURS_TO_SECONDS  = 3600         # [s hour-1]
dt        = 5.0      # time step [s]
Ndays     = 40       # simulation length [days]
Npts      = 200      # number of particles (aquacosms)
Nspecies  = 100      # number of phytoplankton species
EXP       = "SENS02"

### ------------------------------------------------------------------------
### --- Derived constants ---
### ------------------------------------------------------------------------

Nloops    = int(24 * HOURS_TO_SECONDS * Ndays / dt)   # total time steps
Nstore    = int(6 * HOURS_TO_SECONDS / dt)            # store every X hours
Nscalars  = Nspecies + 1                              # species + nutrient

### ------------------------------------------------------------------------
### --- RNG (random generator) ---
### ------------------------------------------------------------------------

rng = np.random.default_rng(441998)

### ------------------------------------------------------------------------
### --- Utility functions ---
### ------------------------------------------------------------------------

def uniform_samples(n, min_val, max_val):
    "Generate uniformly distributed samples."
    return rng.uniform(min_val, max_val, n)

### ------------------------------------------------------------------------
### --- Set-up the physical environment ---
### ------------------------------------------------------------------------

wc = water_column(
    kappa = Gompertz_Eddy_Diffusion_1em1,
    max_depth = 50,                 # [m]
    short_wave_radiation = 100      # [W m-2]
)

### Transport operator (stochastic vertical motion)

Transport = set_up_transport(wc, dt, SODE='Milstein')

### ------------------------------------------------------------------------
### --- Species traits (uniform distributions) ---
### ------------------------------------------------------------------------

Alpha            = uniform_samples(Nspecies, 0.1e-5, 2e-5)
MaxPhotoRate     = uniform_samples(Nspecies, 1.0, 8.0)
BasalMetabolism  = uniform_samples(Nspecies, 0.01, 0.3)
HalfSaturation   = uniform_samples(Nspecies, 0.1, 4.0)

### ------------------------------------------------------------------------
### --- Simulation function (single experiment) ---
### ------------------------------------------------------------------------

def run_simulation(params):

    ### ------------------------------------------------------------------------
    ### Initial conditions  and parameters 
    ### ------------------------------------------------------------------------

    P0, N0, p, idx = params

    ### Create particles (Lagrangian aquacosms)

    Particles = create_particles(Npts, Nscalars, wc)
    Particles[:, 2:-1] = P0
    Particles[:, -1]   = N0

    ### ------------------------------------------------------------------------
    ### --- Biological reaction model ---
    ### ------------------------------------------------------------------------

    React = set_up_reaction(
        wc,
        dt,
        Chemostat_multi_species,
        Alpha           = Alpha,
        MaxPhotoRate    = MaxPhotoRate,
        BasalMetabolism = BasalMetabolism,
        HalfSaturation  = HalfSaturation,
    )

    ### ------------------------------------------------------------------------
    ### --- Particle interaction / diffusion ---
    ### ------------------------------------------------------------------------

    Diffuse = set_up_diffusion(
        Npts,
        Nscalars,
        radius = 10.,
        p      = p,
        dt     = dt,
        wc     = wc
    )

    ### ------------------------------------------------------------------------
    ### --- NetCDF file setup (incremental writing) ---
    ### ------------------------------------------------------------------------

    basename = f"{EXP}_{idx:03d}_P0{P0:.3f}_N0{N0:.3f}_p{p:.2e}.nc"
    ds = Dataset(basename, "w", format="NETCDF4")

    ### Define dimensions

    ds.createDimension("time", None)
    ds.createDimension("particle", Npts)
    ds.createDimension("nutrient", 1)
    ds.createDimension("species", Nspecies)

    ### Coordinate variables

    time_var     = ds.createVariable("time", "f8", ("time",))
    particle_var = ds.createVariable("particle", "i4", ("particle",))
    species_var  = ds.createVariable("species", "i4", ("species",))
    nutrient_var = ds.createVariable("nutrient", "i4", ("nutrient",))

    ID_var = ds.createVariable("ID", "i4", ("particle"))

    depth_var = ds.createVariable("Depth", "f8", ("time","particle"))

    P_var = ds.createVariable(
        "Phytoplankton",
        "f8",
        ("time", "particle", "species"),
        zlib=True,
        complevel=4,
        chunksizes=(1, Npts, Nspecies)
    )

    N_var = ds.createVariable(
        "Nutrient",
        "f8",
        ("time", "particle"),
        zlib=True,
        complevel=4,
        chunksizes=(1, Npts)
    )

    ### ------------------------------------------------------------------------
    ### --- Metadata ---
    ### ------------------------------------------------------------------------

    ds.description = "Aquacosm multi-species Lagrangian simulation"
    ds.dt          = dt
    ds.Ndays       = Ndays
    ds.Nspecies    = Nspecies
    ds.Npts        = Npts

    ds.createVariable("Alpha", "f8", ("species",))[:] = Alpha
    ds.createVariable("MaxPhotoRate", "f8", ("species",))[:] = MaxPhotoRate
    ds.createVariable("BasalMetabolism", "f8", ("species",))[:] = BasalMetabolism
    ds.createVariable("HalfSaturation", "f8", ("species",))[:] = HalfSaturation

    particle_var[:] = np.arange(Npts)
    species_var[:]  = np.arange(Nspecies)
    nutrient_var[:] = np.arange(1)

    P_var.longname = "Phytoplankton carbon"
    P_var.units = "mg C/m3"
    N_var.longname = "Nitrate concentration"
    N_var.units = "mmol N/m3"
    time_var.units = "days"

    ### ------------------------------------------------------------------------
    ### --- Run the simulation ---
    ### ------------------------------------------------------------------------

    frame = 0
    ID_var[:] = Particles[:,0]

    for loop in range(Nloops):

        # transport step
        Transport(Particles)

        # biological reaction
        React(Particles)

        # interaction / diffusion
        Diffuse(Particles)

        # store output
        if loop % Nstore == 0:
            time_var[frame] = loop * dt / SECONDS_PER_DAY
            depth_var[frame, :] = Particles[:,1]
            P_var[frame, :, :]  = Particles[:,2:-1]
            N_var[frame, :]     = Particles[:,-1]
            frame += 1
            if frame % 10 == 0:
                ds.sync() # flush to disk every 10 output ops


    ### ------------------------------------------------------------------------
    ### --- Finalize simulation ---
    ### ------------------------------------------------------------------------

    ds.close()
    #print(f"Completed: {basename}")

### ------------------------------------------------------------------------
### --- Sensitivity analysis setup ---
### ------------------------------------------------------------------------

Nexp = 256   # number of experiments

P0_vals     = uniform_samples(Nexp, 0.01, 1.0)
N0_vals     = uniform_samples(Nexp, 0.1, 4.0)
exponent    = np.floor(uniform_samples(Nexp, 3, 8))
p_vals      = np.power(10.,-exponent)
#radius_vals = uniform_samples(Nexp, 5.0, 20.0)

param_list = [
    (P0_vals[i], N0_vals[i], p_vals[i], i)
    for i in range(Nexp)
]

### ------------------------------------------------------------------------
### --- Parallel execution (X CPUs) ---
### ------------------------------------------------------------------------
CPU = 12
if __name__ == "__main__":

    with Pool(processes=CPU) as pool:
        results = []
        for _ in tqdm(pool.imap(run_simulation, param_list),
                      total=len(param_list),
                      desc="Running simulations"):
            results.append(_)

    print("All simulations complete.")
