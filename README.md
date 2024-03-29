# Aquacosm 1D 
(C) 2019-2021 Francesco Paparella (NYU-AD).
Further development by Giles Fearon and Marcello Vichi (UCT)

This library implements the aquacosm approch described in

_Paparella, F., Vichi, M., 2020. Stirring, Mixing, Growing: 
Microscale Processes Change Larger Scale Phytoplankton Dynamics. 
Front. Mar. Sci. 7. [https://doi.org/10.3389/fmars.2020.00654]_

It is still a work in progress and not all components are functional.
The actual code is found in `aquacosm1D_lib`. It is mostly a pure
python code, but it links to some C code that implements diffusion.
Because of this, it can't be used out of the box: one has to compile
the C code into a shared library, first. On a linux box, with the gcc
compiler, this is as easy as issuing the following command in the
`aquacosm1D_lib` folder:

`gcc -std=c99 -pedantic -Wall -O3 -fPIC aqc1D_diffusion.c -shared -lm -o aqc1D_diffusion.so`

The C code strictly adheres to the c99 standard, and it is linked to
python via the `ctypes` interface, which is supported on all platforms.

The `/configs/` directory contains various configurations for running the example 
published in Paparella and Vichi (2020).
Once the compilation is done, go to the `/configs/PAPA_example` folder and run
the script `PAPA_simulation_p1e-7.py`. That performs a 4-years
simulation using the eddy diffusivity data in `PAPA_6h_1Y_L75_DN.nc`, generated with the 
NEMO ocean model (it repeats 4 times the same year and will take a few hours on a
recent computer).

Most of the action occurs in the "wc" object, which is an instance of
the "water_column_netcdf" class, defined in
`aquacosm1D_lib/aquacosm1D_watercolumn.py`

The biology is found in
`aquacosm1D_lib/aquacosm1D_reactions.py`. Several examples are
present. It should be straightforward to construct others by analogy.

The file `aquacosm1D_lib/aquacosm1D_transport.py` implements the
aquacosms' Brownian motion by means of a stochastic integrator. There
are several integrators available. Use Milstein's. Upon testing, there
are no appreciable advantages with the stiff version by Yin and Gan.
The other integrators are provably inferior.

IMPORTANT NOTE: the stochastic integrators have very strong time step
constraints (see e.g.  [https://doi.org/10.4319/lom.2004.2.289] or
Visser (1997) Marine Ecology Progress Series 158, 275-281). A time
step of 5 seconds should be adequate for open ocean conditions. Too
large a time step will _not_ crash the code, but it will introduce
subtle inhomogeneities in the distribution of the particles, thus
biasing the results. If you're careless, YOU'LL SHOOT YOURSELF IN THE FOOT.

The `/comms/` dir contains various communications (presentations) for 
facilitating collaboration 

## Acknowledgements
This work was partly funded by an anonymous donor as part of the 
Whales and Climate Research Program [https://www.whalesandclimate.org/], 
as well as from the European Union’s Horizon 2020 research and innovation 
programme under the project AtlantECO (grant agreement No 862923).
