import numpy as np
import matplotlib.pyplot as plt

from aquacosm1D_reactions import (
    set_up_reaction,
    NoReactions,
    Sverdrup,
    Sverdrup_incl_K,
    SimpleBFM,
    SimpleBFM_onlyC,
    BioShading,
    BioShading_onlyC,
    Chemostat_multi_species,
)

from aquacosm1D_diffusion import (
    set_up_diffusion,
    sort_by_depth,
)

from aquacosm1D_transport import set_up_transport

from aquacosm1D_watercolumn import (
    water_column,
    water_column_netcdf,
)

from aquacosm1D_utilities import (
    Aquacosm1D_Particles,
    create_particles,
    sort_by_index,
    sort_by_depth,
    average_between,
    gaussian_estimate_field,
    gaussian_average_between,
)

# NumPy floating point behaviour
np.seterr(divide="raise")
np.seterr(invalid="raise")
np.seterr(over="warn")
np.seterr(under="warn")

