import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D_reactions   import set_up_reaction, NoReactions, BioShading_onlyC, Sverdrup, Sverdrup_incl_K
from aquacosm1D_watercolumn import water_column, water_column_netcdf
from eulerian1D_diffusion import set_up_eulerian_diffusion
from aquacosm1D_utilities import Aquacosm1D_Particles
from pylab import *
from plot_eul_aqc_lib import *
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy.interpolate import splev
import xarray as xr
import params

ion()
react_params = params.reactions()
# seed(1234567) moving below to reset seed for each run

#------------------------------------------------------------
dt        = 5. # time step in seconds
Ndays     = 22 #length of the simulation
Nloops    = int(24*3600  *  Ndays  / dt)
Nstore    = int(0.5*3600 / dt) #store the particles every Nshow time steps
Nconsole  = int(6*3600 / dt) # frequency of writing to the console
Nscalars  = 1    #number of scalars carried by each particle

# physical inputs to loop through for sensitivity tests
# (corresponding to CROCO 1D runs)
amplitudes = [0.03] #[0, 0.01, 0.02, 0.03, 0.04]
mean_taus = [0] #[0, 0.05]
mlds = [10] #[10, 25]
Qswmaxs = [250] #250 isthe special case of constant heat flux, otherwise it is diurnal peaks at Qswmax   

for amplitude in amplitudes:
    for mean_tau in mean_taus:
        for mld in mlds:
            for Qswmax in Qswmaxs:  
        
                seed(1234567)
                 
                crocodir='../physics/'
                crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
                crocofile=crocodir+crocofilename                
                                
                wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
                
                Diffuse = set_up_eulerian_diffusion(wc,dt,Nscalars)
                                
                # Reactions
                # BioShading_onlyC
                React = set_up_reaction(wc, dt, BioShading_onlyC,
                                    LightDecay=react_params.LightDecay,
                                    MaxPhotoRate = react_params.MaxPhotoRate, 
                                    BasalMetabolism = react_params.BasalMetabolism,
                                    Chl_C = react_params.Chl_C,
                                    CrowdingMortality = react_params.CrowdingMortality,
                                    CrowdingHalfSaturation = react_params.CrowdingHalfSaturation,
                                    Chl_light_abs =react_params.Chl_light_abs)
                #
                
                # Here's where we initialise the chlorophyll
                data_croco=Dataset(crocofile)
                temp_croco=data_croco.variables['temp'][:,:,0,0]
                zt=data_croco.variables['deptht'][:]
                time_croco = data_croco.variables['time_counter'][:]/86400
                data_croco.close()
                # create a constant chlorophyll ini over surface layer
                z_therm_croco=get_z_therm_croco(time_croco,zt,temp_croco)
                chl_ini=np.zeros(np.shape(zt))+1e-20 # mg/m3
                chl_ini[zt<z_therm_croco[0]]=1
                # add end points - the end points are not used and will be fixed by b.c.
                chl_ini = np.concatenate(([chl_ini[0]],chl_ini,[chl_ini[-1]])) 
                # convert Chl to C using fixed ratio
                c_ini=chl_ini/React.Chl_C
                
                # create the eulerian tracer array in the same format as the aquacosms so we can 
                # use the same reactions library as used in the aquacosms
                Npts=len(Diffuse.zc) 
                Tracers      = Aquacosm1D_Particles(
                        zeros((Npts, Nscalars+2), dtype='float64')
                        )
                Tracers[:,0] = np.arange(Npts)
                Tracers[:,1] = Diffuse.zc
                Tracers[:,2] = c_ini
                
                fname_out='eulerian_r'+str(react_params.MaxPhotoRate)+'_b'+str(react_params.BasalMetabolism)+'_c'+str(react_params.Chl_light_abs)+'_a'+str(react_params.CrowdingMortality)+'_l'+str(react_params.LightDecay)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
                    
                print('working on: ', fname_out)
                
                Cstore = []
                Tstore    = []   #output time stamps
                starttime = wc.times[0]
                wc.set_current_time(starttime)
                for loop in range(Nloops):
                    #-------------------------
                    
                    if loop%Nstore == 0:
                        Cstore.append(Tracers.copy())
                        Tstore.append(wc.get_current_time())
                        
                    if loop % Nconsole == 0:        
                        print("Days elapsed:", loop*dt/(24*3600.),
                              'Mean C:', mean(Tracers[:,2]))
                        
                    Diffuse(Tracers)
                    React(Tracers)
                    
                    wc.increment_current_time(dt)
                
                Cstore = array(Cstore)
                Tstore = array(Tstore)
                
                # write a netcdf output file
                ds = xr.Dataset(data_vars=dict(chl=(["time","depth"],Cstore[:,:,2]*React.Chl_C)),
                                                       coords=dict(depth=(["depth"],Diffuse.zc),time=(["time"],Tstore)))
                ds.time.attrs['units'] = 'seconds since start of simulation' # not a real time period as an analytical experiment
                ds.chl.attrs['long_name'] = 'chlorophyll concentration'
                ds.chl.attrs['units'] = 'mg m^-3'
                ds.depth.attrs['long_name'] = 'depth positive downward'
                ds.depth.attrs['units'] = 'm'
                # add global attributes telling us what the input parameters were
                ds.attrs["Chl_C"] = str(React.Chl_C)
                # write the output
                ds.to_netcdf(fname_out)
                
