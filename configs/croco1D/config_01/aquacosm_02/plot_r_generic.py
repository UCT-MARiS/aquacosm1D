import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
import xarray as xr
from pathlib import Path
from scipy.interpolate import interp1d
import params
from plot_eul_aqc_lib import *
ion()
react_params = params.reactions()
def get_r(React):
    
    # set up aquacosm array in correct format for calling reactions library
    depth=int(wc.max_depth)
    Npts=depth+1 # 1m increments
    Nscalars=1
    Particles = create_particles(Npts, Nscalars, React.wc)
    Particles      = Aquacosm1D_Particles(
        zeros((Npts, Nscalars+2), dtype='float64')
        )
    Particles[:,1] = np.linspace(0, depth,num=Npts)
    sort_by_depth(Particles)
    Particles[:,0] = arange(Npts)
    
    # define Chlorophyll values to compute reaction rates
    # Chls = np.append([0.1],np.linspace(1, 20,num=20))
    Chls = np.array([0.001,0.01,0.1,1,10,100])
    
    RRates_normalised = np.zeros((Npts,len(Chls)))
    for ii,Chl in enumerate(Chls):
        # update particles array with Chl
        Chl_array=np.zeros((Npts,))+Chl
        C=Chl_array/React.Chl_C 
        Particles[:,2]=C
        RRates=React.current_model(Particles, React.wc, 0)
        RRates_normalised[:,ii]=RRates[:,0]/Particles[:,2]
    
    RRates_normalised = RRates_normalised*3600*24 # now days^-1
    
    depths=Particles[:,1]
    # Chls,depths=np.meshgrid(Chls,depths)
    
    return Chls,depths,RRates_normalised 
    
def do_the_plot(Chl, depth, r):
    
    #    
    depth=np.append(depth-0.5,depth[-1]+0.5)
    # Chl=np.append(Chl-0.5,Chl[-1]+0.5)
    Chl=np.append(Chl/3,Chl[-1]*3)
    
    figure(figsize=(5,5))
    ax = subplot(1,1,1) 
    
    ax.set_position(  (0., 0., 1., 0.5))
    # img = ax.scatter(Chl, depth, c=r,
    #                     cmap=cm.viridis, linewidth=0, s=40)
    img = ax.pcolormesh(Chl, depth, r, cmap=cm.RdBu_r)
    img.set_clim(-1, 1)
    ax.set_ylim(20.5, -0.5)
    ax.set_yticks(np.linspace(0,20,num=11))

    # ax.set_xlim(0.003,300)
    ax.set_xscale('log')
    #Let's get rid of the scientific notation on the y-axis
    fmt = mpl.ticker.StrMethodFormatter("{x:g}")
    ax.xaxis.set_major_formatter(fmt)
    ax.minorticks_off()
    
    for C in Chl:
        ax.plot([C,C],[depth[0],depth[-1]],color='k',linewidth=1)
    
    ax.set_ylabel('Depth (m)', fontsize=12,)
    ax.set_xlabel('Chl mg m$^{-3}$)', fontsize=12,)
    
    cbarax = gcf().add_axes((1.02, 0., 0.015, 0.5)) # [left, bottom, width, height] 
    cbar = colorbar(img, cbarax)
    cbar.set_label('R (days$^{-1}$)', fontsize=12)
    
    fname_out='plot_r_generic_'+type(React.current_model).__name__+'_'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'.jpg'
    plt.savefig(fname_out,dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":
    
    # dummy values needed here to create a water column object - not actually used
    dt = 5             
    kappa = 0.001
    max_depth=20
    wc = water_column(kappa, max_depth, short_wave_radiation=250)
    
    React = set_up_reaction(wc, dt, BioShading_onlyC,
                            LightDecay=react_params.LightDecay,
                            MaxPhotoRate = react_params.MaxPhotoRate, 
                            BasalMetabolism = react_params.BasalMetabolism,
                            Chl_C = react_params.Chl_C,
                            CrowdingMortality = react_params.CrowdingMortality,
                            Chl_light_abs = react_params.Chl_light_abs)
    
    Chl, depth, r=get_r(React)
    
    do_the_plot(Chl, depth, r)

    
    
