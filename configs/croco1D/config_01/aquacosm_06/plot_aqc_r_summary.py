import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
import xarray as xr
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
import params # params.py must be in this directory
ion()
        
def do_the_plot(react_params,physics_params,p):
    
    crocodir='../physics/'
    crocofilename="mean"+str(physics_params.mean_tau)+"_mld"+str(physics_params.mld)+"_amp"+str(physics_params.amplitude)+"_flx"+str(physics_params.Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename
    
    dt = 5             
    wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
    React = set_up_reaction(wc, dt, BioShading_onlyC,
                        LightDecay=react_params.LightDecay,
                        AlphaEpsilon=react_params.AlphaEpsilon,
                        MaxPhotoRate = react_params.MaxPhotoRate, 
                        BasalMetabolism = react_params.BasalMetabolism,
                        Chl_C = react_params.Chl_C,
                        CrowdingMortality = react_params.CrowdingMortality,
                        CrowdingHalfSaturation = react_params.CrowdingHalfSaturation,
                        Chl_light_abs = react_params.Chl_light_abs)
       
    fname='aquacosm_p'+"{0:1.0e}".format(p)+'_'+react_params.Name+'_'+physics_params.Name
    diagfile=fname+'_diags.nc'
    print('\n working on ' + diagfile +'\n')
        
    # get the aquacosm data
    aqcfile=fname+'.nc'
    time_aqc,z_aqc,rank_aqc,chl_aqc,_ = get_aqc_output(aqcfile)
    Nt_aqc,Nz_aqc=shape(chl_aqc)
    # repeat time along the z dimension to make arrays all the same size
    time_aqc = np.repeat(time_aqc[:, np.newaxis], Nz_aqc, axis=1)
        
    # get the aquacosm reactions (rhs of reactions equation)
    _,_,_,_,r_aqc,_,_,_,_=get_aqc_diags(diagfile)
    
    # flatten the arrays to make subsetting and computation of summary statistics easier
    time_aqc=time_aqc.flatten()
    z_aqc=z_aqc.flatten()
    chl_aqc=chl_aqc.flatten()
    r_aqc=r_aqc.flatten()
    
    # subset time to only consider data once we've had some mixing and growth away from the initial condition
    # z_aqc=z_aqc[time_aqc>7]
    # chl_aqc=chl_aqc[time_aqc>7]
    # r_aqc=r_aqc[time_aqc>7]
    
    # define the edges for binning the data
    # z
    z_bin_centres = np.linspace(0, 50, num=51, endpoint=True,)
    z_bin_edges = np.append(z_bin_centres-0.5,z_bin_centres[-1]+0.5)
    Nz = len(z_bin_centres)
    # chl
    Chl_bin_centres = np.array([0.001,0.01,0.1,1,10,100])
    Chl_bin_edges = np.append(Chl_bin_centres/3,Chl_bin_centres[-1]*3) # roughly in the middle of the ticks when plotted on log scale
    Nchl = len(Chl_bin_centres)
    
    # compute r_summary over 
    r_summary = np.zeros((Nz,Nchl))
    for z_indx in range(Nz):
        # subset on depth
        r_aqc_sub = r_aqc[(z_aqc>z_bin_edges[z_indx]) & (z_aqc<=z_bin_edges[z_indx+1])] 
        chl_aqc_sub = chl_aqc[(z_aqc>z_bin_edges[z_indx]) & (z_aqc<=z_bin_edges[z_indx+1])]
        for chl_indx in range(Nchl):    
            # subset on Chl
            r_aqc_subsub = r_aqc_sub[(chl_aqc_sub>Chl_bin_edges[chl_indx]) & (chl_aqc_sub<=Chl_bin_edges[chl_indx+1])]
            #
            # so now we have all the r data for this bin over time
            # check if we have enough data to have something meaningful
            if len(r_aqc_subsub)>10:
                r_summary[z_indx,chl_indx]=np.median(r_aqc_subsub)
            else:
                r_summary[z_indx,chl_indx]=np.nan
     
    # do the figure
    figure(figsize=(5,5))
    ax = subplot(1,1,1) 
    
    ax.set_position(  (0., 0., 1., 0.5))
    img = ax.pcolormesh(Chl_bin_edges, z_bin_edges, r_summary, cmap=cm.RdBu_r)
    img.set_clim(-1, 1)
    ax.set_ylim(20.5, -0.5)
    ax.set_yticks(np.linspace(0,20,num=11))

    # ax.set_xlim(0.003,300)
    ax.set_xscale('log')
    #Let's get rid of the scientific notation on the x-axis
    fmt = mpl.ticker.StrMethodFormatter("{x:g}")
    ax.xaxis.set_major_formatter(fmt)
    ax.minorticks_off()
    
    for C in Chl_bin_edges:
        ax.plot([C,C],[z_bin_edges[0],z_bin_edges[-1]],color='k',linewidth=1)
    
    ax.set_ylabel('Depth (m)', fontsize=12,)
    ax.set_xlabel('Chl mg m$^{-3}$)', fontsize=12,)
    
    cbarax = gcf().add_axes((1.02, 0., 0.015, 0.5)) # [left, bottom, width, height] 
    cbar = colorbar(img, cbarax)
    cbar.set_label('R (days$^{-1}$)', fontsize=12)
    
    plt.savefig(fname+'_r_summary.jpg',dpi=500,bbox_inches = 'tight')
    
    
if __name__ == "__main__":
    
    # (could loop through different params settings if we had lots to plot)
    react_params=params.reactions01 # select which reactions params to use
    physics_params=params.physics02 # select which physics params to use
    
    ps= [1e-3,1e-7] #[1e-3,1e-7]
    for p in ps:
        do_the_plot(react_params,physics_params,p)
    
    
