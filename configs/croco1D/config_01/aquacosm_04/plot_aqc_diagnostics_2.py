import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
import xarray as xr
from pathlib import Path
import params
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
ion()
react_params = params.reactions()

def plot_generic(ax,n,time,z,data,label,min_val,max_val,cmap):
    # set the position
    if n<2:
        vrt_loc=0.55-0.55*n
        hzt_loc=0    
    else:
        vrt_loc=0.55-0.55*(n-2)
        hzt_loc=1
    ax[n].set_position(  (hzt_loc, vrt_loc, 0.8, 0.5))
    
    img = ax[n].scatter(time, z, c=data,
                        cmap=cmap, linewidth=0, s=40)
    img.set_clim(min_val, max_val)
    ax[n].set_ylim(50, 0)
    ax[n].set_xlim(0,22)
    ax[n].set_xticks(range(0,23))
    ax[n].set_ylabel('Depth (m)', fontsize=15,)
    # add a colorbar
    cbarax = gcf().add_axes((hzt_loc+0.82, vrt_loc, 0.015, 0.5)) # [left, bottom, width, height] 
    cbar = colorbar(img, cbarax)
    cbar.set_label(label, fontsize=15)
    
    if n<4:
        ax[n].set_xticklabels([])
    else:
        ax[n].set_xlabel('Time (days)', fontsize=15)
        
def do_the_plot(mld,amplitude,mean_tau,Qswmax,p):
        
    fname='aquacosm_p'+"{0:1.0e}".format(p)+'_r'+str(react_params.MaxPhotoRate)+'_b'+str(react_params.BasalMetabolism)+'_c'+str(react_params.Chl_light_abs)+'_a'+str(react_params.CrowdingMortality)+'_l'+str(react_params.LightDecay)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)
    diagfile=fname+'_diags.nc'
    print('\n working on ' + diagfile +'\n')
    
    # get the croco input data
    crocodir='../physics/'
    crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename
    time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
    #
    z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco)
    Nt_croco,Nz_croco=shape(temp_croco)
    # repeat time along the z dimension for plotting
    time_croco = np.repeat(time_croco[:, np.newaxis], Nz_croco, axis=1)
    # repeat depth along the time dimension for plotting
    z_croco = np.repeat(z[np.newaxis,:], Nt_croco, axis=0)
    
    # get the aquacosm data
    aqcfile=fname+'.nc'
    time_aqc,z_aqc,rank_aqc,chl_aqc,_ = get_aqc_output(aqcfile)
    
    # get the aquacosm diagnostics
    time_aqc,z_diag,z_aqc,rank_aqc,r,stdev_r,stdev_r_norm,stdev_chl,stdev_chl_norm=get_aqc_diags(diagfile)
    
    # repeat time along the z dimension for plotting
    Nt,Nz_aqc=shape(r)
    Nz_diag=len(z_diag)
    time_aqc_plt = np.repeat(time_aqc[:, np.newaxis], Nz_aqc, axis=1)
    time_diag_plt = np.repeat(time_aqc[:, np.newaxis], Nz_diag, axis=1)
    # repeat depth along the time dimension for plotting
    z_diag_plt = np.repeat(z_diag[np.newaxis,:], Nt, axis=0)
    
    figure(figsize=(10,5))
    ax = [subplot(5,1,i+1) for i in range(5)] # locations are set in plot_generic() - just need the right number of subplots here (excluding the timeseries)
    
    max_chl=12
    plot_generic(ax,0,time_aqc_plt,z_aqc,chl_aqc,'Chl mg m$^{-3}$)',0,max_chl,cm.viridis)
    plot_generic(ax,1,time_diag_plt,z_diag_plt,stdev_chl,'$\sigma_{Chl}$ (mg m$^{-3}$)',0,8,cm.viridis)
    plot_generic(ax,2,time_aqc_plt,z_aqc,r,'R (days$^{-1}$)',-0.8,0.8,cm.RdBu_r)
    plot_generic(ax,3,time_diag_plt,z_diag_plt,stdev_r,'$\sigma_{R}$ (days$^{-1}$)',0,0.25,cm.viridis)
    plot_generic(ax,4,time_diag_plt,z_diag_plt,stdev_r_norm,'$\sigma_{R}/\overline{R}$ (-)',-1,1,cm.RdBu_r)
    # add the thermocline
    for n in range(5):
        cnt=ax[n].contour(time_croco,z_croco,temp_croco,[11],colors='w',linewidths=2.5,linestyles='dashed')
    
    # add a title
    ax[0].text(10, -2, 'p = '+"{0:1.0e}".format(p),fontsize=15)
    
    # get the eulerian data for time-series plot
    eulfile='eulerian_r'+str(react_params.MaxPhotoRate)+'_b'+str(react_params.BasalMetabolism)+'_c'+str(react_params.Chl_light_abs)+'_a'+str(react_params.CrowdingMortality)+'_l'+str(react_params.LightDecay)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
    time_eul,z_eul,chl_eul=get_eul_output(eulfile)
    # interpolate z_therm from croco output onto eulerian time axis (just in case different)
    z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco[:,0])),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
    # get average chl in surface layer for eulerian run
    chl_eul_avg=get_Cs_eulerian(time_eul,z,zw,chl_eul,z_therm_eul)
    # get average chl in surface layer for aquacosm run
    chl_aqc_avg=get_Cs_aquacosm(time_croco[:,0],z_therm_croco,time_aqc,z_aqc,chl_aqc)
    #
    # time-series plot
    ts=gcf().add_axes((0., 0.55-0.55*2, 0.8, 0.5))
    ts.plot(time_eul,chl_eul_avg, 'k', linewidth=2, label='Eulerian')
    ts.plot(time_aqc,chl_aqc_avg, linewidth=2, label='Aquacosms, p = '+"{0:1.0e}".format(p))
    ts.set_ylabel('average surface Chl (mg m$^{-3}$)', fontsize=15)
    ts.set_xlabel('Time (days)', fontsize=15)
    ts.set_ylim(0,max_chl)
    ts.set_xlim(0,22)
    ts.set_xticks(range(0,23))
    ts.grid(linestyle=':', linewidth=0.5)
    ts.legend(fontsize=12, loc="upper left")
    
    plt.savefig(fname+'_diags_2.jpg',dpi=500,bbox_inches = 'tight')
    
    
if __name__ == "__main__":
    
    amplitudes = [0.03] #[0, 0.01, 0.02, 0.03, 0.04]
    mlds = [10] #[10, 25]
    mean_taus = [0] #[0, 0.05]
    Qswmaxs = [250] #[0, 250, 800]
    ps= [1e-3,1e-7] #[1e-3,1e-7]
    for amplitude in amplitudes:
        for mld in mlds:
            for Qswmax in Qswmaxs:
                for mean_tau in mean_taus:
                    crocodir='../physics/'
                    crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
                    crocofile=crocodir+crocofilename
                
                    
                    for p in ps:
                        do_the_plot(mld,amplitude,mean_tau,Qswmax,p)
    
    
