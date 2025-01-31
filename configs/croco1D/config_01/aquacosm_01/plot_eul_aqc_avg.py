import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
import params
ion()

react_params = params.reactions01()

def do_the_plot(crocofile,eulfile,aqcfile):
    
    # get the croco output
    time_croco, z, zw, temp_croco, chl_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
    z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco)
    Cs_croco=get_Cs_eulerian(time_croco,z,zw,chl_croco,z_therm_croco)
    Nt_croco,Nz_croco=shape(temp_croco)
    # repeat time along the z dimension for plotting
    time_croco = np.repeat(time_croco[:, np.newaxis], Nz_croco, axis=1)
    # repeat depth along the time dimension for plotting
    z_croco = np.repeat(z[np.newaxis,:], Nt_croco, axis=0)
    
    # get the eulerian output
    time_eul,z_eul,chl_eul=get_eul_output(eulfile)
    # interpolate z_therm onto eulerian time axis (just in case different)
    z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco[:,0])),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
    Cs_eul=get_Cs_eulerian(time_eul,z,zw,chl_eul,z_therm_eul)
    Nt_eul,Nz_eul=shape(chl_eul)
    # repeat time along the z dimension for plotting
    time_eul = np.repeat(time_eul[:, np.newaxis], Nz_eul, axis=1)
    # repeat depth along the time dimension for plotting
    z_eul = np.repeat(z_eul[np.newaxis,:], Nt_eul, axis=0)
    
    # get the aquacosm output
    time_aqc,z_aqc,z_rank,chl_aqc,p=get_aqc_output(aqcfile)
    Cs_aqc=get_Cs_aquacosm(time_croco[:,0],z_therm_croco,time_aqc,z_aqc,chl_aqc)
    Nt_aqc,Nz_aqc=shape(chl_aqc)
    # repeat time along the z dimension for plotting
    time_aqc = np.repeat(time_aqc[:, np.newaxis], Nz_aqc, axis=1)
    
    if 'BioShading_onlyC' in aqcfile:
        max_chl = 30
    elif 'Sverdrup' in aqcfile:
        max_chl = chl_aqc.max()
    elif 'NoReactions' in aqcfile:
        max_chl = 1
    else:
        max_chl = chl_aqc.max()
    
    def plot_C(ax,n,time,z,data,label,t,cmap,min_value,max_value,label_cbar):
        #ax[n].set_position(  (0.08, 0.86-0.14*n, 0.8, 0.13))
        ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
        img = ax[n].scatter(time, z, c=data,
                            cmap=cmap, linewidth=0, s=40)
        img.set_clim(min_value, max_value)
        ax[n].set_ylim(50, 0)
        ax[n].set_xlim(0,21)
        ax[n].set_xticks(range(0,22))
        ax[n].set_ylabel('Depth (m)', fontsize=15,)
        if n==2:
            ax[n].set_xlabel('Time (days)', fontsize=15,)
        else:
            ax[n].set_xticklabels([])
        
        cbarax = gcf().add_axes((0.82, 0.55-0.55*n, 0.015, 0.5)) # [left, bottom, width, height] 
        cbar = colorbar(img, cbarax)
        cbar.set_label(label_cbar, fontsize=15)
        
        ax[n].text(0.2, 45, label,
                   fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
        ax[n].plot([t, t], [0, 50], '--k', linewidth=1)
    
    def plot_C_OLD(n,time,z,carbon,label,t):
        ax[n].set_position(  (0., 0.55-0.55*n, 0.8, 0.5))
        img = ax[n].scatter(time, z, c=carbon,
                            cmap=cm.viridis, linewidth=0, s=40)
        img.set_clim(0, max_chl)
        ax[n].set_ylim(50, 0)
        ax[n].set_xlim(0,21)
        ax[n].set_xticks(range(0,22))
        ax[n].set_ylabel('Depth (m)', fontsize=15,)
        if n==2:
            ax[n].set_xlabel('Time (days)', fontsize=15,)
            cbarax = gcf().add_axes((0.82, -0.5, 0.015, 1.5)) # [left, bottom, width, height] 
            cbar = colorbar(img, cbarax)
            cbar.set_label("Chlorophyll  (mg m$^{-3}$)", fontsize=15)
        else:
            ax[n].set_xticklabels([])
        
        ax[n].text(0.2, 45, label,
                   fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
        ax[n].plot([t, t], [0, 50], '--k',
                            linewidth=1)
        
    figure(figsize=(10,5))
    #ax = [subplot(4,1,i+1) for i in range(4)]
    ax = [subplot(3,1,i+1) for i in range(3)]   
        
    # time to plot
    t=21 #days
    tindx_croco = (np.abs(time_croco[:,0] - t)).argmin()
    tindx_eul = (np.abs(time_eul[:,0] - t)).argmin()
    tindx_aqc = (np.abs(time_aqc[:,0] - t)).argmin()
    
    # do the scatter plots
    plot_C(ax,0,time_croco,z_croco,u_croco,'CROCO',t,cm.RdBu,-u_croco.max(),u_croco.max(),'u (m s$^{-1}$)')
    # plot_C(0,time_croco,z_croco,chl_croco,'CROCO',t)
    plot_C(ax,1,time_eul,z_eul,chl_eul,'Eulerian offline',t,cm.viridis,0,20,"Chlorophyll  (mg m$^{-3}$)")
    plot_C(ax,2,time_aqc,z_aqc,chl_aqc,'Aquacosms, p = '+p,t,cm.viridis,0,20,"Chlorophyll  (mg m$^{-3}$)")
    # add the thermocline
    for n in range(0,3):
        cnt=ax[n].plot(time_croco,z_therm_croco,'w',linestyle='dashed')
        
    # compare the profiles at a given time
    bx=gcf().add_axes((0.0, -1.3, 0.3, 0.6))
    #
    z_croco_t = z_croco[tindx_croco,:]
    chl_croco_t = chl_croco[tindx_croco,:]
    #bx.plot(chl_croco_t, z_croco_t, '0.65', linewidth=6, label="CROCO")
    #
    z_eul_t = z_eul[tindx_eul,:]
    chl_eul_t = chl_eul[tindx_eul,:]
    bx.plot(chl_eul_t, z_eul_t, '-r', linewidth=3, label="Eulerian offline")
    #
    z_aqc_t = z_aqc[tindx_aqc,:]
    chl_aqc_t = chl_aqc[tindx_aqc,:]
    Particles_for_gaussian=array([z_rank,z_aqc_t,chl_aqc_t]).transpose() # reconstructing the particle format for parsing to the gaussian function
    z_aqc_gaus, chl_aqc_gaus = gaussian_estimate_field(Particles_for_gaussian,2, 50.0, 1.25) # default stddev = 2.5
    bx.plot(chl_aqc_gaus, z_aqc_gaus, '-k', linewidth=3, label="Coarse-grained")
    imgr = bx.scatter(chl_aqc_t, z_aqc_t, c=chl_aqc_t, cmap=cm.viridis, s=15,zorder=3,label="Aquacosms")
    imgr.set_clim(0, max_chl)
    bx.set_ylim(51, -1)
    bx.set_xlim(0-0.025*max_chl,1.025*max_chl)
    bx.set_xlabel('Chlorophyll  (mg m$^{-3}$)', fontsize=15)
    bx.set_ylabel('Depth (m)', fontsize=15)
    bx.text(0,-2,'Time = '+str(t)+' days',ha='left',fontsize=12)
    bx.grid(linestyle=':', linewidth=0.5)
    bx.legend(fontsize=12, loc="lower right") #, bbox_to_anchor=(0.1, 0.75),framealpha=0.99)
    
    # compare Cs, the chlorophyll averaged over the surface layer
    ts=gcf().add_axes((0.4, -1.3, 0.45, 0.6))
    #ts.plot(time_croco[:,0],Cs_croco, '0.65', linewidth=6, label='CROCO')
    ts.plot(time_eul[:,0],Cs_eul, '-r', linewidth=3, label='Eulerian')
    ts.plot(time_aqc[:,0],Cs_aqc, '-k', linewidth=3, label='Aquacosms, p = '+p)
    ts.set_ylabel('average surafce Chl (mg m$^{-3}$)', fontsize=15)
    ts.set_xlabel('Time (days)', fontsize=15)
    if 'BioShading_onlyC' in aqcfile:
        ts.set_ylim(0, 30)
    if 'Sverdrup_incl_K' in aqcfile:
        ts.set_ylim(0,max(Cs_aqc.max(),Cs_eul.max())+0.2)
    if 'NoReactions' in aqcfile:
        ts.set_ylim(0, 1)
    ts.set_xlim(0,22)
   # ts.set_ylim(0,11)
    ts.set_xticks(range(0,22))
    ts.grid(linestyle=':', linewidth=0.5)
    if 'NoReactions' in aqcfile:
        ts.legend(fontsize=12, loc="upper right")
    else:
        ts.legend(fontsize=12, loc="upper left")
        
    plt.savefig('plot_eul_aqc_3_'+Path(aqcfile).stem+'.jpg',dpi=500,bbox_inches = 'tight')
    
if __name__ == "__main__":
    mean_tau= 0
    amplitude = 0.03
    Qswmax = 800
    mld = 10
    crocodir='../physics/'
    crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc"
    crocofile=crocodir+crocofilename 
    
    p=1e-7
   
    aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_'+ react_params.Name+'_l' + str(react_params.LightDecay)+ '_K' + str(react_params.CarryingCapacity)+'_r'+str(react_params.BasePhotoRate)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'   
    eulfile='eulerian_'+react_params.Name+'_l' + str(react_params.LightDecay)+ '_K' + str(react_params.CarryingCapacity)+'_r'+str(react_params.BasePhotoRate)+ '_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
    do_the_plot(crocofile,eulfile,aqcfile)
    
    
