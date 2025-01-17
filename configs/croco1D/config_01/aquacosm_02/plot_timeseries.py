import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
import matplotlib.pyplot as plt
import params
ion()
react_params = params.reactions()

Qswmaxs=[800]
mean_taus=[0]

for mean_tau in mean_taus:
    for Qswmax in Qswmaxs:
        
        fig, ax = plt.subplots(figsize=(10,20))

        #ax = [subplot(3,1,i+1) for i in range(3)]
    
        
        # amplitudes = [0.01, 0.03]
        amplitudes = [0.03]
        mld = 10
        for amplitude in amplitudes:
            crocodir='../physics/'
            crocofilename = 'mean'+str(mean_tau)+'_mld10_amp'+str(amplitude)+'_flx'+str(Qswmax)+'_lat30_T016_hmax50.nc'
            crocofile=crocodir+crocofilename 
            
            dt = 5             
            wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
            React = set_up_reaction(wc, dt, BioShading_onlyC,
                                    LightDecay=react_params.LightDecay,
                                    MaxPhotoRate = react_params.MaxPhotoRate, 
                                    BasalMetabolism = react_params.BasalMetabolism,
                                    Chl_C = react_params.Chl_C,
                                    CrowdingMortality = react_params.CrowdingMortality,
                                    Chl_light_abs = react_params.Chl_light_abs)
            
            time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
            z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco)
            # get kappa onto the rho axis so each value will have a cell depth
            kappa_croco_r=w2rho(time_croco,zw,z,kappa_croco)
            
            # get the eulerian data
            eulfile='eulerian_r'+str(react_params.MaxPhotoRate) +'_b'+str(react_params.BasalMetabolism)+'_c'+str(react_params.Chl_light_abs)+'_a'+str(react_params.CrowdingMortality)+'_l'+str(react_params.LightDecay)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
            time_eul,z_eul,chl_eul=get_eul_output(eulfile)
            # interpolate z_therm onto eulerian time axis (just in case different)
            p = 1e-7
            aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_r'+str(react_params.MaxPhotoRate)+'_b'+str(react_params.BasalMetabolism)+'_c'+str(react_params.Chl_light_abs)+'_a'+str(react_params.CrowdingMortality)+'_l'+str(react_params.LightDecay)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
            z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
            time_aqc,z_aqc,rank_aqc,chl_aqc,_ = get_aqc_output(aqcfile)
            # interpolate kappa_croco_r onto eulerian time axis
            kappa_eul_r=np.zeros_like(chl_eul)
            for ii,zi in enumerate(z):
                kappa_eul_r[:,ii]=interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([kappa_croco_r[0,ii]],kappa_croco_r[:,ii])),kind='linear')([time_eul])
            # interpolate kappa_croco_r onto eulerian time axis
            kappa_eul_r=np.zeros_like(chl_eul)
            for ii,zi in enumerate(z):
                kappa_eul_r[:,ii]=interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([kappa_croco_r[0,ii]],kappa_croco_r[:,ii])),kind='linear')([time_eul])
        
              
            kappa_eul_euphotic = get_Cs_eulerian(time_eul,z,zw,kappa_eul_r,z_therm_croco)
            #
            r = get_aqc_reactions(time_eul*86400, z_aqc, rank_aqc, chl_aqc, React)
            # compute epsilon
            sumr = 0
            raux = np.zeros(1008)
            for i in range(1008):
                raux[i] = sum(r[i,:])/200
                i+= 1
            r = raux/(86400*(time_eul[1]-time_eul[0]))
            #print(r,z_therm_eul,kappa_eul_euphotic)
            eps=(r)*np.square(z_therm_eul)/kappa_eul_euphotic
            
            #ax[0].plot(time_eul,kappa_eul_euphotic,label=str(amplitude))
            #ax[1].plot(time_eul,z_therm_eul,label=str(amplitude))
            #ax.plot(time_eul,eps,label=str(amplitude))
        
       # for n in range(3):
        #    ax[n].set_xlim(1,21)
         #   ax[n].set_xticks(range(0,21))
        
        
        # sx=gcf().add_axes((0.0, 1.1, 0.8, 0.4))
        # sx.plot(time_croco,s_flux,'k')
        # sx.set_xticklabels([])
        # sx.set_ylabel('Surface radiation (W m$^{-2}$)', fontsize=15,)
        # sx.set_ylim(0,800)
        
        #ax[0].set_ylabel("$\kappa_s$ (m$^2$ s$^{-1}$)", fontsize=15)
        #ax[0].set_xticklabels([])
        #ax[0].set_ylim(0, 0.0025)
        
        #ax[1].set_ylabel("$\ell$ (m)", fontsize=15)
        #ax[1].set_xticklabels([])
        #ax[1].set_ylim(10, 25)
        
        #ax[2].set_ylabel("$\epsilon$ (-)", fontsize=15)
        #ax[2].set_xlabel("Time (days)", fontsize=15)
        #ax.set_ylim(-1, 8)
        #ax.set_yticks(range(-1,8,90))
        #ax.set_xticks(range(0,22,44))
        plt.plot(time_eul,eps,'.')
        #plt.ylim(-3,8)
        #plt.xticks(np.arange(0, 22, 0.5))
        #plt.yticks(np.arange(-3,8,0.1))
        plt.show()
        
        # ax[0].legend(fontsize=10, loc="upper left",title="$\\tau^{ac0}$ (N m$^{-2}$)")
            
        plt.savefig('plot_timeseries_r'+str(React.MaxPhotoRate*(60.*60.*24.))+'_c'+str(React.Chl_light_abs)+'_a'+str(React.CrowdingMortality*(60.*60.*24.))+'_l'+str(React.LightDecay)+'_mean'+str(mean_tau)+'_mld10_flx'+str(Qswmax)+'.jpg',dpi=500,bbox_inches = 'tight')
