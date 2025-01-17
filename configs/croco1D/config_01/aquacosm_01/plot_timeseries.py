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

react_params = params.reactions01
Qswmaxs=[800]
mean_taus=[0]
mld = 10
#reactions = [Sverdrup_incl_K]
#tau_mean = 0
def eps_daily_ave(kappa,mld):
    kappa_avg= np.zeros(22)
    mld_avg = np.zeros(22)
    for d in arange(1,23):
        indx_after = d*48
        indx_before = indx_after-48
        kappa_avg[d-1] = mean(kappa[indx_before:indx_after])
        mld_avg[d-1] = mean(mld[indx_before:indx_after])

    return kappa_avg,mld_avg
    

for mean_tau in mean_taus:
    for Qswmax in Qswmaxs:
            fig, ax = plt.subplots(figsize=(10,20))
    
            ax = [subplot(3,1,i+1) for i in range(3)]
            
            amplitudes = [0.01,0.02,0.03]
            for amplitude in amplitudes:
                crocodir='../physics/'
                crocofilename = 'mean'+str(mean_tau)+'_mld10_amp'+str(amplitude)+'_flx'+str(Qswmax)+'_lat30_T016_hmax50.nc'
                crocofile=crocodir+crocofilename 
                
                time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
                z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco)
                # get kappa onto the rho axis so each value will have a cell depth
                kappa_croco_r=w2rho(time_croco,zw,z,kappa_croco)
                
                # get the eulerian data
                eulfile='eulerian_'+react_params.Name+'_l' + str(react_params.LightDecay)+ '_K' + str(react_params.CarryingCapacity)+'_r'+str(react_params.BasePhotoRate)+ '_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'
                time_eul,z_eul,chl_eul=get_eul_output(eulfile)
                # interpolate z_therm onto eulerian time axis (just in case different)
                z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
                
                wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
                dt = 5
                React = set_up_reaction(wc, dt, Sverdrup_incl_K, 
                                            LightDecay = react_params.LightDecay,
                                            BasePhotoRate = react_params.BasePhotoRate,
                                            RespirationRate = react_params.RespirationRate,
                                            CarryingCapacity = react_params.CarryingCapacity)
                
                # interpolate kappa_croco_r onto eulerian time axis
                kappa_eul_r=np.zeros_like(chl_eul)
                for ii,zi in enumerate(z):
                    kappa_eul_r[:,ii]=interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([kappa_croco_r[0,ii]],kappa_croco_r[:,ii])),kind='linear')([time_eul])
            
                  
                kappa_eul_euphotic = get_Cs_eulerian(time_eul,z,zw,kappa_eul_r,z_therm_croco)
                #
                p = 1e-7
                aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_'+ react_params.Name+'_l' + str(react_params.LightDecay)+ '_K' + str(react_params.CarryingCapacity)+'_r'+str(react_params.BasePhotoRate)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'  
                z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco)),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul]))
                time_aqc,z_aqc,rank_aqc,chl_aqc,_ = get_aqc_output(aqcfile)
                r = get_aqc_reactions(time_eul*86400, z_aqc, rank_aqc, chl_aqc, React)
                # compute epsilon
                sumr = 0
                raux = np.zeros(1056)
                for i in range(1056):
                    raux[i] = sum(r[i,:])/200
                    i+= 1
                r = raux
                kappa_avg,therm_avg = eps_daily_ave(kappa_eul_euphotic,z_therm_eul)
                eps = (react_params.BasePhotoRate/86400)*(therm_avg**2)/(kappa_avg)
                ax[0].plot(np.linspace(1, 23,22),kappa_avg,label = str(amplitude))
                ax[1].plot(np.linspace(1,23,22),therm_avg,label = str(amplitude))
                ax[2].plot(np.linspace(1,23,22),eps,label = str(amplitude))
                #eps using non-averaged values
                #eps=(react_params.BasePhotoRate/86400)*(z_therm_eul**2)/kappa_eul_euphotic
                
                #ax[0].plot(time_croco[1:],kappa_eul_euphotic,label=str(amplitude))
                #ax[1].plot(time_croco[1:],z_therm_eul,label=str(amplitude))
                #ax[2].plot(time_croco[1:],eps,label=str(amplitude))
            
            for n in range(3):
                ax[n].set_xlim(0,22)
                ax[n].set_xticks(range(0,23))
                #ax[n].set_xticklabels([])
            
            ax[0].set_ylabel("$\kappa_s$ (m$^2$ s$^{-1}$)", fontsize=15)
            ax[0].set_ylim(0, 0.0012)
            
            ax[1].set_ylabel("$\ell$ (m)", fontsize=15)
            ax[1].set_ylim(14, 25)
            
            ax[2].set_ylabel("$\epsilon$ (-)", fontsize=15)
            ax[2].set_xlabel("Time (days)", fontsize=15)
            ax[2].set_ylim(0, 160)
            
            ax[0].legend(fontsize=10, loc="upper left",title="$\\tau^{ac0}$ (N m$^{-2}$)")
                
            plt.savefig('plot_timeseries_'+react_params.Name+'_mean'+str(mean_tau)+'_mld10_flx'+str(Qswmax)+'.jpg',dpi=500,bbox_inches = 'tight')
