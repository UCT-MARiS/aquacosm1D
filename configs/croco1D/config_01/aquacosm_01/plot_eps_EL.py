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

if __name__ == "__main__":
    mean_taus = [0]#[0,0.05]#[0,0.05]
    amplitudes = [0.02] 
    Qswmax = 800#[250,800]
    mld = 10
    growth = [1]
    Ls = [5.,2.5]
    Ks = [10,20]
    j=0
    crocodir='../physics/' 
    figure(figsize=(10,5))
    #ts=gcf().add_axes((0.4, -1.3, 0.45, 0.6))
    n = len(amplitudes)*len(growth)*len(mean_taus)*len(Ls)*len(Ks)
    tindx = 1055
    
    for amplitude in amplitudes:
        for mean_tau in mean_taus:
            for L in Ls:
                for K in Ks:
                    for r in growth:
                        p=1e-7
                        crocofilename="mean"+str(mean_tau)+"_mld"+str(mld)+"_amp"+str(amplitude)+"_flx"+str(Qswmax)+"_lat30_T016_hmax50.nc" 
                        crocofile=crocodir+crocofilename 
                        aqcfile='aquacosm_p'+"{0:1.0e}".format(p)+'_'+ react_params.Name+'_l' + str(L)+ '_K' + str(K)+'_r'+str(r)+'_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc'    
                        eulfile='eulerian_'+react_params.Name+'_l' + str(L)+ '_K' + str(K)+'_r'+str(r)+ '_mean'+str(mean_tau)+"_amp"+str(amplitude)+"_mld"+str(mld)+"_flx"+str(Qswmax)+'.nc' 
                        time_croco, z, zw, temp_croco, chl_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile) 
                        ime_aqc,z_aqc,rank_aqc,chl_aqc,_ = get_aqc_output(aqcfile)
                        z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco) 
                        wc = water_column_netcdf(DatasetName=crocofile, max_depth=50)
                        dt = 5
                        React = set_up_reaction(wc, dt, Sverdrup_incl_K, 
                                                    LightDecay = L,
                                                    BasePhotoRate = r,
                                                    RespirationRate = react_params.RespirationRate,
                                                    CarryingCapacity = K)
                        time_eul,z_eul,chl_eul=get_eul_output(eulfile) 
                        z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco[:])),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul])) 
                        Cs_eul=get_Cs_eulerian(time_eul,z,zw,chl_eul,z_therm_eul) 
                        time_aqc,z_aqc,z_rank,chl_aqc,p=get_aqc_output(aqcfile) 
                        Cs_aqc=get_Cs_aquacosm(time_croco[:],z_therm_croco,time_aqc,z_aqc,chl_aqc) 
                        #tindx_croco = (np.abs(time_croco[:] - t)).argmin() 
                        kappa_croco_surf=get_Cs_eulerian(time_croco,z,zw,kappa_croco[:,1:],z_therm_croco) 
                        #kappa_croco_surf = kappa_croco_surf[1:]
                        #z_therm_croco = z_therm_croco[1:]
                        rrate = get_aqc_reactions(time_eul*86400, z_aqc, rank_aqc, chl_aqc, React)
                        # compute epsilon
                        sumr = 0
                        raux = np.zeros(1056)
                        for i in range(1056):
                            raux[i] = sum(rrate[i,:])/200
                            i+= 1
                        rrate = raux
                        label = 'A'+str(amplitude)+'l'+str(L)+'r'+str(r)+'tau'+str(mean_tau)+'K'+str(K)
                        print(label)
                        eps= (r/86400)*(z_therm_croco[tindx]**3)/(kappa_croco_surf[tindx]*L)
                        diff = (abs(Cs_eul[tindx]-Cs_aqc[tindx]))/(max(Cs_eul[tindx],Cs_aqc[tindx])) ###normalize diff
                        plt.scatter(np.log(eps),diff,cmap='viridis_r',label = label)
                        j += 1
    
    plt.legend(loc='best')
    plt.ylabel('E-L')
    plt.xlabel('$\epsilon$')
    #ts.set_xlim(-1,2)
    plt.savefig('plot_eps_EL'+'24p.jpg',dpi=500,bbox_inches = 'tight')
    plt.show()