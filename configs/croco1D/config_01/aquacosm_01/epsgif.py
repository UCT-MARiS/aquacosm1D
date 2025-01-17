import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from scipy.interpolate import interp1d
from plot_eul_aqc_lib import *
import math
import matplotlib.pyplot as plt
import params
ion()

react_params = params.reactions01()

mean_taus = [0]#[0,0.05]#[0,0.05]
amplitudes = [0.01] 
Qswmax = 800#[250,800]
mld = 10
Ls = [2.5,5.]#[5.,2.5]
Ks = [10,20]#[10,20]
growth =[0.5,1]# [0.5,1]

def eps_daily_ave(kappa,mld,eul,lag):
    kappa_avg= np.zeros(22)
    mld_avg = np.zeros(22)
    eul_avg = np.zeros(22)
    lag_avg = np.zeros(22)
    for d in arange(1,23):
        indx_after = d*48
        indx_before = indx_after-48
        kappa_avg[d-1] = mean(kappa[indx_before:indx_after])
        mld_avg[d-1] = mean(mld[indx_before:indx_after])
        eul_avg[d-1] = mean(eul[indx_before:indx_after])
        lag_avg[d-1] = mean(lag[indx_before:indx_after])
    diff_avg = abs(eul_avg - lag_avg)
    return kappa_avg,mld_avg,diff_avg
    
def get_kappa_surface(time,z,zw,kappa,z_therm):
    kappa_s=np.zeros(len(time))  
    
    z_thickness=zw[1:]-zw[0:-1];
    
    for t in range(0,len(time)):
        #could put a check here to make sure temp_thermocline is in the range of temp[t,:]
        # find the z indices above and below the interpolated depth
        indx_above=np.where(zw<z_therm[t])[-1][-1]
        indx_below=indx_above+1
        # compute the weights to apply to the indices either side of the depth of max strat
        weight_above=(z[indx_above]-z_therm[t])/(z[indx_above]-z[indx_below]);
        weight_below=1-weight_above;
        #
        # integrate kappa from the surface to z_therm
        for k in range(0,indx_above): # using range means we don't include indx_above in the loop, which is what we want
            kappa_s[t]=kappa_s[t] + kappa[t,k]*z_thickness[k]
        # add half of the layer above
        kappa_s[t]=kappa_s[t]+0.5*kappa[t,indx_above]*z_thickness[indx_above]
        # add the mean of the layer above and below, weighted by where the
        # thermocline is relative to the two layers
        kappa_s[t]=kappa_s[t]+weight_above*np.mean(kappa[t,indx_above:indx_below])*np.mean(z_thickness[indx_above:indx_below])
    
    kappa_s[kappa_s == 0] = np.NaN 
    z_therm[z_therm == 0] = np.NaN   

    kappa_s = kappa_s/z_therm # now the mean over the surface layer    
    return kappa_s

def do_the_anim():
    t = 0
    fig = figure(figsize=(10,5))
    bx = gcf().add_axes((0.08, 0.05, 0.9, 0.9))
    #bx.set_ylim(0,0.4)
    #bx.set_xlim(0,35)
    n = len(amplitudes)*len(growth)*len(mean_taus)*len(Ls)*len(Ks)
    eps,diff,time_croco,label = get_eps()
    print(label)
    eps_scat = bx.scatter(eps[t,:],diff[t,:],c=np.arange(0,n),cmap='viridis')
    bx.legend(labels=label,loc='best')
    bx.set_ylabel('E-L')
    bx.set_xlabel('$\epsilon$')
    bx.set_xlim(0,5)
    bx.set_ylim(0,0.35)
    def animate(i):
        f = int(i)
       # f = int(i/22*1056)
        eps_scat.set_offsets(np.c_[eps[f,:],diff[f,:]])
        #bx.set_xlim(min(eps[f,:]-0.5),max(eps[f,:])+0.5)
        #bx.set_ylim(0,max(diff[f,:]+0.05))
        bx.set_title('Time = {0:.2f} days'.format(i))
    anim = FuncAnimation(fig, animate, frames=np.arange(0,22,1))
    my_writer=PillowWriter(fps=10) 
    anim.save('eps.gif', writer=my_writer)
   
    
    
    
def get_eps():
    # get the croco output
    n = len(amplitudes)*len(growth)*len(mean_taus)*len(Ls)*len(Ks)
    eps = np.zeros((22,n))
    diff = np.zeros((22,n))
    labels = np.chararray(n,itemsize=20)
    j = 0
    crocodir='../physics/' 
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
                        kappa_croco_surf=get_kappa_surface(time_croco,z,zw,kappa_croco,z_therm_croco) 
                        #kappa_croco_surf = kappa_croco_surf[1:]
                        #z_therm_croco = z_therm_croco[1:]
                        rrate = get_aqc_reactions(time_eul*86400, z_aqc, z_rank, chl_aqc, React)
                        # compute epsilon
                        sumr = 0
                        raux = np.zeros(1056)
                        for i in range(1056):
                            raux[i] = sum(rrate[i,:])/200
                            i+= 1
                        rrate = raux
                        time_eul,z_eul,chl_eul=get_eul_output(eulfile) 
                        z_therm_eul = np.squeeze(interp1d(concatenate(([time_eul[0]], time_croco[:])),concatenate(([z_therm_croco[0]], z_therm_croco)),kind='linear')([time_eul])) 
                        Cs_eul=get_Cs_eulerian(time_eul,z,zw,chl_eul,z_therm_eul) 
                        time_aqc,z_aqc,z_rank,chl_aqc,p=get_aqc_output(aqcfile) 
                        Cs_aqc=get_Cs_aquacosm(time_croco[:],z_therm_croco,time_aqc,z_aqc,chl_aqc) 
                        #tindx_croco = (np.abs(time_croco[:] - t)).argmin() 
                        kappa_croco_surf=get_kappa_surface(time_croco,z,zw,kappa_croco,z_therm_croco) 
                        kappa_croco_surf = kappa_croco_surf[1:]
                        z_therm_croco = z_therm_croco[1:]
                        kappa_avg,therm_avg,diff[:,j] = eps_daily_ave(kappa_croco_surf, z_therm_eul, Cs_eul, Cs_aqc)
                        eps[:,j] = (r/86400)*(therm_avg**2)/(kappa_avg)
                        diff[:,j] = diff[:,j]/(max(np.max(Cs_eul),np.max(Cs_aqc)))
                        #eps[:,j] = (r/86400)*(z_therm_eul**2)/(kappa_croco_surf)
                        #labels[j] = 'A'+str(amplitude)+'l'+str(L)+'r'+str(r)+'tau'+str(mean_tau)
                        #diff[:,j] = (np.abs(Cs_eul-Cs_aqc))/(max(np.max(Cs_eul),np.max(Cs_aqc))) ###normalize diff
                        j += 1
    return np.log(eps),diff,time_croco,labels
if __name__ == "__main__":
    do_the_anim()
