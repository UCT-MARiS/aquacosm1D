import sys
sys.path.insert(0, '../../../../aquacosm1D_lib')
from aquacosm1D import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pathlib import Path
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from scipy.interpolate import interp1d
import math
import matplotlib.pyplot as plt
ion()
from plot_croco_lib import *

Qswmax=800
tau_mean=0.05
amplitude = 0

fig, ax = plt.subplots(figsize=(10,20))
crocofile = 'mean'+str(tau_mean)+'_mld10_amp'+str(amplitude)+'_flx'+str(Qswmax)+'_lat30_T016_hmax50.nc'
time_croco, z, zw, temp_croco, tpas_croco, kappa_croco, u_croco, v_croco, s_flux, tau_x, tau_y, dzetadx = get_croco_output(crocofile)
z_therm_croco=get_z_therm_croco(time_croco,z,temp_croco)
# get kappa onto the rho axis so each value will have a cell depth
kappa_croco_r=w2rho(time_croco,zw,z,kappa_croco)
kappa_croco_r[:,:] = np.flip(kappa_croco_r,1)
kappa_croco_surf=get_Cs_eulerian(time_croco,z,zw,kappa_croco_r,z_therm_croco)
#ax[1].plot(time_croco,kappa_croco_surf,label=str(amplitude))
bx = gcf().add_axes(ax)
kappa, = bx.plot(z,kappa_croco_r[0,:])
bx_title='Time = 0'
def animate(i):
    tindx = (np.abs(time_croco[:] - i)).argmin()
    kappa.set_ydata(kappa_croco_r[tindx,:])
    # bx_title.set_text('Time = '+str(i)+' days')
    bx.set_title('Time = {0:.2f} days'.format(i))

anim = FuncAnimation(fig, animate, frames=np.arange(0,22,6/24))

my_writer=PillowWriter(fps=10) 
anim.save('test.gif', writer=my_writer)