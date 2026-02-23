from pylab import *

#--- Must match the parameters in co_existence_1000.py ---
dt        = 5. # time step in seconds
Ndays     = 100 #4*365 #length of the simulation
Nloops    = int(24*3600  *  Ndays  / dt)
Nstore    = int(1*3600 / dt) #store the particles every Nshow time steps
Npts      = 200  #number of particles
Nspecies  = 1000 #number of distinct plankton species
Nscalars  = Nspecies+1   #number of scalars carried by each particle
p         = 1.e-5

days = arange(0, Nloops, Nstore)*dt/86400 #time index of Pstore in days
#---------------------------------------------------------

basename = f"Species_{Nspecies}.p{p:1.0e}.Gompertz_Ktop_1.e-1.Kbt_1.e-4"
print("Loading")
AlphaEpsilon    = fromfile(basename+".AlphaEpsilon.pydat")
MaxPhotoRate    = fromfile(basename+".MaxPhotoRate.pydat")
BasalMetabolism = fromfile(basename+".BasalMetabolism.pydat")
HalfSaturation  = fromfile(basename+".HalfSaturation.pydat")
Pstore          = fromfile(basename+".Pstore.pydat")
Pstore = reshape(Pstore, (Nloops//Nstore, Npts, Nspecies+3))

print("Loading done")

#-----------------------------------------------------------

semilogy(days,  mean(Pstore[:,:,-1], axis=1), '-y',
         linewidth=3, label="Nutrient")
for i in range(Nspecies):
    semilogy(days,  mean(Pstore[:,:,2+i], axis=1), '-k',
             linewidth=0.25)
xlabel("Days", fontsize=18)
ylabel("Concentration", fontsize=18)
legend(['Nutrient', 'C content of the i-th species'], loc="upper left")
tight_layout()

#-----------------------------------------------------------
#-- Find which species reach the highest concentration

SpeciesConc = [] 
# Each row contains (species n; last conc; conc change over last day, max/rms concentration; std of concentration)
for i in range(2, Nspecies+2):
    SpeciesConc.append([
        i, mean(Pstore[-1,:,i]),
        (mean(Pstore[-1,:,i]) - mean(Pstore[-1-int(86400/(dt*Nstore)),:,i])),
        amax(Pstore[-1,:,i])/sqrt(sum(Pstore[-1,:,i]**2)),
        std(Pstore[-1,:,i])/mean(Pstore[-1,:,i])
    ])
SpeciesConc = array(SpeciesConc)
idx = argsort(SpeciesConc[:,1])[::-1]

print(f"Concentration    Daily con. change    Max/L2       c.var.      alpha*1e5        r           mu           K")
for i in idx[:30]:
    print(f"{SpeciesConc[i,1]:>10.4f}      {SpeciesConc[i,2]:>14.6f}    {SpeciesConc[i,3]:10.4f}    {SpeciesConc[i,4]:10.6f}    {AlphaEpsilon[i]*1.e5*0.217/0.4:>9.6f}    {MaxPhotoRate[i]:>9.6f}    {BasalMetabolism[i]:>9.6f}    {HalfSaturation[i]:>9.6f}")
