#import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal



def reflec_en(a0):
#a0=np.linspace(1,140,140)
    beta_s=(0.006*a0)**0.5/(1+(0.006*a0)**0.5)
#    beta_s=0.38
#    beta_s=0.05*a0**0.5
    gg_s=1/(1-beta_s**2)**0.5
    phi=0.7*a0/1836
    B=phi+1/gg_s
    print(B)
    p=(beta_s*B+(B**2+beta_s**2-1)**0.5)/(1-beta_s**2)
#    p=(beta_s*B)/(1-beta_s**2)
    ek=1836*0.51*((p**2.0+1.)**0.5-1)
    return ek

if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
  ######## Constant defined here ########
  pi        =     3.1415926535897932384626
  q0        =     1.602176565e-19 # C
  m0        =     9.10938291e-31  # kg
  v0        =     2.99792458e8    # m/s^2
  kb        =     1.3806488e-23   # J/K
  mu0       =     4.0e-7*np.pi       # N/A^2
  epsilon0  =     8.8541878176203899e-12 # F/m
  h_planck  =     6.62606957e-34  # J s
  wavelength=     0.8e-6
  frequency =     v0*2*pi/wavelength
  
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 25,  
          }  
  font_size = 25
##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####



  to_path='./'
  set_relativistic =0 

  a0 = np.linspace(1,150,150)
  ne = a0/np.pi/(0.3e-6/wavelength)
  n_proton = ne*2./8.
  n_carbon = ne*1./8.
  
  PP        = a0**2/(n_proton*1836.*1+n_carbon*1836.*12)
  PP_no_c   = a0**2/(n_proton*1836.*1)
  beta_HB   = PP**0.5/(1+PP**0.5)
  ek_proton = 2*PP/(1+2*PP**0.5)*0.51*1836. 
  ek_proton_no_c = 1.1*2*PP_no_c/(1+2*PP_no_c**0.5)*0.51*1836. 
  ek_reflec = reflec_en(a0)

  sim_a0         = np.array([10.2, 15.3, 30.6, 61.2, 122.4])
  sim_ek_proton  = np.array([31,   49,   130,  227,  335])
  sim_ek_proton_m  = np.array([22,   37,   110,  205,  315])
  sim_ek_proton_sc = np.array([6.6, 15, 38, 95, 174])
#  hote_t   = np.array([1.07, 0.69, 0.53,  0.43,  0.316,  0.264, 0.17, 0.099])

  plt.plot(a0,ek_proton,':',color='crimson',linewidth=4, label='HB theory',zorder=0)
  plt.plot(a0,ek_proton_no_c,'-',color='crimson',linewidth=4, label='no carbon',zorder=0)
  plt.plot(a0,ek_reflec,'--',color='crimson',linewidth=4, label='Reflected proton',zorder=0)
  plt.scatter(sim_a0,sim_ek_proton_m,c='red',marker='o',s=200, label='PIC proton main', edgecolors='black', linewidth='2.5',alpha=1,zorder=2)
  plt.scatter(sim_a0,sim_ek_proton,c='yellow',marker='o',s=200, label='PIC proton', edgecolors='black', linewidth='2.5',alpha=1,zorder=2)
  plt.scatter(sim_a0,sim_ek_proton_sc,c='dodgerblue',marker='o',s=200, label='proton_sheet_cr', edgecolors='black', linewidth='2.5',alpha=1,zorder=1)

  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$a_0$',fontdict=font)
  plt.ylabel(r'$\varepsilon_p$'+' [MeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle=':', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)
#  plt.xlim(0,5) 
#  plt.ylim(0,1.25) 
  plt.legend(loc='best',fontsize=18,framealpha=0.5)
  plt.subplots_adjust(left=0.16, bottom=0.15, right=0.95, top=0.95,
                wspace=None, hspace=None)
    #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(8.4, 7.0)
  fig.savefig('./theory_HB.png',format='png',dpi=160)
  plt.close("all")
