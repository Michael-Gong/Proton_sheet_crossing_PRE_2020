#import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
  
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
  theory_1 = a0*0.51
  theory_2 = 5*a0**0.5*0.51
  theory_3 = 1.67*a0**0.75*0.51
  theory_4 = 2.5*a0**0.667*0.51


  sim_a0         = np.array([10.2, 15.3, 30.6, 61.2, 122.4])
  sim_Te         = np.array([8.7,  15.7, 25,   33,   40])
  sim_Te_max     = np.array([20,   28,   38,   39,   40.5])
  sim_Te_fit     = np.array([4.9,  7.6,  14.3, 22.4, 30.2])
  sim_Te_fit_pre = np.array([1.3,  2.4,  4.7,  7.7,  15.6])
#  hote_t   = np.array([1.07, 0.69, 0.53,  0.43,  0.316,  0.264, 0.17, 0.099])

  plt.plot(a0,0.25*theory_1,'--',color='crimson',linewidth=4, label='$T_e=0.25a_0m_ec^2$',zorder=0)
#  plt.plot(a0,theory_1,'-',color='k',linewidth=4, label='$T_e=a0m_ec^2$',zorder=0)
#  plt.plot(a0,theory_2,'--',color='k',linewidth=4, label='$T_e=5a0^{1/2}m_ec^2$',zorder=0)
#  plt.plot(a0,theory_3,':',color='k',linewidth=4, label='$T_e=1.7a0^{3/4}m_ec^2$',zorder=0)
  plt.plot(a0,theory_4,'--',color='mediumblue',linewidth=4, label='$T_e=2.5a_0^{2/3}m_ec^2$',zorder=0)
#  plt.plot(a0,ek_proton_no_c,'-',color='crimson',linewidth=4, label='no carbon',zorder=0)
#  plt.scatter(sim_a0,sim_Te,c='lime',marker='o',s=200, label='PIC $T_e$', edgecolors='black', linewidth='2.5',alpha=1,zorder=2)
#  plt.scatter(sim_a0,sim_Te_max,c='orange',marker='^',s=200, label='PIC Max{$T_e$}', edgecolors='black', linewidth='2.5',alpha=1,zorder=1)
  plt.scatter(sim_a0,sim_Te_fit_pre,c='orange',marker='o',s=200, label='PIC $T_e$ pre' , edgecolors='black', linewidth='2.5',alpha=1,zorder=1)
  plt.scatter(sim_a0,sim_Te_fit,c='deepskyblue',marker='o',s=200, label='PIC $T_e$ fit' , edgecolors='black', linewidth='2.5',alpha=1,zorder=1)

  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$a_0$',fontdict=font)
  plt.ylabel(r'$T_e$'+' [MeV]',fontdict=font)
  plt.xticks([0,25,50,75,100,125],fontsize=font_size); 
  plt.yticks([0,10,20,30,40],fontsize=font_size);
  plt.grid(which='major',color='k', linestyle=':', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)
  plt.xlim(0,140) 
  plt.ylim(0,40) 
  plt.legend(loc='best',fontsize=18,framealpha=0.5)
  plt.subplots_adjust(left=0.16, bottom=0.15, right=0.95, top=0.95,
                wspace=None, hspace=None)
    #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(8.4, 7.0)
  fig.savefig('./theory_Te.png',format='png',dpi=160)
  plt.close("all")
