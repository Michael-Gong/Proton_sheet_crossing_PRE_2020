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
    print('beta_s:',beta_s)
#    beta_s=0.38
#    beta_s=0.05*a0**0.5
    gg_s=1/(1-beta_s**2)**0.5
    print('gg_s:',gg_s)
    phi=3.14*3./8.*a0/1836
    B=phi+1/gg_s
    print(B)
    p=(beta_s*B+(B**2+beta_s**2-1)**0.5)/(1-beta_s**2)
#    p=(beta_s*B)/(1-beta_s**2)
    ek=1836*0.51*((p**2.0+1.)**0.5-1)
    beta_i=2.*beta_s/(1+beta_s**2)
    gg_i=1/(1-beta_i**2)**0.5
    ek_hb = 2*(0.006*a0)/(1+2*(0.006*a0)**0.5)*1836*0.51 # 1836.*0.51*(gg_i-1) 
    ek_2 = 1836*0.51*(0.5*0.006*a0+phi/(1-0.006*a0)+(2*phi)**0.5/(1-0.006*a0)/(0.006*a0)**0.5)
    return ek, ek_hb, ek_2

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
  PP_no_c   = 0.8*a0**2/(n_proton*1836.*1)
  beta_HB   = PP**0.5/(1+PP**0.5)
  ek_proton = 2*PP/(1+2*PP**0.5)*0.51*1836. 
  ek_proton_no_c = 2*PP_no_c/(1+2*PP_no_c**0.5)*0.51*1836. 
  ek_reflec,ek_hb,ek_ref2 = reflec_en(a0)

  sim_a0         = np.array([10.2, 15.3, 30.6, 61.2, 122.4])
  sim_ek_proton  = np.array([31,   49,   130,  227,  335])
#  sim_ek_proton_m  = np.array([22,   37,   110,  205,  315])
  sim_ek_proton_m  = np.array([31,   45.5,   121,  210,  335])
  sim_ek_proton_sc = np.array([6.6, 15, 38, 95, 174])
#  hote_t   = np.array([1.07, 0.69, 0.53,  0.43,  0.316,  0.264, 0.17, 0.099])

  plt.plot(a0,ek_reflec,'--',color='red',linewidth=4, label='Theoretic equation',zorder=0)
  F = 0.0385*a0
  theory_LS = 918.*F**2/2/(1+F)
  plt.plot(a0,theory_LS,'--',color='green',linewidth=4, label='LS model',zorder=0)
#  plt.plot(a0,2e-3*a0**2*30.0,'--',color='green',linewidth=4, label='LS model',zorder=0)
  theory_CE_mlt = 230*1.0+np.zeros_like(a0)
  plt.plot(a0,theory_CE_mlt,'--',color='gray',linewidth=4, label='CE MLT model',zorder=0)
  Te_1 = 0.25*a0*0.51
  Te_2 = 2.5*a0**0.667*0.51
  theory_st_1 = Te_1*(np.exp(7.5)*6.5+1)/(np.exp(7.5)-1)
  theory_st_2 = Te_2*(np.exp(7.5)*6.5+1)/(np.exp(7.5)-1)
  plt.fill_between(a0,theory_st_1,theory_st_2,color='blue',alpha=.25,label='Thermal model',zorder=0)
  plt.plot(a0,0.5*2*a0,'--',color='blue',linewidth=4, label='TNSA model',zorder=0)
#  plt.plot(a0,ek_ref2,'--',color='red',linewidth=4, label='Theoretic reflected proton',zorder=0)
#  plt.plot(a0,ek_hb,'--',color='mediumseagreen',linewidth=4, label='Theoretic HB',zorder=0)
  plt.scatter(sim_a0,sim_ek_proton_m,c='tomato',marker='o',s=300, label='PIC ', edgecolors='black', linewidth='2.5',alpha=1,zorder=2)
#  plt.scatter(sim_a0,sim_ek_proton,c='yellow',marker='o',s=200, label='PIC proton', edgecolors='black', linewidth='2.5',alpha=1,zorder=2)
#  plt.scatter(sim_a0,sim_ek_proton_sc,c='dodgerblue',marker='o',s=200, label='proton_sheet_cr', edgecolors='black', linewidth='2.5',alpha=1,zorder=1)
  plt.xlabel('$a_0$',fontdict=font)
  plt.ylabel(r'$\varepsilon_p$'+' [MeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle=':', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)
  plt.xlim(0,140) 
  plt.ylim(0,500) 
  plt.legend(loc='best',fontsize=18,framealpha=0.5)
  plt.subplots_adjust(left=0.16, bottom=0.17, right=0.95, top=0.95,
                wspace=None, hspace=None)
  fig = plt.gcf()
  fig.set_size_inches(10., 8.0)
  fig.savefig('./wrap_proton_theory.png',format='png',dpi=160)
  plt.close("all")


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
  plt.scatter(sim_a0,sim_Te_fit_pre,c='orange',marker='o',s=200, label='PIC $T_e$ at $t=67$ fs' , edgecolors='black', linewidth='2.5',alpha=1,zorder=1)
  plt.scatter(sim_a0,sim_Te_fit,c='deepskyblue',marker='o',s=200, label='PIC $T_e$ at $t=93$ fs' , edgecolors='black', linewidth='2.5',alpha=1,zorder=1)

  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$a_0$',fontdict=font)
  plt.ylabel(r'$T_e$'+' [MeV]',fontdict=font)
  plt.xticks([0,25,50,75,100,125],fontsize=font_size); 
  plt.yticks([0,10,20,30,40],fontsize=font_size);
  plt.grid(which='major',color='k', linestyle=':', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)
  plt.xlim(0,140) 
  plt.ylim(0,40) 
  plt.legend(loc='best',fontsize=18,framealpha=1)
  plt.subplots_adjust(left=0.16, bottom=0.17, right=0.95, top=0.96,
                wspace=None, hspace=None)
    #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(10., 5.0)
  fig.savefig('./wrap_electron_theory.png',format='png',dpi=160)
  plt.close("all")
