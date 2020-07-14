#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage
  
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
  wavelength=     1.0e-6
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

  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          } 

  font_size  = 25
  font_size2 = 20

##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('crimson')
  c_blue= matplotlib.colors.colorConverter.to_rgba('mediumblue')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128) 
##end for transparent colorbar##
 
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

  ######### Parameter you should set ###########
  start   =  50 # start time
  stop    =  700  # end time
  step    =  5  # the interval or step
  from_path='./PW_w020_trace2/'
  to_path  ='./PW_w020_trace2_fig/'

  marker_size=8
  space_2 = 1

  e_xx = np.loadtxt(from_path+'electron_xx.txt') 
  e_yy = np.loadtxt(from_path+'electron_yy.txt') 
  e_px = np.loadtxt(from_path+'electron_xx.txt') 
  e_py = np.loadtxt(from_path+'electron_yy.txt') 
  e_gg = (e_px**2+e_py**2+1.0)**0.5
#  for n in range(start,stop+step,step):
  if 1>0:
    #### header data ####
    n=699
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 

    name = 'ey'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    eee = 60.0
    levels = np.linspace(-eee, eee, 41)
    ex.T[ex.T < -eee]=-eee
    ex.T[ex.T >  eee]= eee
    plt.contourf(X, Y, ex.T, levels=levels, cmap=cmap_br)
    #### manifesting colorbar, changing label and axis properties ####
    cbar=plt.colorbar(pad=0.01,ticks=np.linspace(-eee, eee, 5))
    cbar.ax.tick_params(labelsize=font2['size'])
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
    cbar.set_label('$E_y$ [$m_ec\omega_0/|e|$]',fontdict=font2)

    data = sdf.read(from_path+'trace'+str(n).zfill(4)+".sdf",dict=True)
    ion_x  = data['Grid/Particles/subset_trace_p/t_proton'].data[0]/1e-6
    ion_y  = data['Grid/Particles/subset_trace_p/t_proton'].data[1]/1e-6
    ion_px = data['Particles/Px/subset_trace_p/t_proton'].data/(1836.*m0*v0)
    ion_py = data['Particles/Py/subset_trace_p/t_proton'].data/(1836.*m0*v0)
    ion_gg = (ion_px**2+ion_py**2+1.)**0.5
    ion_ek = (ion_gg-1)*0.51*1836.
#          ion_px = data['Particles/Px/proton'].data/(1836*m0*v0)
    plt.scatter(ion_x[::space_2], ion_y[::space_2], c=ion_ek[::space_2], norm=colors.Normalize(vmin=0, vmax=100), cmap='cividis', s=marker_size, marker='.',alpha=0.8,label='proton',zorder=3,lw=0)
    cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 100, 5))
    cbar.ax.tick_params(labelsize=font2['size'])
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
    cbar.set_label(r'$\varepsilon_i$'+' [MeV]',fontdict=font2)
    
    #cbar.set_label('Normalized electric field',fontdict=font)        
    if n-0 >= 100:
        print('n-0=',n-0)
        index=np.arange(n-0-100,n-0,1)
        print('index is',index)
        select = np.where(e_gg[:,n-0] > 1.0)
        print('select is',select)
        print('e_xx[:,index][select,:] shape is',e_xx[:,index][select,:].shape)
        #print('yy_2d_x[select,index] shape is',(yy_2d_x[select,index]).shape)
        #print('xx_2d_x[select,0].size is',xx_2d_x[select,0].size)
        plt.scatter(e_xx[:,index][select,:].reshape(np.array(select).size,index.size), e_yy[:,index][select,:].reshape(np.array(select).size,index.size), c=np.tile(np.arange(index.size)*5, (e_xx[select,0].size, 1)), s=np.tile(np.arange(index.size)*5, (e_xx[select,0].size, 1))*0.005, cmap='winter', edgecolors='None',zorder=4)

    plt.ylim(-8,8)
    plt.xlim(-4,16)
    plt.xlabel('$x$ [$\mu m$]',fontdict=font)
    plt.ylabel('$y$ [$\mu m$]',fontdict=font)
    plt.xticks([0,5,10,15],fontsize=font_size); 
    plt.yticks([-8,-4,0,4,8],fontsize=font_size);
    plt.title('At '+str(round(time/1.0e-15,4))+' fs',fontdict=font)
   
    fig = plt.gcf()
    fig.set_size_inches(14, 8)
    fig.savefig('wrap_shining.png',format='png',dpi=160)
    plt.close("all")
    print('wrap_shining.png')

