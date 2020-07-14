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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import multiprocessing as mp
  
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
  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 16,  
          } 
  space_1 = 4
  space_2 = 4
  font_size = 25
  marker_size=0.05
##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('red')
  c_blue= matplotlib.colors.colorConverter.to_rgba('blue')
  c_yellow= matplotlib.colors.colorConverter.to_rgba('yellow')
  c_cyan= matplotlib.colors.colorConverter.to_rgba('cyan')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128)
  cmap_ryb= matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_yellow,c_blue],128)
  cmap_yc = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_yellow,c_white_trans,c_cyan],128) 

  c_brown = matplotlib.colors.colorConverter.to_rgba('gold')
  c_green = matplotlib.colors.colorConverter.to_rgba('springgreen')
  cmap_bg = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_brown,c_white_trans,c_white_trans,c_green],128)
  c_black = matplotlib.colors.colorConverter.to_rgba('black')
  c_white = matplotlib.colors.colorConverter.to_rgba('white')
  cmap_bw = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_white_trans,c_white,c_black],128)
   
##end for transparent colorbar##
  upper = matplotlib.cm.jet(np.arange(256))
  lower = np.ones((int(256/4),4))
  for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
  cmap = np.vstack(( lower, upper ))
  mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

 
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

  color_list = ['blue','limegreen','red'] 
  
def processplot(n): 
    #from_path = './new_single_a30/'
    #to_path   = './new_single_a30_fig/'
    from_path = './CP_w020/'
    to_path   = './CP_w020_fig/'
    ######### Script code drawing figure ################
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 
    den_e = data['Derived/Number_Density/electron'].data/denunit
    den_p = data['Derived/Number_Density/proton'].data/denunit
    den_c = data['Derived/Number_Density/carbon'].data/denunit

    X=X[::space_1,::space_2]
    Y=Y[::space_1,::space_2]
    den_e=den_e[::space_1,::space_2]
    den_p=den_p[::space_1,::space_2]
    den_c=den_c[::space_1,::space_2]
    
    plt.subplot(3,1,1)
    levels = np.logspace(-2, 2, 41)
#    den_e.T[den_e.T < ]=-eee
    plt.contourf(X, Y, den_e.T, levels=levels, cmap=mycolor_jet, norm=colors.LogNorm(vmin=0.01,vmax=100))
    cbar=plt.colorbar(pad=0.01,ticks=[0.01,0.1,1,10,100])
    cbar.ax.tick_params(labelsize=font2['size']) 
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
    cbar.set_label('$n_e$ [$n_c$]',fontdict=font2)        
    #plt.text(2.,10,str(round(time/1.0e-15,0))+' fs',fontdict=font2,color='k')
    plt.ylim(-12,12)
    plt.xlim(-10,70)
    plt.xlabel('$x$ [$\mu m$]',fontdict=font)
    plt.ylabel('$y$ [$\mu m$]',fontdict=font)
    plt.xticks(fontsize=font_size); 
    plt.yticks(fontsize=font_size);
#    plt.legend(loc='best',fontsize=font_size,framealpha=0.5)
    plt.title('electron number density',fontdict=font)

    plt.subplot(3,1,2)
    levels = np.logspace(-2, 2, 41)
#    den_e.T[den_e.T < ]=-eee
    plt.contourf(X, Y, den_p.T, levels=levels, cmap='gist_ncar_r', norm=colors.LogNorm(vmin=0.01,vmax=100))
    cbar=plt.colorbar(pad=0.01,ticks=[0.01,0.1,1,10,100])
    cbar.ax.tick_params(labelsize=font2['size']) 
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
    cbar.set_label('$n_p$ [$n_c$]',fontdict=font2)        
    #plt.text(2.,10,str(round(time/1.0e-15,0))+' fs',fontdict=font2,color='k')
    plt.ylim(-12,12)
    plt.xlim(-10,70)
    plt.xlabel('$x$ [$\mu m$]',fontdict=font)
    plt.ylabel('$y$ [$\mu m$]',fontdict=font)
    plt.xticks(fontsize=font_size); 
    plt.yticks(fontsize=font_size);
#    plt.legend(loc='best',fontsize=font_size,framealpha=0.5)
    plt.title('proton number density',fontdict=font)

    plt.subplot(3,1,3)
    levels = np.logspace(-2, 2, 41)
#    den_e.T[den_e.T < ]=-eee
    plt.contourf(X, Y, den_c.T, levels=levels, cmap='terrain_r', norm=colors.LogNorm(vmin=0.01,vmax=100))
    cbar=plt.colorbar(pad=0.01,ticks=[0.01,0.1,1,10,100])
    cbar.ax.tick_params(labelsize=font2['size']) 
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
    cbar.set_label('$n_c$ [$n_c$]',fontdict=font2)        
    #plt.text(2.,10,str(round(time/1.0e-15,0))+' fs',fontdict=font2,color='k')
    plt.ylim(-12,12)
    plt.xlim(-10,70)
    plt.xlabel('$x$ [$\mu m$]',fontdict=font)
    plt.ylabel('$y$ [$\mu m$]',fontdict=font)
    plt.xticks(fontsize=font_size); 
    plt.yticks(fontsize=font_size);
#    plt.legend(loc='best',fontsize=font_size,framealpha=0.5)
    plt.title('carbon number density',fontdict=font)

    plt.subplots_adjust(left=0.12, bottom=0.1, right=0.98, top=0.97, wspace=0.011, hspace=0.16)

    fig = plt.gcf()
    fig.set_size_inches(12, 12.)
    fig.savefig(to_path+'den_comb_'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    print(to_path+'den_comb_'+str(n).zfill(4)+'.png')
    return 0

if __name__ == '__main__':
  start   =  1 # start time
  stop    =  19  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=1)
  results = pool.map(processplot,inputs)
  print(results)
