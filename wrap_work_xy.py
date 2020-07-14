import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
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


######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 28,  
        }  

font2 = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 24,  
        }  

font_size =28

upper = matplotlib.cm.viridis(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_viridis = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


if __name__ == '__main__':
  start   =  1 # start time
  stop    =  19  # end time
  step    =  1  # the interval or step
    
  from_path = './PW_w020/'
  to_path   = './PW_w020_fig/'
  part_name = 'electron'  
  mass      = 1.0

  n=9
  data   = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  header = data['Header']
  time1  = header['time']
  px     = data['Particles/Px/'+part_name].data/(mass*m0*v0)
  py     = data['Particles/Py/'+part_name].data/(mass*m0*v0)
  work_x = data['Particles/Time_Integrated_Work_x/'+part_name].data*0.51
  work_y = data['Particles/Time_Integrated_Work_y/'+part_name].data*0.51
  gg     = (px**2+py**2+1)**0.5
  theta  = np.arctan2(py,px)*180.0/np.pi
  grid_x = data['Grid/Particles/'+part_name].data[0]/1.0e-6      
  grid_y = data['Grid/Particles/'+part_name].data[1]/1.0e-6      
  temp_id= data['Particles/ID/'+part_name].data
  weight = data['Particles/Weight/'+part_name].data*4e-6

  condition = (abs(grid_y)<12) & (px>0)    
  work_x = work_x[condition]
  work_y = work_y[condition]
  weight = weight[condition] 
  
  fig,host = plt.subplots()
#      plt.hist2d(work_x, work_y, bins=(200, 200), range=[[-200,800],[-200,800]], cmap='cubehelix_r', weights=weight, normed=False, norm=colors.Normalize(vmin=0, vmax=1e9))
#      plt.hist2d(work_x, work_y, bins=(200, 200), range=[[-200,800],[-200,800]], cmap='cubehelix_r', weights=weight, normed=False)
  plt.hist2d(work_x, work_y, bins=(200, 200), range=[[-700,300],[-300,700]], cmap=mycolor_viridis, weights=weight, normed=False,norm=colors.LogNorm(vmin=1e5,vmax=1e11))
  cbar=plt.colorbar(pad=0.005,ticks=[1e5,1e7,1e9,1e11])
  cbar.set_label('dN/(dW$_x$dW$_y$)'+' [A.U.]',fontdict=font)
  cbar.ax.tick_params(labelsize=font2['size'])
  plt.xlim(-700,300)
  plt.ylim(-300,700)
  plt.xlabel('W$_x$ [MeV]',fontdict=font)
  plt.ylabel('W$_y$ [MeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',linestyle=':', linewidth=0.5, color='k')
    #  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
  plt.subplots_adjust(left=0.19, bottom=0.14, right=0.94, top=None,
                    wspace=None, hspace=None)
    
  fig = plt.gcf()
  fig.set_size_inches(9.6, 8)
  fig.savefig('./wrap_work_xy.png',format='png',dpi=160)
  plt.close("all")
  print('./wrap_work_xy.png')
