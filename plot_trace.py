import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os

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
        'style'  : 'normal',
        'color'  : 'black',  
    'weight' : 'normal',  
        'size'   : 14,  
       }  


directory = './PW_w020_trace2/'
to_path   = './PW_w020_trace2_fig/'
px = np.loadtxt(directory+'electron_px.txt')
py = np.loadtxt(directory+'electron_py.txt')
xx = np.loadtxt(directory+'electron_xx.txt')
yy = np.loadtxt(directory+'electron_yy.txt')
wx = np.loadtxt(directory+'electron_wx.txt')
wy = np.loadtxt(directory+'electron_wy.txt')

#ey_averaged = -8.0/3.2*yy
#bz_averaged = -8.0/3.2*yy

#laser_ey = fieldey-ey_averaged
#laser_bz = fieldbz-bz_averaged


gg = (px**2+py**2+1)**0.5
R = gg-px
theta = np.arctan2(py,px)

#number=400

scatter_size = 40
tt = np.linspace(5.0,89.9,850)
#tt = tt[:,np.newaxis]

for index in range(gg[:,0].size):
    plt.subplot(2,2,1)
    norm = matplotlib.colors.Normalize(vmin=1.,vmax=150.)
    plt.plot(xx[index,:], yy[index,:],':k',linewidth=0.5)
    plt.scatter(xx[index,:], yy[index,:], c=gg[index,:], s=scatter_size, cmap='rainbow', edgecolors='None')
#    cbar=plt.colorbar( ticks=np.linspace(np.min(1), np.max(601), 5) )
    cbar=plt.colorbar()
    cbar.set_label(r'$\gamma$',fontdict=font)
    #   plt.legend(loc='upper right')
    #plt.xlim(-500,900)
#    plt.ylim(-4.9,4.9)
    plt.xlabel(r'$x\ [\mu m]$',fontdict=font)
    plt.ylabel(r'$y\ [\mu m]$',fontdict=font)
    #plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)


    plt.subplot(2,2,2)
    norm = matplotlib.colors.Normalize(vmin=1.,vmax=150.)
    plt.plot(px[index,:], py[index,:],':k',linewidth=0.5)
    plt.scatter(px[index,:], py[index,:], c=gg[index,:], s=scatter_size, cmap='rainbow', edgecolors='None')
#    cbar=plt.colorbar( ticks=np.linspace(np.min(1), np.max(601), 5) )
    cbar=plt.colorbar()
    cbar.set_label(r'$\gamma$',fontdict=font)
    #   plt.legend(loc='upper right')
    #plt.xlim(-500,900)
#    plt.ylim(-4.9,4.9)
    plt.xlabel(r'$p_x\ [m_ec]$',fontdict=font)
    plt.ylabel(r'$p_y\ [m_ec]$',fontdict=font)
    #plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)



    plt.subplot(2,2,3)
#    norm = matplotlib.colors.Normalize(vmin=1.,vmax=50)
    plt.plot(xx[index,:],px[index,:],':k',linewidth=0.5)
    plt.scatter(xx[index,:],px[index,:], c=R[index,:], s=scatter_size, cmap='rainbow', edgecolors='None')
    cbar=plt.colorbar()
    cbar.set_label(r'$R$',fontdict=font)
    #   plt.legend(loc='upper right')
    #plt.xlim(-500,900)
    #plt.ylim(-4.9,4.9)
    plt.xlabel(r'$x\ [\mu m]$',fontdict=font)
    plt.ylabel(r'$p_x\ [m_ec]$',fontdict=font)
    #plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
 
    

    plt.subplot(2,2,4)
    norm = matplotlib.colors.Normalize(vmin=1.,vmax=150)
    plt.plot(wx[index,:],wy[index,:],':k',linewidth=0.5)
    plt.scatter(wx[index,:],wy[index,:], c=gg[index,:],  s=scatter_size, cmap='rainbow', edgecolors='None')
    cbar=plt.colorbar()
    cbar.set_label(r'$\gamma$',fontdict=font)
    #   plt.legend(loc='upper right')
    #plt.xlim(-500,900)
    #plt.ylim(-4.9,4.9)
    plt.xlabel(r'$Wx\ [m_ec]$',fontdict=font)
    plt.ylabel(r'$Wy\ [m_ec]$',fontdict=font)
    #plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
    
    #plt.show()
    #lt.figure(figsize=(100,100))
    fig = plt.gcf()
    fig.set_size_inches(24, 15)
    fig.savefig(to_path+'e_trace'+str(index).zfill(4)+'.png',format='png',dpi=80)
    plt.close("all")
    print(to_path+'e_trace'+str(index).zfill(4)+'.png')
