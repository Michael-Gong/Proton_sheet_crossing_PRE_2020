import sdf
#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
#from matplotlib import colors, ticker, cm
#from matplotlib.mlab import bivariate_normal
#from optparse import OptionParser
#import os
#from colour import Color

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
#print 'electric field unit: '+str(exunit)
#print 'magnetic field unit: '+str(bxunit)
#print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
    	'weight' : 'normal',  
        'size'   : 20,  
       }  

from_path = './PW_w020_trace2/'
to_path   = './PW_w020_trace2/'
part_name='subset_trace_p/t_electron'
part_mass=1

data = sdf.read(from_path+"trace0450.sdf",dict=True)
grid_x = data['Grid/Particles/'+part_name].data[0]/wavelength
grid_y = data['Grid/Particles/'+part_name].data[1]/wavelength
px = data['Particles/Px/'+part_name].data/(part_mass*m0*v0)
py = data['Particles/Py/'+part_name].data/(part_mass*m0*v0)
theta = np.arctan2(py,px)*180.0/np.pi
gg = (px**2+py**2+1)**0.5
Ek = (gg-1)*part_mass*0.51

part13_id = data['Particles/ID/'+part_name].data
print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))
part13_id = part13_id[ (Ek>80) & (px>0.0)]
#part13_id = part13_id[abs(grid_y)<0.5]
#choice = np.random.choice(range(part13_id.size), 20000, replace=False)
#part13_id = part13_id[choice]
print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))

######### Parameter you should set ###########
start   =  1  # start time
stop    =  700  # end time
step    =  1  # the interval or step

#  if (os.path.isdir('jpg') == False):
#    os.mkdir('jpg')
######### Script code drawing figure ################
part_id = part13_id
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+'trace'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    part_id = np.intersect1d(data['Particles/ID/'+part_name].data, part_id)
    print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))

print('After intersecting with all.sdf')
print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))


######### Parameter you should set ###########
start   =  1  # start time
stop    =  700  # end time
step    =  1  # the interval or step

px_d = np.zeros([part_id.size,stop-start+1])
py_d = np.zeros([part_id.size,stop-start+1])
xx_d = np.zeros([part_id.size,stop-start+1])
yy_d = np.zeros([part_id.size,stop-start+1])
wx_d = np.zeros([part_id.size,stop-start+1])
wy_d = np.zeros([part_id.size,stop-start+1])
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+'trace'+str(n).zfill(4)+".sdf",dict=True)
    px = data['Particles/Px/'+part_name].data/(part_mass*m0*v0)
    py = data['Particles/Py/'+part_name].data/(part_mass*m0*v0)
    grid_x = data['Grid/Particles/'+part_name].data[0]/wavelength
    grid_y = data['Grid/Particles/'+part_name].data[1]/wavelength
    work_x = data['Particles/Time_Integrated_Work_x/'+part_name].data
    work_y = data['Particles/Time_Integrated_Work_y/'+part_name].data
    temp_id =  data['Particles/ID/'+part_name].data

    px = px[np.in1d(temp_id,part_id)]
    py = py[np.in1d(temp_id,part_id)]
    grid_x = grid_x[np.in1d(temp_id,part_id)]
    grid_y = grid_y[np.in1d(temp_id,part_id)]
    work_x = work_x[np.in1d(temp_id,part_id)]
    work_y = work_y[np.in1d(temp_id,part_id)]
    temp_id = temp_id[np.in1d(temp_id,part_id)]

    for ie in range(part_id.size):
        px_d[ie,n-start] = px[temp_id==part_id[ie]]
        py_d[ie,n-start] = py[temp_id==part_id[ie]]
        xx_d[ie,n-start] = grid_x[temp_id==part_id[ie]]
        yy_d[ie,n-start] = grid_y[temp_id==part_id[ie]]
        wx_d[ie,n-start] = work_x[temp_id==part_id[ie]]
        wy_d[ie,n-start] = work_y[temp_id==part_id[ie]]
    print('finised '+part_name+' '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

np.savetxt(to_path+'electron_px.txt',px_d)
np.savetxt(to_path+'electron_py.txt',py_d)
np.savetxt(to_path+'electron_xx.txt',xx_d)
np.savetxt(to_path+'electron_yy.txt',yy_d)
np.savetxt(to_path+'electron_wx.txt',wx_d)
np.savetxt(to_path+'electron_wy.txt',wy_d)
