import numpy as np
from pymech.neksuite import readnek
from pymech import open_dataset
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

### case parameters setting
datapath = '/scratch/shiyud/nekoexamples/Grad_Jump_Penal_test/advecting_cone/'
fieldname = 'data_orig/field0'
ipl = 30 # snapshot index
nelemx = 15 # number of elements in z direction
nelemy = nelemx # number of elements in x direction

### solver parameters setting
GLL_order = 4 # GLL points inside an element

### array initialisation
cone = np.zeros([nelemx*GLL_order,nelemy*GLL_order])

# establish a mapping from 1D field into 3D field
field1d_map = np.array(range(nelemx*nelemy))
field2d_map = np.reshape(field1d_map,[nelemx,nelemy])

ds = open_dataset(fieldname+".f00000")
#%% coordinates
x = ds['x'].data
y = ds['y'].data

filename = datapath+fieldname+'.f'+str(ipl).zfill(5)
dsi = readnek(filename)
element = dsi.elem

for ix_elem in range(nelemx):
    for iy_elem in range(nelemy):
        cone[ix_elem*GLL_order:(ix_elem+1)*GLL_order,\
             iy_elem*GLL_order:(iy_elem+1)*GLL_order] =\
            np.squeeze(element[field2d_map[ix_elem,iy_elem]].temp[0,0,:,:])

### plotting
X, Y = np.meshgrid(x,y) 
fig = plt.figure(figsize = (8,8))
ax = plt.axes(projection='3d')
ax.plot_surface(X,Y,cone,cmap = plt.cm.cividis)
ax.set_xlim([-2,2])
ax.set_ylim([-2,2])
ax.set_zlim([-0.5,1])
plt.show()

