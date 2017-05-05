import matplotlib.pyplot as plt
import numpy as np
from coordsys import CoordSys
from spec_func import *

# num_cells=4;
# T=80;

# num_cells = 50
# T = 6

num_cells = 50
T = 1


N_states = 2

#      n  l  m      E        p      sigma
qn = [[1, 0, 0, -2.56853, 1.12186, 0.24792],      #1s_sigma_g
      [2, 0, 0, -0.78975, 0.62208, 1.25053],      #2s_sigma_g
      [3, 0, 0, -0.37851, 0.430664, 2.25079],     #3s_sigma_g
      [2, 1, 0, -1.22415, 0.77449, 0.80764],      #2p_sigma_u
      [3, 1, 0, -0.49729, 0.493633, 1.83612],     #3p_sigma_u
      [4, 1, 0, -0.27065, 0.364171, 2.84435],     #4p_sigma_u
      [3, 2, 0, -0.45660, 0.473004, 1.95981],     #3d_sigma_g
      [4, 3, 0, -0.25158, 0.351105, 2.98741],     #4f_sigma_u
      [2, 1, 1, -0.91265, 0.66873, 0.09352],      #2p_pi_u
      [3, 2, 1, -0.44939, 0.469258, 0.98343]]     #3d_pi_g

coorsys = CoordSys(num_cells,T,'au')
coorsys.set_origin_cells(num_cells/2+1)
x=coorsys.x()
[X1,Y1,Z1] = np.meshgrid(x,x,x)
s=X1.shape

# Conversion cartesian coordinates to prolate spheroidal ones

Ra = [0, 0, -0.7]
Rb = [0, 0,  0.7]
R = np.sqrt((Rb[0]-Ra[0])**2+(Rb[1]-Ra[1])**2+(Rb[2]-Ra[2])**2)
ra = np.sqrt((X1-Ra[0])**2+(Y1-Ra[1])**2+(Z1-Ra[2])**2)
rb = np.sqrt((X1-Rb[0])**2+(Y1-Rb[1])**2+(Z1-Rb[2])**2)

lam = np.divide((ra+rb), (R))
mu = np.divide((ra-rb), (R))
phi = np.arctan2(Y1, X1)

bas_fun = [x]
# X1, Y1, Z1 = X1.flatten(), Y1.flatten(), Z1.flatten()

for jj in range(N_states):
    if qn[jj][2]==0:
       bas_fun.append(LargeLambda(qn[jj], jj, lam)*LargeM(qn[jj],jj,mu)*np.cos(qn[jj][2]*phi))
    else:
       bas_fun.append(LargeLambda(qn[jj], jj, lam)*LargeM(qn[jj],jj,mu)*np.sin(qn[jj][2]*phi))

    bas_fun[jj + 1] = np.nan_to_num(bas_fun[jj + 1])

    ME = np.trapz(np.trapz(np.trapz(np.abs(bas_fun[jj+1])**2, x, axis=2), x, axis=1), x, axis=0)
    print(ME)
    bas_fun[jj+1]=bas_fun[jj+1]/np.sqrt(ME)


# a=plt.contour((bas_fun[1][:,125,:]), 30); plt.colorbar()
# plt.imshow((bas_fun[1][:,125,:]))
# plt.show()

from fci.gauss_fit import GFit

qn=1

# data = np.vstack((X1, Y1, Z1, bas_fun[qn].flatten())).T

wf = GFit(init='fit',
          sn=10,
          qn=qn,
          mf=2,
          num_fu=12,
          psave='./',
          pload='./')

wf.do_fit(bas_fun[1], x=X1, y=Y1, z=Z1)

x = np.linspace(-6.5, 6.5, 300)
y = np.linspace(4.0, 8.9, 300)

XXX = np.vstack((x, x * 0.0 + 6.45, 0.0 * x))

xi, yi = np.meshgrid(x, y)
x, y = xi.flatten(), yi.flatten()
z = x * 0.0
XX = np.vstack((x, y, z))

# wf.save()
print(wf._gf)
# wf.draw_func(x,y,par='2d')
g = wf.show_func(XX)
g1 = wf.show_gf(XXX)
AA = wf.get_value(XX.T)

fig = plt.figure()
# for j in range(0,wf._num_fu):
#     plt.plot(XXX[0,:].T,g1[:,j])

# plt.plot(xi[150, :].T, g.reshape(xi.shape)[150, :])
# plt.plot(xi[150, :].T, AA.reshape(xi.shape)[150, :])

# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(xi,yi,g1.reshape(xi.shape), cmap=cm.jet, linewidth=0.2)

plt.contour(xi, yi, -AA.reshape(xi.shape), colors='red')
plt.contour(xi, yi, -g.reshape(xi.shape), colors='blue')

plt.hold(True)
plt.show()



