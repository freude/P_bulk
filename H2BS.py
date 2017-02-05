import numpy as np
from coordsys import CoordSys
from spec_func import *

# num_cells=4;
# T=80;

#num_cells=60;
#T=8;

N_states = 10

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

#N_states = 8;

#qn=[1 0  0 -13.6;...
#    2 0  0 -3.4;...
#    3 0  0 -1.511;...
#    4 0  0 -0.85;...
#    5 0  0 -0.544;...
#    6 0  0 -0.3778;...
#    7 0  0 -0.2776;...
#    8 0  0 -0.2125];

num_cells = 40
T = 4

coorsys = CoordSys(num_cells,T,'au')
coorsys.set_origin_cells(num_cells/2+1)
x=coorsys.x()

[X1,Y1,Z1] = np.meshgrid(x,x,x)
s=X1.shape

# Conversion cartesian coordinates to prolate spheroidal ones

Ra = [0, 0, -0.7]
Rb = [0, 0,  0.7]
R = np.sqrt((Rb[0]-Ra[0])**2+(Rb[1]-Ra[1])**2+(Rb[2]-Ra[2])**2);
print R
ra = np.sqrt((X1-Ra[0])**2+(Y1-Ra[1])**2+(Z1-Ra[2])**2)
rb = np.sqrt((X1-Rb[0])**2+(Y1-Rb[1])**2+(Z1-Rb[2])**2)

lam = np.divide((ra+rb), (R))
mu = np.divide((ra-rb), (R))
phi = np.arctan2(Y1, X1)

bas_fun = [x]

for jj in range(N_states):

   if qn[jj][2]==0:
       bas_fun.append(LargeLambda(qn[jj], jj, lam)*LargeM(qn[jj],jj,mu)*np.exp(1j*qn[jj][2]*phi))
   else:
       bas_fun.append(np.imag(LargeLambda(qn[jj], jj, lam)*LargeM(qn[jj],jj,mu)*np.exp(1j*qn[jj][2]*phi)))


   ME = np.trapz(np.trapz(np.trapz(np.abs(bas_fun[jj+1])**2, x, axis=2), x, axis=1), x, axis=0)
   bas_fun[jj+1]=bas_fun[jj+1]/np.sqrt(ME)


# save(strcat(pwd,'/dis_scr/bas_fun.mat'), 'bas_fun', '-v7.3');
#
# a1=qn(1:N_states,4).*40;
# Nbands=N_states;
# save(strcat(pwd,'/dis_scr/M.mat'), 'Nbands', 'a1');
