import numpy as np
import matplotlib.pyplot as plt
from transform_to_uc import transform_to_uc
from read_cube import read_cube



data2d_1 = read_cube('/home/mk/qe_si/results/si_wf.cube')
data2d_1 = transform_to_uc(data2d_1, 30)

data2d_1 = np.concatenate((data2d_1, data2d_1), axis=0)
data2d_1 = np.concatenate((data2d_1, data2d_1), axis=1)
data2d_1 = np.concatenate((data2d_1, data2d_1), axis=2)
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.imshow(data2d_1[15, :, :], cmap='jet')
ax2.imshow(data2d_1[:, 15, :], cmap='jet')
ax3.imshow(data2d_1[:, :, 15], cmap='jet')
plt.show()

from mayavi import mlab
src = mlab.pipeline.scalar_field(data2d_1)
mlab.pipeline.iso_surface(src, contours=[data2d_1.max()-0.3*data2d_1.ptp(), ],)
mlab.pipeline.iso_surface(src, contours=[data2d_1.max()-0.8*data2d_1.ptp(), ],
                          opacity=0.3)

mlab.show()
