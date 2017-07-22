import numpy as np
from invdisttree import Invdisttree
import silicon_params as si


def transform_to_uc(wf1, L):
    """
    The function converts a wave functions computed on a grid in a primitive cell to
    the wave functions specified on a user-defined grid in a unit cell.

    :param wf1:  3D real-space wave function computed by ABINIT for a primitive cell
    :param L:    number of points along each of dimenssions
    :return:     3D real-space wave function computed by ABINIT for a unit cell
    """

    a0 = si.a_Si
    num_points = wf1.shape[0]       # number of points along each of dimenssions in wf1 array
    num_cells = 3                   # number of primitive cells needed to form unit cell

    xx = np.linspace(0, num_cells, num_points*num_cells, endpoint=False) - 1.0
    x1, y1, z1 = np.meshgrid(xx,xx,xx)

    wf=np.zeros((num_cells*num_points, num_cells*num_points, num_cells*num_points))

    for j1 in xrange(num_cells):
        for j2 in xrange(num_cells):
            for j3 in xrange(num_cells):

                wf[j1*num_points:((j1+1)*num_points),
                   j2*num_points:((j2+1)*num_points),
                   j3*num_points:((j3+1)*num_points)] = wf1

    x = (y1 + z1) * 0.5 * a0
    y = (x1 + z1) * 0.5 * a0
    z = (x1 + y1) * 0.5 * a0

    F = Invdisttree(np.vstack((x.flatten(), y.flatten(), z.flatten())).T, wf.flatten())

    lin = np.linspace(0, a0, L, endpoint=False)
    x1 ,y1, z1 = np.meshgrid(lin, lin, lin)

    wf = F(np.vstack((x1.flatten(), y1.flatten(), z1.flatten())).T, nnear=11, eps=0, p=1)

    return wf.reshape(x1.shape)