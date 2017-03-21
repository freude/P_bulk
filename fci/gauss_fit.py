"""
    The class is for fitting a grid function by a set of gaussian functions
"""

from IPython.core.debugger import Tracer
import glob
import os
import sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import Rbf
import math
# -----------------------------------------------------------
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
# -----------------------------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from fci.invdisttree import Invdisttree
# -----------------------------------------------------------
import numpy as np
from scipy.ndimage.filters import maximum_filter, minimum_filter
from scipy.ndimage import binary_dilation
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import matplotlib.pyplot as pp
from fci.aux_fun import mat2table

class GFit(object):

    def __init__(self, data = None, **kw):

        # path to load/save fitting parameters
        self._save = kw.get('psave', '/data/users/mklymenko/work_py/mb_project/')

        self._num_fu = kw.get('num_fu', 55)  # number of Gaussian primitive functions
        self._sn = kw.get('sn', 0)           # id of the system
        self._qn = kw.get('qn', 0)           # id of the wave function
        self._flag = kw.get('mf', 2)         # choice of the model function

        self._save = os.path.join(self._save, str(self._sn), '_', str(self._qn), '.dat')

        self._gf = []
        self.peaks = []
        self.coords = []

        print('\n---------------------------------------------------------------------')
        print('The fitting is provided for the function')
        print(str(self._sn) + '_' + str(self._qn))
        print('The number of primitive Gaussians is {}'.format(self._num_fu))
        print('---------------------------------------------------------------------\n')

        if data is not None:
            self.do_fit(data)

    @staticmethod
    def resample(data, xx, yy, zz):

        # xx, yy, zz = np.mgrid[-6.5:6.5:160j, 4.0:8.9:160j, -7.5:7.5:160j]
        xx, yy, zz = xx.flatten(), yy.flatten(), zz.flatten()
        xx = np.vstack((xx, yy, zz))

        invdisttree = Invdisttree(data[:, :3], data[:, 3:], leafsize=10, stat=1)
        return xx, invdisttree(xx.T, nnear=130, eps=0, p=1)

    @staticmethod
    def detect_min_max(arr):
        """That's a very nice function detecting all local minima and maxima
        and computing their coordinates.
        The method is based on derivatives.
        """

        x1, x2, y1, y2, z1, z2 = GFit.coords_of_cube(arr)
        steps = 160
        xi, yi, zi = np.mgrid[x1:x2:100j, y1:y2:100j, z1:z2:100j]
        _, arr = GFit.resample(arr, xi, yi, zi)
        arr = arr.reshape(xi.shape)

        return detect_peaks(xi, yi, zi, arr)

        # max_value = max(np.absolute(np.reshape(arr, -1)))
        # peaks_max = []
        # peaks_min = []
        # x_max = []
        # y_max = []
        # z_max = []
        # x_min = []
        # y_min = []
        # z_min = []
        #
        # for j1 in range(10, arr.shape[0]-10):
        #     for j2 in range(10, arr.shape[1]-10):
        #         for j3 in range(10, arr.shape[2]-10):
        #             if (np.absolute(arr[j1, j2, j3]) > 0.3*max_value):
        #
        #                 aaaa = [
        #                     arr[j1, j2, j3 + 1], arr[j1, j2 + 1, j3],
        #                     arr[j1 + 1, j2, j3], arr[j1, j2, j3 - 1],
        #                     arr[j1, j2 - 1, j3], arr[j1 - 1, j2, j3],
        #                     arr[j1 + 1, j2 + 1, j3 + 1],
        #                     arr[j1 - 1, j2 - 1, j3 - 1],
        #                     arr[j1 - 1, j2 + 1, j3 + 1], arr[j1, j2 + 1, j3 + 1],
        #                     arr[j1, j2 - 1, j3 - 1], arr[j1, j2 - 1, j3 + 1],
        #                     arr[j1, j2 + 1, j3 - 1], arr[j1 + 1, j2, j3 + 1],
        #                     arr[j1 - 1, j2, j3 - 1], arr[j1 - 1, j2, j3 + 1],
        #                     arr[j1 + 1, j2, j3 - 1], arr[j1 + 1, j2 + 1, j3],
        #                     arr[j1 - 1, j2 - 1, j3], arr[j1 + 1, j2 - 1, j3],
        #                     arr[j1 - 1, j2 + 1, j3], arr
        #                     [j1 + 1, j2 - 1, j3 + 1], arr
        #                     [j1 + 1, j2 + 1, j3 - 1], arr
        #                     [j1 - 1, j2 - 1, j3 + 1], arr
        #                     [j1 + 1, j2 - 1, j3 - 1], arr
        #                     [j1 - 1, j2 + 1, j3 - 1]]
        #                 bbbb = [
        #                     arr[j1, j2, j3 + 9], arr[j1, j2 + 9, j3],
        #                     arr[j1 + 9, j2, j3], arr[j1, j2, j3 - 9],
        #                     arr[j1, j2 - 9, j3], arr[j1 - 9, j2, j3]]
        #
        #                 if ((arr[j1, j2, j3] > max(aaaa)) and (max(aaaa) > max(bbbb))):
        #                     peaks_max = np.append(peaks_max, arr[j1, j2, j3])
        #                     x_max = np.append(x_max, xi[j1, j2, j3])
        #                     y_max = np.append(y_max, yi[j1, j2, j3])
        #                     z_max = np.append(z_max, zi[j1, j2, j3])
        #
        #                 if ((arr[j1, j2, j3] < min(aaaa)) and (min(aaaa) < min(bbbb))):
        #                     peaks_min = np.append(peaks_min, arr[j1, j2, j3])
        #                     x_min = np.append(x_min, xi[j1, j2, j3])
        #                     y_min = np.append(y_min, yi[j1, j2, j3])
        #                     z_min = np.append(z_min, zi[j1, j2, j3])
        #
        # return peaks_min, np.vstack(
        #     (x_min, y_min, z_min)), peaks_max, np.vstack(
        #     (x_max, y_max, z_max))

    @staticmethod
    def coords_of_cube(data):

        x_min = np.min(data[:, 0])
        x_max = np.max(data[:, 0])
        y_min = np.min(data[:, 1])
        y_max = np.max(data[:, 1])
        z_min = np.min(data[:, 2])
        z_max = np.max(data[:, 2])

        return x_min, x_max, y_min, y_max, z_min, z_max

    def save(self):
        """Save Gaussian functions coefficients to the file"""

        if self._save != '0':
            p = self._save+self._path[-3:-1]+'_'+str(self._qn)+'.dat'
            np.savetxt(p, self._gf)
        else:
            sys.exit("Wrong path to save")

# -------------------------------------------------------------------------
# -------------------------------model functions---------------------------
# -------------------------------------------------------------------------

#    def modelfun(self, x, *par):
#
#        g = np.zeros(len(x[0]))
#
#        for j in range(0, len(self._sites)):
#            r1 = pow((x[0]-self._sites[j, 0]), 2) + \
#                pow((x[1]-self._sites[j, 1]), 2)+pow((x[2]-self._sites[j, 2]), 2)
#            for j1 in range(0, (len(par)/2/len(self._sites))):
#                g = g+par[j*len(par)/len(self._sites)/2+j1] * \
#                    np.exp(-r1/par[len(par)/2+j*len(par)/len(self._sites)/2+j1])
#
#        return g

    def modelfun1(self, x, *par):
        """
        The model function represented by a sum of
        the Gaussian functions with variable positions, widths and
        amplitudes
        """

        g = np.zeros(len(x[0]))

        for j in range(len(par)//5):
            x1 = par[j*5]
            x2 = par[j*5+1]
            x3 = par[j*5+2]
            w = par[j*5+3]
            a = par[j*5+4]
            r1 = pow((x[0]-x1), 2)+pow((x[1]-x2), 2)+pow((x[2]-x3), 2)
            # if ((a > 1.1) or (a < -1.1)): a=0
            g = g+a*np.exp(-r1/abs(w))

        return g

# -------------------------------------------------------------------------
# --------------Tool for extracting values of the wave function------------
# -------------------------------------------------------------------------

    def show_func(self, x):
        """
        Computes the value of the wave function in points stored in
        the vector x using fitting parameters and the model functions.
        """

        if (self._flag == 1):
            g = self.modelfun(x, *self._gf)
        elif (self._flag == 2):
            g = self.modelfun1(x, *self._gf)
        elif ((self._flag == 0) & (self._load != '0')):
            pass
        else:
            # pass
            sys.exit("Wrong flag in do_fit")

        return g

    def show_gf(self, x):
        """Same as show_func(self,x) but returns
        decomposed primitive Gaussian functions
        """
        g = np.zeros((len(x[0]), self._num_fu), dtype=np.float64)
        for j in range(self._num_fu):
            x1 = self._gf[j*5]
            x2 = self._gf[j*5+1]
            x3 = self._gf[j*5+2]
            w = self._gf[j*5+3]
            a = self._gf[j*5+4]
            r1 = pow((x[0]-x1), 2)+pow((x[1]-x2), 2)+pow((x[2]-x3), 2)
            g[:, j] = a*np.exp(-r1/abs(w))

        return g

# ----------------------------------------------------------------------

    @staticmethod
    def moments(data):
        """Returns (height, x, y, width_x, width_y)
         the gaussian parameters of a 2D distribution by calculating its
        moments """

        data = np.absolute(data)
        total = data.sum()
        X = np.indices(data.shape)
        x = (X*data).sum()/total
        width = np.sqrt((((X-x)**2)*data).sum()/data.sum())
        m_max = data.max()
        m_min = data.min()
        if np.absolute(m_max) >= np.absolute(m_min):
            height = m_max
        else:
            height = m_min
        return height, x, width

#-----------------------------------------------------------------------

    def do_fit(self, data, x=None, y=None, z=None):
        """ The function does the fitting procedure"""

        if (self._flag == 1):

            # The initial values of parameters are stored in self._gf

            popt, pcov = curve_fit(self.modelfun, data[:, :3], data[:, :3], p0=self._gf, maxfev=5000)
            self._gf = popt

        elif (self._flag == 2):

            if not self.peaks:
                print('\n---------------------------------------------------------------------')
                print('Detecting maxima and minima of target function...',)

                if x is None:
                    peaks_min, min_coord, peaks_max, max_coord = self.detect_min_max(data)
                else:
                    peaks_min, min_coord, peaks_max, max_coord = detect_peaks(x, y, z, data)

                print('done')
                print('Number of the min peaks: {}'.format(len(peaks_min)))
                print('Number of the max peaks: {}'.format(len(peaks_max)))
                print('---------------------------------------------------------------------\n')

                if peaks_max == []:
                    self.peaks=np.insert(peaks_min, np.arange(len(peaks_max)), peaks_max)
                    self.coords=np.insert(min_coord, np.arange(max_coord.shape[1]), max_coord, axis=1)
                else:
                    self.peaks = np.insert(peaks_max, np.arange(len(peaks_min)), peaks_min)
                    self.coords = np.insert(max_coord, np.arange(min_coord.shape[1]), min_coord, axis=1)

            if x is not None:
                data = mat2table(x, y, z, data)

            par = [0.0]*(self._num_fu*5)
            j1 = 0
            aaaa = 1
            for j in range(self._num_fu):
                if (j > aaaa*self.coords.shape[1]-1):
                    j1 = 0
                    aaaa += 1
                par[j*5] = self.coords[0, j1]
                par[j*5+1] = self.coords[1, j1]
                par[j*5+2] = self.coords[2, j1]
                # par[j*5+3] = 0.1003+0.1000*math.copysign(1, (pow(-1, j)))
                par[j*5+3] = 0.0001
#                if j < 15:
#                    par[j*5+3] = 0.00001
#                else:
#                    par[j*5+3] = 0.0005
                par[j*5+4] = self.peaks[j1]
#                print(coords[0, j1], coords[1, j1], coords[2, j1])
                j1 += 1
            popt, pcov = curve_fit(
                self.modelfun1, data[:, :3].T, np.squeeze(data[:, 3:]), p0=par)
            self._gf = popt
#             self.error=np.diagonal(pcov, offset=0)
#             print(pcov)
        else:
            sys.exit("Wrong flag in do_fit")


def detect_peaks(x, y, z, image):
    """
    Takes an image and detect the peaks usingthe local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define an 8-connected neighborhood
    neighborhood = generate_binary_structure(3, 3)
    a = np.zeros((5, 5, 5))
    a[2, 2, 2] = 1
    a = binary_dilation(a, structure=neighborhood).astype(a.dtype)
    neighborhood = binary_dilation(a, structure=neighborhood).astype(bool)

    #apply the local maximum filter; all pixel of maximal value
    #in their neighborhood are set to 1
    local_max = maximum_filter(image, footprint=neighborhood)==image

    # #local_max is a mask that contains the peaks we are
    # #looking for, but also the background.
    # #In order to isolate the peaks we must remove the background from the mask.
    #
    # #we create the mask of the background
    # background = (image==0)
    #
    # #a little technicality: we must erode the background in order to
    # #successfully subtract it form local_max, otherwise a line will
    # #appear along the background border (artifact of the local maximum filter)
    # eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
    #
    # #we obtain the final mask, containing only peaks,
    # #by removing the background from the local_max mask (xor operation)
    # detected_peaks = local_max ^ eroded_background

    x_peaks_max = x[local_max]
    y_peaks_max = y[local_max]
    z_peaks_max = z[local_max]
    peaks_max = image[local_max]

    if peaks_max.any():
        norm = np.max(peaks_max)
        x_peaks_max = x_peaks_max[peaks_max > 0.3 * norm]
        y_peaks_max = y_peaks_max[peaks_max > 0.3 * norm]
        z_peaks_max = z_peaks_max[peaks_max > 0.3 * norm]
        peaks_max = peaks_max[peaks_max > 0.3 * norm]

    local_max = minimum_filter(image, footprint=neighborhood) == image

    x_peaks_min = x[local_max]
    y_peaks_min = y[local_max]
    z_peaks_min = z[local_max]
    peaks_min = image[local_max]

    if peaks_min.any():
        norm = np.min(peaks_min)
        x_peaks_min = x_peaks_min[peaks_min < 0.3 * norm]
        y_peaks_min = y_peaks_min[peaks_min < 0.3 * norm]
        z_peaks_min = z_peaks_min[peaks_min < 0.3 * norm]
        peaks_min = peaks_min[peaks_min < 0.3 * norm]

    return peaks_min, np.vstack((x_peaks_min, y_peaks_min, z_peaks_min)),\
           peaks_max, np.vstack((x_peaks_max, y_peaks_max, z_peaks_max))
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

if __name__ == "__main__":

    x = np.linspace(-6.5, 6.5, 300)
    y = np.linspace(4.0, 8.9, 300)

    XXX = np.vstack((x, x*0.0+6.45, 0.0*x))

    xi, yi = np.meshgrid(x, y)
    x, y = xi.flatten(), yi.flatten()
    z = x*0.0
    XX = np.vstack((x, y, z))

    wf = GFit(
        init='fit',
        sn=10,
        qn=1,
        mf=2,
        num_fu=14,
        psave='/data/users/mklymenko/work_py/mb_project/',
        pload='/data/users/mklymenko/work_py/mb_project/')

    wf.save()
    print(wf._gf)
    # wf.draw_func(x,y,par='2d')
    g = wf.show_func(XX)
    g1 = wf.show_gf(XXX)
    AA = wf.get_value(XX.T)

    fig = plt.figure()
    # for j in range(0,wf._num_fu):
    #     plt.plot(XXX[0,:].T,g1[:,j])

    #plt.plot(xi[150, :].T, g.reshape(xi.shape)[150, :])
    #plt.plot(xi[150, :].T, AA.reshape(xi.shape)[150, :])

    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(xi,yi,g1.reshape(xi.shape), cmap=cm.jet, linewidth=0.2)

    plt.contour(xi, yi, -AA.reshape(xi.shape), colors='red')
    plt.contour(xi, yi, -g.reshape(xi.shape), colors='blue')

    plt.hold(True)
    plt.show()
