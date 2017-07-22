import numpy as np
from abi_read import abi_read
from per_freq import per_freq
from smoother_2 import smoother_2, find
from pot_mat import pot_mat
from scipy.integrate import simps
from coordsys import CoordSys


def me2(k1, k2, flag):
    """
    If 'me', the function computes the integral
    over a coupling filtered potential

    If 'pot', the function computes the filtered potential itself

    :param k1:
    :param k2:
    :param flag:
    :return:
    """

    if np.array_equal(k1, k2) and flag == 'mes':
    
        ME = 0    # matrix element is usually needed only for nondiagonal term
    
    else:
        # ----- compute Fourier harmonics of the periodic Bloch fucntion -----
    
        # if (exist(strcat(pwd, '/dis_scr/G.mat'), 'file') == 2)
            # G = load(strcat(pwd, '/dis_scr/G.mat'));
        # G = G.G;
        # else
        num_cells = 1
        T = 70
        wf = np.conj(abi_read(num_cells, T, k1)) * (abi_read(num_cells, T, k2))
        G = per_freq(wf)
        # save(strcat(pwd, 'dis_scr/G.mat'), 'G');
        # end;

        # -----------------------apply filter function ----------------------

        if np.array_equal(k1, k2):
            num_cells = 120
            coorsys = CoordSys(num_cells, 2, 'au')
            coorsys.set_origin_cells(num_cells / 2 + 1)
        else:
            num_cells = 30
            coorsys = CoordSys(num_cells, 16, 'au')
            coorsys.set_origin_cells(num_cells / 2 + 1)

        V, _, _, _ = pot_mat(coorsys, k1, k2)
        x = coorsys.x()

        # -----------------------apply filterfunction - --------------------

        # fprintf('apply filter function...')
        R_tr = 7
        V1sm = smoother_2(x, G, k1, k2, R_tr)
        # fprintf('Done!\n')

        if flag == 'pot':

            jj1 = np.argmin(np.abs(x + 0.65 * R_tr))
            jj2 = np.argmin(np.abs(x - 0.65 * R_tr))

            # delta = V1sm(jj2, jj2, jj2) - V(jj2, jj2, jj2);
            # V = V + delta;

            V[jj1:jj2, jj1:jj2, jj1:jj2] = -V1sm[jj1:jj2, jj1:jj2, jj1:jj2]
            V1sm = -V

    # ----------------------------------------------------------------
    
    M2 = 1
    
    if (((k1[find(k1)] > 0) and (k2[find(k2)] < 0)) or
            ((k2[find(k2)] > 0) and (k1[find(k1)] < 0))) and\
            (find(k1) == find(k2)):
        M1 = 1
    else:
        M1 = 1

    
    [X, Y, Z] = np.meshgrid(x, x, x);
    kk = 0 * k2 - 0 * k1
    
    if flag == 'mes':
        ME = 2.7 * simps(simps(simps(V1sm * M1 * M2 * np.exp(1j * (kk[0] * X + kk[1] * Y + kk[2] * Z)), x), x), x)
    elif flag == 'pot':
        ME = (x, V1sm)
    else:
        pass

    return ME

