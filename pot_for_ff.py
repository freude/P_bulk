import os
import cPickle as pickle
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from me2 import me2
import silicon_params as si


def read_txt_between_delimiters(file, del1, del2):

    block = ""

    with open(file) as input_data:
        # Skips text before the beginning of the interesting block:
        for line in input_data:
            if line.strip() == del1:  # Or whatever test is needed
                break
        # Reads text until the end of the block:
        for line in input_data:  # This keeps reading the file
            if line.strip() == del2:
                break
            block+=line

    return block


def parse_mesh_file():

    pwd = os.getcwd()
    pwd = os.path.join(pwd, 'dis_scr/mesh_sample.mesh')
    ind = read_txt_between_delimiters(pwd, 'Vertices', 'Tetrahedra')
    data_lines = ind.splitlines()
    data_lines = data_lines[1:-1]
    b = np.array([np.fromstring(j, dtype=float, sep=' ') for j in data_lines])
    return b[:, :3]


def pot_for_ff(k1, k2, file_ind):

    # -----------------------apply filter function - --------------------

    pwd = os.getcwd()

    print('Computing the smoothed potential on a regular grid...'),
    x, V1sm = me2(k1, k2, 'pot')
    print('Done!\n')

    # ----------------------------------------------------------------

    print('Read the mesh...'),
    b = parse_mesh_file()
    print('Done!\n')

    print('Build the interpolant for the potential...'),

    # prepare data

    if (k1 == k2).all():
        x = x[::2]
        V1sm = np.real(V1sm[::2, ::2, ::2])

    if not np.iscomplexobj(V1sm):
        if os.path.isfile(os.path.join(pwd, 'dis_scr/F' + str(file_ind) + '.pkl')):
            with open(os.path.join(pwd, 'dis_scr/F' + str(file_ind) + '.pkl'), 'rb') as input:
                F = pickle.load(input)
        else:
            V1sm = np.transpose(V1sm, (1, 0, 2))
            F = RegularGridInterpolator((x, x, x), V1sm)
            with open(os.path.join(pwd, 'dis_scr/F'+ str(file_ind) + '.pkl'), 'wb') as output:
                pickle.dump(F, output)

        print('Done!\n')

        print('Compute the potetnial on the mesh...'),
        pot = F(b)
        pot = np.nan_to_num(pot)
        print('Done!\n')

        print('Save the potential and the mesh...'),
        np.save(os.path.join(pwd, 'dis_scr/pot' + str(file_ind) + '.npy'), pot)
        np.save(os.path.join(pwd, 'dis_scr/mesh.dat'), b)
        print('Done!\n')

    else:

        if os.path.isfile(os.path.join(pwd, 'dis_scr/F'+ str(file_ind) + 'r.pkl')) and\
                os.path.isfile(os.path.join(pwd, 'dis_scr/F'+ str(file_ind) + 'i.pkl')):

            with open(os.path.join(pwd, 'dis_scr/F'+ str(file_ind) + 'r.pkl'), 'rb') as input:
                Fr = pickle.load(input)
            with open(os.path.join(pwd, 'dis_scr/F'+ str(file_ind) + 'i.pkl'), 'rb') as input:
                Fi = pickle.load(input)

        else:
            V1sm_r = np.transpose(np.real(V1sm), (1, 0, 2))
            V1sm_i = np.transpose(np.imag(V1sm), (1, 0, 2))
            Fr = RegularGridInterpolator((x, x, x), V1sm_r)
            Fi = RegularGridInterpolator((x, x, x), V1sm_i)

            with open(os.path.join(pwd, 'dis_scr/F'+ str(file_ind) + 'r.pkl'), 'wb') as output:
                pickle.dump(Fr, output)
            with open(os.path.join(pwd, 'dis_scr/F'+ str(file_ind) + 'i.pkl'), 'wb') as output:
                pickle.dump(Fi, output)

        print('Done!\n')

        print('Compute the potetnial on the mesh...'),
        pot_r = Fr(b[:, 0], b[:, 1], b[:, 2])
        pot_i = Fi(b[:, 0], b[:, 1], b[:, 2])
        pot_r = np.nan_to_num(pot_r)
        pot_i = np.nan_to_num(pot_i)
        print('Done!\n')

        print('Save the potential and the mesh...'),
        np.save(os.path.join(pwd, 'dis_scr/pot' + str(file_ind) + '.npy'), pot_r)
        np.save(os.path.join(pwd, 'dis_scr/pot' + str(file_ind) + '.npy'), pot_i)
        np.save(os.path.join(pwd, 'dis_scr/mesh.dat'), b)
        print('Done!\n')


if __name__ =='__main__':

    k0 = si.k0 * si.ab

    kk = k0 * np.array([[1, 0, 0],
                        [-1, 0, 0],
                        [0, 1, 0],
                        [0, -1, 0],
                        [0, 0, 1],
                        [0, 0, -1]])

    pot_for_ff(kk[0, :], kk[0, :], '1')
    pot_for_ff(kk[1, :], kk[1, :], '2')
    pot_for_ff(kk[2, :], kk[2, :], '3')