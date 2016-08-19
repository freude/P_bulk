
"""The module contains the XStr class."""

import numpy as np


class XStr(object):

    """Single slater determinant class."""

    def __init__(self, **kw):

        self._id = kw.get('id', np.array([1, 0, 0, 0]))  # id of the string
        self.K = len(self._id)                           # number of orbitals
        self.N = np.sum(self._id)                        # number of electrons in the string
        self._index = self.get_index()                   # index the string

    @property
    def id(self):
        """I'm the 'x' property."""
        return self._id

    @property
    def index(self):
        """I'm the 'x' property."""
        return self._index

    @id.setter
    def id(self, value):
        self._id = value

    def get_index(self):
        """This function assigns the unique index to the configuration."""
        auxiliary_mat = np.ones((self.K - self.N + 1, self.N + 1))

        for orb_ind in xrange(1, (self.K - self.N + 1)):
            for el_ind in xrange(1, (self.N + 1)):
                auxiliary_mat[orb_ind, el_ind] = auxiliary_mat[orb_ind - 1, el_ind] + auxiliary_mat[orb_ind, el_ind - 1]

        YY = np.zeros((self.K - self.N + 1, self.N))
        YY[1:, :] = auxiliary_mat[:-1, 1:]

        j1=0; j2=0; ind=0

        for orb_ind in xrange(self.K):
            if (self._id[orb_ind] == 1):
                ind = ind + YY[j2, j1]
                j1 += 1
            else:
                j2 += 1

        return int(ind)

    def get_energy(self, E):
        return sum([self._id[j] * E[j] for j in xrange(self.K)])
