#!/usr/bin/python

import math
import numpy as np
from confset import ConfSet
from int2e import CoulombInt
import const
import sys

class FCI(CoulombInt):

    def __init__(self,**kw):

        self.Nel=kw.get('Nel',2)
        self.Norb=kw.get('Norb',4)
        self.spin=kw.get('M',0)

        if (((0.5*self.Nel+self.spin) % 1 != 0)|((0.5*self.Nel-self.spin) % 1 != 0)):
            sys.exit("Cann't get desired magnitude of the spin projection")

        self.str_a=ConfSet(Nel = 0.5*self.Nel+self.spin, Norb = self.Norb)
        self.str_b=ConfSet(Nel = 0.5*self.Nel-self.spin, Norb = self.Norb)

        if 'Nel' in kw:
            del kw['Nel']

        if 'Norb' in kw:
            del kw['Norb']

        if 'M' in kw:
            del kw['M']

        CoulombInt.__init__(self,self.Norb,**kw)

        p=kw.get('p','/data/users/mklymenko/science/H2/programing/15d2/')
        sn=kw.get('sn',0)

        try:
            with open(p+"results/EE_"+str(sn)+".dat"):
                self.E=np.loadtxt(p+"results/EE_"+str(sn)+".dat")[0:self.Norb]
        except IOError:
            sys.exit("\n Cann't load energies at"+p+"results/EE_"+str(sn)+".dat")

        self._M=np.zeros((self.str_a.num_conf*self.str_b.num_conf,self.str_a.num_conf*self.str_b.num_conf))

    def ic(self,a,b):
        return self.str_a.num_conf*b+a


    def diagM(self,a1,b1):
        self._M[self.ic(a1.index,b1.index),self.ic(a1.index,b1.index)]=a1.get_energy(self.E)+b1.get_energy(self.E)

    def sigma1(self,a1,b1):

        for i in [i1 for i1, e in enumerate(self.str_b.excitab[b1.index,:]) if e != 0]:

            orb1, orb2 = ConfSet.get_exci_orb(self.str_b.conf[i],b1)

            self._M[self.ic(a1.index,b1.index),self.ic(a1.index,i)] += \
                -0.5*self.sum2(orb1,orb2)*self.str_b.excitab[b1.index,i]

            for j in [j1 for j1, e in enumerate(self.str_b.excitab[self.str_b.conf[i].index,:]) if e != 0]:

                orb3, orb4 = ConfSet.get_exci_orb(self.str_b.conf[j],self.str_b.conf[i])

                self._M[self.ic(a1.index,b1.index),self.ic(a1.index,j)] += \
                    0.5*self.getInt([orb1,orb2,orb3,orb4])*self.str_b.excitab[b1.index,i]*self.str_b.excitab[i,j]


    def sigma2(self,a1,b1):

        for i in [i1 for i1, e in enumerate(self.str_a.excitab[a1.index,:]) if e != 0]:

            orb1, orb2 = ConfSet.get_exci_orb(self.str_a.conf[i],a1)

            self._M[self.ic(a1.index,b1.index),self.ic(i,b1.index)] +=\
                -0.5*self.sum2(orb1,orb2)*self.str_b.excitab[a1.index,i]

            for j in [j1 for j1, e in enumerate(self.str_a.excitab[self.str_a.conf[i].index,:]) if e != 0]:

                orb3, orb4 = ConfSet.get_exci_orb(self.str_a.conf[j],self.str_a.conf[i])

                self._M[self.ic(a1.index,b1.index),self.ic(j,b1.index)] +=\
                    0.5*self.getInt([orb1,orb2,orb3,orb4])*self.str_a.excitab[a1.index,i]*self.str_a.excitab[i,j]

    def sigma3(self,a1,b1):

        for i in [i1 for i1, e in enumerate(self.str_a.excitab[a1.index,:]) if e != 0]:

            orb1,orb2=ConfSet.get_exci_orb(self.str_a.conf[i],a1)

            for j in [j1 for j1, e in enumerate(self.str_b.excitab[b1.index,:]) if e != 0]:

                orb3,orb4=ConfSet.get_exci_orb(self.str_b.conf[j],b1)

                self._M[self.ic(a1.index,b1.index),self.ic(i,j)] +=\
                    self.getInt([orb1,orb2,orb3,orb4])*self.str_a.excitab[a1.index,i]*self.str_a.excitab[b1.index,j]

    def clean_M(self):
        self._M=np.zeros((self.str_a.num_conf*self.str_b.num_conf,self.str_a.num_conf*self.str_b.num_conf))

if __name__ == "__main__":

    paths={'pload' : '/data/users/mklymenko/work_py/mb_project/',\
           'psave' : '/data/users/mklymenko/work_py/mb_project/',\
           'p' : '/data/users/mklymenko/science/H2/programing/15d2/'}

    c=FCI(Nel=2,Norb=3,M=0,sn=0,mf=2,**paths)

    for j1 in c.str_a.conf:
        for j2 in c.str_b.conf:
            c.diagM(j1,j2)
            #c.sigma1(j1,j2)   #nonherm
            #c.sigma2(j1,j2)   #nonherm
            #c.sigma3(j1,j2)

    E,wf=np.linalg.eig(c._M)

    np.savetxt(paths['psave']+"M.dat",c._M)
    np.savetxt(paths['psave']+"EE.dat",E)

    c.str_a.print_conf()
