#!/usr/bin/python

from numpy import *
from collections import defaultdict

import UNIFACGroup
import Molecule
from Phase import Phase

R = 8.314

################################################################################
################################################################################
class UNIFAC(Phase):

    def __init__(self, name, group_list, A, B):
        self.name = name
        self.group_list = group_list
        self.A = A
        self.B = B
        self.molecule_dict = {}
        self.conc = {}
        self.J = {}
        self.L = {}
        self.x = {}
        self.T = 298.15
        self.p = 1.0e5

    def methodName(self):
        print 'UNIFAC'

    def get_tau(self, g1, g2):
        a = self.A[g1][g2] + self.B[g1][g2]/self.T
        return exp(-a/self.T)

    def print_group_list(self):
        for name in self.group_list.keys():
            print name

    def get_s(self, k, i):
        mol = self.molecule_dict[i]
        G = mol.get_G()
        s = 0.0
        for m in G:
            s += G[m]*self.get_tau(m,k)
        return s

    def print_s(self):
        for k in self.group_list:
            for i in self.molecule_dict:
                print k, i, self.get_s(k, i)

    def set_theta(self):
        self.set_x()
        self.theta = {}
        for k in self.group_list.keys():
            self.theta[k] = 0.0
            for i in self.molecule_dict:
                mol = self.molecule_dict[i]
                G = mol.get_G()
                if (k in G):
                    self.theta[k] += G[k]*self.x[i]

    def set_eta(self):
        self.set_x()
        self.eta = {}
        for k in self.group_list:
            self.eta[k] = 0.0
            for i in self.molecule_dict:
                self.eta[k] += self.get_s(k,i)*self.x[i]

    def print_eta(self):
        for k in self.group_list:
            print k, self.eta[k]

    def print_theta(self):
        for k in self.group_list:
            print k, self.theta[k]

    def get_lngamma(self):
        self.set_x()
        self.set_eta()
        self.set_theta()
        sum_J = 0.0
        sum_L = 0.0
        J = {}
        L = {}
        term = {}
        for i in self.molecule_dict:
            mol = self.molecule_dict[i]
            r = mol.get_r()
            q = mol.get_q()
            sum_J += r*self.x[i]
            sum_L += q*self.x[i]
            J[i] = r
            L[i] = q
            G = mol.get_G()
            tmp = 0.0
            for k in self.group_list:
#                print i, k, self.get_s(k,i)/self.eta[k], self.theta[k]
                tmp += self.theta[k]*self.get_s(k,i)/self.eta[k] 
                if (k in G):
                    tmp -= G[k]*log(self.get_s(k,i)/self.eta[k])
#                    print i, k, G[k]
            term[i] = tmp

        for i in self.molecule_dict:
            J[i] /= sum_J
            L[i] /= sum_L
#            print i, J[i], L[i]

        ln_gamma_C = {}
        ln_gamma_R = {}
        ln_gamma = {}
        for i in self.molecule_dict:
            mol = self.molecule_dict[i]
            qi = mol.get_q() 
            ln_gamma_C[i] = 1.0 - J[i] + log(J[i]) \
                            - 5.0*qi*(1.0-J[i]/L[i]+log(J[i]/L[i]))
            ln_gamma_R[i] = qi*(1.0-log(L[i])) - term[i]
            ln_gamma[i] = ln_gamma_C[i] + ln_gamma_R[i]
#            print i, ln_gamma_R[i]
#            print self.x[i], J[i], L[i], term[i], ln_gamma_R[i], self.x[i]
        return ln_gamma

    def chemicalPotential(self):
        mu = {}
        mu_ref = {}
        for i in self.molecule_dict:
            pvap = self.molecule_dict[i].vaporPressure(self.T)
            mu_ref[i] = R*self.T*log(pvap/self.p)
            mu[i] = mu_ref[i] + R*self.T*log(self.x[i])
        return mu

