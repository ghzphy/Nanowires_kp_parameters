import torch
import torch.nn
import numpy as np
import os 
import re

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


def delta(m,n):
    if m == n:
        return 1
    return 0


def delta_odd(m,n):
    if (m + n) %2 == 1:
        return 1
    return 0


class Fit_KP():
    def __init__(self):
        """
        initialize Fit_KP
        """
        super(TBHCNN, self).__init__()
        self.E0 = torch.tensor(-0.562259615384615, requires_grad=True)
        self.aa = 3.80993
        self.l = torch.tensor(-5.88, requires_grad=True)
        self.m = torch.tensor(-2.16, requires_grad=True)
        self.n = torch.tensor(-7.26, requires_grad=True)
        self.KL = 100
        self.KX = 39
        self.KY = 39

        dis = [0 for i in range(self.KX * self.KY)]
        num = []
        for j in range(self.KX):
            for k in range(self.KY):
                index = j * self.KX + k
                dis[index] = (j ** 2 + k ** 2) ** 0.5
                num.append([j + 1, k + 1])
        sort_id = sorted(range(len(dis)), key=lambda x: dis[x], reverse=False)
        num_new = [[] for i in range(self.KX * self.KY)]
        for i in range(len(num)):
            num_new[i] = num[sort_id[i]]
        self.num_new = num_new[:self.KL]

    def read_training_set(self, references, kvectors):
        numk, numb = references.shape
        self.numb = numb
        self.numk = numk
        self.H_size_init = numb
        self.H_size = numb
        # self.references = torch.tensor(references, torch.float64)
        self.references = torch.tensor(references)
        self.K = kvectors


    def reinitialize(self):
        """
        initialize/reinitialize the real-space Hamiltonian matrices
    
        """
        self.L = self.aa * self.l
        self.M = self.aa * self.m
        self.N = self.aa * self.n

        dx = 1
        self.Ham00 = torch.zeros((3 * self.KL, 3 * self.KL))
        self.Ham01 = torch.zeros((3 * self.KL, 3 * self.KL))

        for k in range(self.KL):
            p = self.num_new[k][0]
            q = self.num_new[k][1]
            kp = p * np.pi / self.KX
            kq = q * np.pi / self.KY
            for k1 in range(self.KL):
                p1 = self.num_new[k1][0]
                q1 = self.num_new[k1][1]
                kp1 = p1 * np.pi / self.KX
                kq1 = q1 * np.pi / self.KY
                # H11, H22, H33
                H00_11 = (self.L * 2 / (dx ** 2) + self.M * kp ** 2 + self.M * kq ** 2 + self.E0) * delta(p, p1) * delta(q, q1)
                H00_22 = (self.M * 2 / (dx ** 2) + self.L * kp ** 2 + self.M * kq ** 2 + self.E0) * delta(p, p1) * delta(q, q1)
                H00_33 = (self.M * 2 / (dx ** 2) + self.M * kp ** 2 + self.L * kq ** 2 + self.E0) * delta(p, p1) * delta(q, q1)
                # H12=N*kx*ky; H13=N*kx*kz;
                H00_12 = 0
                H00_21 = H00_12
                H00_13 = 0
                H00_31 = H00_13
                # H23=N*ky*kz
                if (q == q1) or (p == p1):
                    H00_23 = 0
                else:
                    H00_23 = -self.N * (4 * kp1 / np.pi * p / (p ** 2 - p1 ** 2) * delta_odd(p, p1)) * (
                    4 * kq1 / np.pi * q / (q ** 2 - q1 ** 2) * delta_odd(q, q1))
                H00_32 = H00_23
                # Ham
                self.Ham00[3 * k, 3 * k] = H00_11
                self.Ham00[3 * k+1, 3 * k+1] = H00_22
                self.Ham00[3 * k+2, 3 * k+2] = H00_33
                self.Ham00[3 * k+1, 3 * k+2] = H00_23
                self.Ham00[3 * k+2, 3 * k+1] = H00_32


                H01_11 = -self.L / (dx ** 2) * delta(p, p1) * delta(q, q1)
                H01_22 = -self.M / (dx ** 2) * delta(p, p1) * delta(q, q1)
                H01_33 = -self.M / (dx ** 2) * delta(p, p1) * delta(q, q1)
                # H12=N*kx*ky; H13=N*kx*kz;
                if p != p1:
                    H01_12 = -self.N / (2 * dx) * (4 * kp1 / np.pi * p / (p ** 2 - p1 ** 2) * delta_odd(p, p1)) * delta(q,
                                                                                                                   q1)
                else:
                    H01_12 = 0
                H01_21 = H01_12

                if q != q1:
                    H01_13 = -self.N / (2 * dx) * (4 * kq1 / np.pi * q / (q ** 2 - q1 ** 2) * delta_odd(q, q1)) * delta(p,
                                                                                                                   p1)
                else:
                    H01_13 = 0
                H01_31 = H01_13
                # H23=N*ky*kz
                H01_23 = 0
                H01_32 = H01_23

                self.Ham01[3 * k, 3 * k] = H01_11
                self.Ham01[3 * k + 1, 3 * k + 1] = H01_22
                self.Ham01[3 * k + 2, 3 * k + 2] = H01_33
                self.Ham01[3 * k, 3 * k + 1] = H01_12
                self.Ham01[3 * k + 1, 3 * k] = H01_21
                self.Ham01[3 * k, 3 * k + 2] = H01_13
                self.Ham01[3 * k +2, 3 * k] = H01_31

        self.Ham10 = self.Ham01.T


    def compute_bands(self):
        """
        Using the real-space Hamiltonians considered to compute the k-space 
        Hamiltonians to compute the energy bands
        
        """
        self.reproduces = torch.zeros([self.numk, self.H_size], dtype=torch.float64)
        for i in range(self.numk):
            HH = self.Ham00 + np.exp(1j * self.K[i, 2]) * self.Ham01 + np.exp(-1j * self.K[i, 2]) * self.Ham10
            a = torch.linalg.eigvals(HH).type(torch.float64)
            sorted_a, _ = torch.sort(a, descending=True)
            self.reproduces[i] = sorted_a[:self.H_size]
        
        return self.reproduces
