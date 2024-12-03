from model_cas import Fit_KP
import torch
from torch import nn
import numpy as np
import scipy.io as scio


data = scio.loadmat("matlab.mat")

# Energy bands data as references 
references = data['energies_ref']
kvectors = data['kpoints_ref']
# print(references)

def adjust_lr(optimizer, lr):
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr


# Hyperparameters
lr = 0.1
criterion = nn.MSELoss()
threshold = 1e-5
max_training_steps = 10000


def fitting(Optimizer, loss_threshold, max_train_steps, fit_kp):
    train_steps = 1
    reproduces = fit_kp.compute_bands() 
    loss = criterion(reproduces, fit_kp.references)

    print("initial loss value: %.8f" % loss)
    while (loss > loss_threshold) and (train_steps < max_train_steps):

        reproduces = fit_kp.compute_bands()
        loss = criterion(reproduces, fit_kp.references)
        Optimizer.zero_grad()
        loss.backward(retain_graph=True)
        Optimizer.step()
        fit_kp.reinitialize()

        train_steps += 1
        if train_steps % 10 == 0:
            print("steps:%s, loss:%.8f, lr: %.3f" % (train_steps, loss, Optimizer.param_groups[0]['lr']))


def main():
    fit_kp = Fit_KP()
    fit_kp.read_training_set(references, kvectors)
    fit_kp.reinitialize()

    Optimizer = torch.optim.Adam([fit_kp.l, fit_kp.m, fit_kp.n, fit_kp.E0], lr=lr)
    
    fitting(Optimizer, threshold , max_training_steps, fit_kp)


if __name__ == '__main__':
    main()
