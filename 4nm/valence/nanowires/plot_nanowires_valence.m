clear all
tic

global N_bands Ly Lz
N_bands = 40;
Ly = 29; Lz = 29;
%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoints = generate_k([0,0,-3.1415916/8],[0,0,3.1415916/8], 201);

L = -5.53;M = -3.64;N = -8.32;E0 = -0.562259615;
% L = -5.88; M = -2.16; N = -7.26;E0 = -0.562259615;
% L = -6.28;M = -2.26;N = -7.26; E0 = -0.562259615;
% L = -5.4351; M = -3.7377; N = -8.4889; E0 = -0.5514;

params = [L,M,N,E0];
%params = [-5.0171   -2.4918   -6.3461   -0.5245];
%params = [-5.8983   -2.9748   -7.7642   -0.4851];
energies = E_K_reduce_100(params,kpoints);
%energies = E_K_reduce_110(params,kpoints);

kpath = kpoints2kpath(kpoints);
figure;
plot(kpath,energies,'r');
ylim([-1,-0.5]);
title(['L:', num2str(L),' M:', num2str(M),' N:', num2str(N)])
toc