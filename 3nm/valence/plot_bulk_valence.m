clear all;
tic

L = [0.5745, 0.5745, 0.5745];
G = [0,0,0];
K = [0.8618, 0.8618, 0.0000];
len_L = 300;
kpoints1 = generate_k(L,G,len_L);
kpoints2 = generate_k(G,K,len_L);
kpoints = [kpoints1;kpoints2];


L = -5.53;
M = -3.64;
N = -8.32;
E0 = -0.562259615;
params = [L,M,N,E0];
energies = bulk_valence(params,kpoints);

kpath = kpoints2kpath(kpoints);
figure;
plot(kpath,energies);
%ylim([-5,0]);
toc