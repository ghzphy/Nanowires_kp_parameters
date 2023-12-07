clear all
tic
kpoints = generate_k([0,0,-3.1415916/2],[0,0,3.1415916/2], 161);

global N_bands Lx Ly
N_bands = 30;
Lx = 39; Ly = 39;

mt = 0.196;
ml = 0.916;
k0 = 0.15*2*3.1415926/5.431;
E0 = 0.562259615;
params = [mt,ml,k0,E0];
%params = [0.2036 0.7679 0.1701 0.5811]; %optimized bulk params
energies = E_K_reduce_001(params,kpoints);

kpath = kpoints2kpath(kpoints);
figure;
%hold on
plot(kpath,energies);
ylim([0,2]);
toc