clear all
tic

global N_bands Lx Ly
N_bands = 1;
Lx = 29; Ly = 29;

%%%
load Si_nanowires_3nm_hse_valence.mat
k_range = 301:2:900;
k_range = 601:2:850;
b_range = 577:576+N_bands;
kpoints_ref = kpoints(k_range,:);
energies_ref = energies(b_range,:);
energies_ref = energies_ref'; 
energies_ref = energies_ref(k_range,:);
%%%% Sort the band eigenvalues 'descend'
[m,~] = size(energies_ref);
for i = 1:m
    energies_ref(i,:) = sort(energies_ref(i,:));
end
figure
kpath = kpoints2kpath(kpoints_ref);
plot(kpath,energies_ref);



%%%
func = @E_K_reduce_001;
options = optimset('Display','iter','TolFun',1e-6,'PlotFcns','optimplotresnorm');
mt = 0.196;
ml = 0.916;
k0 = 0.15*2*3.1415926/5.431;
E0 = 0.562259615-0.0867;
a_init = [mt,ml,k0,E0];
[a,resnorm,residual,exitflag,output] = lsqcurvefit(func,a_init,kpoints_ref,energies_ref,[],[],options)

params = a;
energies_fit = E_K_reduce_001(params,kpoints_ref);
figure;
kpath = kpoints2kpath(kpoints_ref);
plot(kpath,energies_ref,'k',kpath,energies_fit,'r')
title('Lsqcurvefit')
toc