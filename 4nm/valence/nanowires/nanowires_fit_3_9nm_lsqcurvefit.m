clear all
tic

global N_bands Lx Ly
N_bands = 1;
Lx = 39; Ly = 39;

%%%
load Si_nanowires_4nm_hse_valence.mat
k_range = 601:2:750;
b_range = 901-N_bands:900;
kpoints_ref = kpoints(k_range,:);
energies_ref = energies(b_range,:);
energies_ref = energies_ref'; 
energies_ref = energies_ref(k_range,:);
%%%% Sort the band eigenvalues 'descend'
[m,~] = size(energies_ref);
for i = 1:m
    energies_ref(i,:) = sort(energies_ref(i,:),'descend');
end
figure
kpath = kpoints2kpath(kpoints_ref);
plot(kpath,energies_ref);



%%%
func = @nanowires_valence;
options = optimset('Display','iter','TolFun',1e-6,'PlotFcns','optimplotresnorm');
%options = optimset('Display','iter','TolFun',1e-8);
%a_init = [-5.53  -3.64  -8.32 -0.562259615];
%a_init = [-5.88  -2.16  -7.26 -0.562259615];
a_init = [-5.1559   -2.5487   -6.5469   -0.5201];
a_int = [-2.2631   -2.0711   -3.4293   -0.5721];
a_int = [-3.5131   -1.4958   -4.3062   -0.5976];
a_int = [-5.88  -2.16  -7.26 -0.562259615];
[a,resnorm,residual,exitflag,output] = lsqcurvefit(func,a_init,kpoints_ref,energies_ref,[],[],options)
%%%

% a = [-5.1559   -2.5487   -6.5469   -0.5201];
% a = [-3.5131   -1.4958   -4.3062   -0.5976];
params = a;
energies_fit = nanowires_valence(params,kpoints_ref);
figure;
kpath = kpoints2kpath(kpoints_ref);
plot(kpath,energies_ref,'k',kpath,energies_fit,'r')
title('Lsqcurvefit')
toc