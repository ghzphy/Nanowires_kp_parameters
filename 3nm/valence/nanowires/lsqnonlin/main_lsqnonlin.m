clear all
tic

global N_bands Lx Ly weight
N_bands = 5;
Lx = 30; Ly = 30;


%%%
%load SiNW1_3nm_vogl_valence.mat
load SiNW1_3nm_HSE_valence.mat
k_range = 601:2:850;
b_range = 577-N_bands:576;
N_k = length(k_range);
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


b_weight = ones(N_bands,1);
b_weight = [1,1.0,1.0,0.0,0.0];
k_weight = ones(N_k,1);
k_weight = [ones(100,1); zeros(N_k-100,1)];
weight = zeros(N_k,N_bands);
for i_k = 1: N_k
    for i_b = 1:N_bands
        weight(i_k,i_b) = k_weight(i_k)*b_weight(i_b);
    end
end

%%%
func = @(params)compared_energies(params,kpoints_ref,energies_ref);
%a_init = [-5.53  -3.64  -8.32 -0.562259615];
%a_init = [-5.88  -2.16  -7.26 -0.562259615];
params_init = [-5.1559   -2.5487   -6.5469   -0.5201];
params_init = [-5.53  -3.64  -8.32 -0.562259615];
lb = [];
ub = [];
%options = optimset('Display','iter','TolFun',1e-8);
options = optimset('Display','iter','TolFun',1e-8,'PlotFcns','optimplotresnorm');
[params_end,resnorm,residual,exitflag,output] = lsqnonlin(func,params_init,lb,ub,options)
%%%

energies_fit = nanowires_valence(params_end,kpoints_ref);
figure;
kpath = kpoints2kpath(kpoints_ref);
plot(kpath,energies_ref,'k',kpath,energies_fit,'r')
title('Lsqcurvefit')
toc