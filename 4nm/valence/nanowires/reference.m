function [kpoints,energies,kpoints_ref,energies_ref] = reference()
global N_bands
%N_bands = 2;

load SiNW1_3nm_vogl_valence.mat kpoints  energies;

%%% filter the data
k_range = 601:2:800;
b_range = 577-N_bands:576;


kpoints_ref = kpoints(k_range,:);
energies_ref = energies(b_range,:);
energies_ref = energies_ref'; 
energies_ref = energies_ref(k_range,:);


%%%% Sort the band eigenvalues 'descend'
[m,~] = size(energies_ref);
for i = 1:m
    energies_ref(i,:) = sort(energies_ref(i,:),'descend');
end

% figure
% kpath = kpoints2kpath(kpoints);
% plot(kpath,energies);
% ylim([-1.5,1])
end