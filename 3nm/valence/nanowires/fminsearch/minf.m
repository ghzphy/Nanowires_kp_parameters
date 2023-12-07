function f = minf(params,kpoints_ref, energies_ref)
%global N_bands Lx Ly
[m,n] = size(energies_ref);
energies_fit = nanowires_valence(params,kpoints_ref);
f = 0;
for i = 1:m
    for j = 1:n
        f = f + (energies_fit(i,j)-energies_ref(i,j))^2;
    end
end
end