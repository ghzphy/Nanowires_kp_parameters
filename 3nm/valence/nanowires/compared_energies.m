% Input : params  : [L,M,N,E0]
%         kpoints : [:,3]
%         energies_ref: ab-initio calculations
% Output: lsq: [:,N_E] Only return the difference between fit and DFT.
function lsq = compared_energies(params,kpoints,energies_ref)
   global weight;
   energies = nanowires_valence(params,kpoints);
   lsq = energies - energies_ref;
   lsq = lsq.*weight;
end