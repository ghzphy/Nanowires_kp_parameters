load Si_bulk_LGK.mat
kpoints_origin = kpoints;
kpoints = kpoints(580:620,:);
energies = energies';
energies_origin = energies; 
energies = energies(580:620,2:4);


kpath = kpoints2kpath(kpoints_origin);
figure;
plot(kpath,energies_origin);
ylim([-1,0])

kpath = kpoints2kpath(kpoints);
figure;
plot(kpath,energies,'k.');

func = @bulk_valence;
options = optimset('Display','iter','TolFun',1e-8);
L = -5.53;M = -3.64;N = -8.32;
E0 = -0.562259615;
E0 = -0.551355;
[a,resnorm] = lsqcurvefit(func,[L,M,N,E0],kpoints,energies,[],[],options)
%  -5.4351   -3.7377   -8.4889   -0.5514 resnorm 1.9539e-07
energies_fitted = bulk_valence(a,kpoints_origin);
kpath = kpoints2kpath(kpoints_origin);
figure;
plot(kpath,energies_origin,'k',kpath,energies_fitted,'r')
ylim([-1,0])