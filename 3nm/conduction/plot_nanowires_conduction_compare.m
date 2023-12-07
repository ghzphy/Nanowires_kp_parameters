clear all
tic

global N_bands Lx Ly Lz
N_bands = 12; Lx = 29; Ly = 29; Lz = 29;

%%%%%%%%%%%%%%%%%%%%%%%%%% Reference Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Si_nanowires_3nm_hse_valence.mat
kpoints = kpoints(201:2:1000,:);
energies = energies(577:576+N_bands,:); % 573£º576
energies = energies'; 
energies = energies(201:2:1000,:);

h = figure;
%set(h,'visible','off')
%figure('visible','off')
kpath = kpoints2kpath(kpoints);
plot(kpath,energies,'k');
%ylim([-2,-0.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mt = 0.196;
ml = 0.916;
k0 = 0.15*2*3.1415926/5.431;
E0 = 0.562259615;%-0.0867;
%E0 = 0.5621;
params = [mt,ml,k0,E0];

params = [0.2438    0.7803    0.1720    0.5127];

energies = E_K_reduce_001(params,kpoints);
kpath = kpoints2kpath(kpoints);
hold on
plot(kpath,energies,'r');
ylim([0.5,1.2]);
toc