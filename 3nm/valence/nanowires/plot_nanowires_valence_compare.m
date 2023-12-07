clear all
tic

global N_bands Lx Ly Lz
N_bands = 12; Lx = 29; Ly = 29; Lz = 29;

%%%%%%%%%%%%%%%%%%%%%%%%%% Reference Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Si_nanowires_3nm_hse_valence.mat
kpoints = kpoints(201:1000,:);
energies = energies(563:576,:); % 573£º576
energies = energies'; 
energies = energies(201:1000,:);

h = figure;
kpath = kpoints2kpath(kpoints);
plot(kpath,energies,'k');
ylim([-2,-0.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kpoints = generate_k([0,0,0],[0,0,3.1415916/3], 61);

L = -5.53;M = -3.64;N = -8.32;E0 = -0.562259615;
% TED Data
%L = -6.28;M = -2.26;N = -7.26; E0 = -0.562259615;

params = [L,M,N,E0];
%params = [-5.0171   -2.4918   -6.3461   -0.5245];
%params = [-3.5026   -1.3352   -4.5891   -0.5995];
params = [-3.5131   -1.4958   -4.3062   -0.5976];
%params = [-3.5202   -1.4993   -4.3169   -0.5972];
energies = E_K_reduce_100(params,kpoints);
kpath = kpoints2kpath(kpoints);
hold on
plot(kpath,energies,'r');
ylim([-1.2,-0.6])
%print(gcf,'-dpng','abc.png')  
%saveas(gcf, filename)
toc