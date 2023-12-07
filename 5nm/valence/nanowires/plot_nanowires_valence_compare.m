clear all
tic

global N_bands Lx Ly Lz
N_bands = 12; Lx = 49; Ly = 49; Lz = 49;

%%%%%%%%%%%%%%%%%%%%%%%%%% Reference Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Si_nanowires_5nm_hse_valence.mat
kpoints = kpoints(201:1000,:);
energies = energies(1445-N_bands:1444,:); % 573£º576
energies = energies'; 
energies = energies(201:1000,:);

h = figure;
%set(h,'visible','off')
%figure('visible','off')
kpath = kpoints2kpath(kpoints);
plot(kpath,energies,'k');
ylim([-2,-0.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kpoints = generate_k([0,0,0],[0,0,3.1415916/3], 61);

L = -5.53;M = -3.64;N = -8.32;E0 = -0.562259615;

params = [L,M,N,E0];
% % Band1 601:2:800
% params = [-2.4809   -2.8570   -4.4130   -0.5632];
% % Band1 601:2:780
% params = [-3.8846   -1.1828   -4.1587   -0.5935];
% % Band1 601:2:770
% params = [ -3.6666   -1.2130   -4.0277   -0.5937];
% % Band3 601:2:800
% params = [-2.4590   -2.3970   -3.9901   -0.5717];
% % Band3 601:2:780
% params = [-2.5230   -2.4426   -4.1190   -0.5707];
% 
% params = [-3.5030   -1.2676   -3.9899   -0.5931];

energies = E_K_reduce_100(params,kpoints);

kpath = kpoints2kpath(kpoints);
hold on
%figure('visible','off')
plot(kpath,energies,'r');
ylim([-0.9,-0.6])
toc