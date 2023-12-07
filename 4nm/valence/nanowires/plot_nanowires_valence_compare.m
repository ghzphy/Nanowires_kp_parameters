clear all
tic
set(groot, ...
'DefaultTextInterpreter', 'LaTeX',...
'DefaultAxesTickLabelInterpreter', 'LaTeX',...
'DefaultLegendInterpreter', 'LaTeX',...
'DefaultAxesFontName', 'LaTeX',...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 1.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 20, ...
'DefaultAxesFontName', 'Times New Roman', ...
'DefaultLineLineWidth', 2.0, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 20, ...
'DefaultTextFontName', 'Times New Roman', ...
'DefaultAxesBox', 'on', ...
'defaultfigureposition',[50 50 400 500],...
'defaultAxesTitleFontSize',.3);


global N_bands Lx Ly Lz
N_bands = 12; Lx = 39; Ly = 39; Lz = 39;

%%%%%%%%%%%%%%%%%%%%%%%%%% Reference Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Si_nanowires_4nm_hse_valence.mat
kpoints = kpoints(201:1000,:);
energies = energies(563:900,:); % 573£º576
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
%params = [-5.0171   -2.4918   -6.3461   -0.5245];
%params = [-3.5026   -1.3352   -4.5891   -0.5995];
%params = [-2.5365   -3.0109   -4.6064   -0.5578];
% params = [-3.7414   -1.3122   -4.2528   -0.6068];
% params = [-3.5091   -1.2730   -4.0057   -0.6097];
params = [ -3.0377   -1.0342   -3.2740   -0.6218];

energies = E_K_reduce_100(params,kpoints);

kpath = kpoints2kpath(kpoints);
hold on
plot(kpath,energies,'r');
xlim([kpath(1),kpath(end)])
ylim([-0.9,-0.6])
toc