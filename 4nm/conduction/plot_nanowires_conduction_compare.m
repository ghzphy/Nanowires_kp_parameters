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
kpoints = kpoints(201:2:1000,:);
energies = energies(901:900+N_bands,:); % 573��576
energies = energies'; 
energies = energies(201:2:1000,:);

h = figure;
%set(h,'visible','off')
%figure('visible','off')
kpath = kpoints2kpath(kpoints);
plot(kpath,energies,'k');
%ylim([-2,-0.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kpoints = generate_k([0,0,0],[0,0,3.1415916/3], 61);

mt = 0.196;
ml = 0.916;
k0 = 0.15*2*3.1415926/5.431;
E0 = 0.562259615-0.0188;
params = [mt,ml,k0,E0];
params = [0.2182    0.8580    0.1725    0.5531];

%energies = nanowires_valence(params,kpoints);
energies = E_K_reduce_001(params,kpoints);
kpath = kpoints2kpath(kpoints);
hold on
plot(kpath,energies,'r');
ylim([0.5,1.2]);
toc