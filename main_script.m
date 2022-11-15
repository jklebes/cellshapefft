%% Initialization

clear all 
close all
addpath('FTgui')

Param = struct('siz',[],'tleng',[],'pas2',[],'name',[],'fres',[]...
,'cut',[],'propor',[],'pathin',[],'pathout',[],'timestep',[]...
,'sigma',[],'strel',[],'register',[],'regsize',[],'nbpoints',[],'sizeview',[],'contour',[],'rec',[],'tsart',[]);


%% Output folder location

Param.pathin = {'C:\Users\Melina\Documents\DOCTORAT\Matlab\2A\CAMILLE\Probleme\3\1\'}; 

%% List of folders to treat/folder to treat

%When multiple folders, write them like this:

 Param.name = {'C:\Users\Melina\Documents\DOCTORAT\Matlab\2A\CAMILLE\Probleme\1\';
      'C:\Users\Melina\Documents\DOCTORAT\Matlab\2A\CAMILLE\Probleme\2\'};
% When one folder only is required, write it like this:

 Param.name = { 'C:\Users\Melina\Documents\DOCTORAT\Matlab\2A\CAMILLE\Probleme\2\'};
 
%% Parameters you should change and that are NOT tunable with the user interface

Param.rec = 0.5;                % Overlap between the sliding boxes 
Param.tsart = 1;                % Beginning of the analysis
Param.tleng = 5;               % End of the analysis
Param.timestep = 1;         % Temporal average

%% Parameters that ARE tunable with the user interface

Param.pas2= 100;                % Size of the subwindows
Param.cut = 5;                       % Cut of the spectra
Param.propor =  0.02;          % Proportion of points to keep
Param.sigma = 0.8;              % Gaussian filter
Param.register = 0;             %   Save spectra if 1, 0 otherwise
Param.regsize = 0;              % Analysis in size required ? 
Param.nbpoints = 30;          % Number of points to keep in the ellipse fit
Param.strel = 4;                % Strel of the spectra
Param.def_ellipse = 1;      % Use ellipse method for cell deformation

%% Parameters that can be changed if desired

Param.scale = 150;                    % scaling for cell deformation 
Param.fres = 'results';               % name of the output folder created

%% Never change this 

Param.contour=[];

%% Creation of the output folder

 Param.pathout = [Param.pathin{1}, Param.fres ,'_averageon_', num2str(Param.timestep) ,'_averageon_' num2str(size(Param.name,1)) '/']; % path to result folder
 if ~exist(Param.pathout,'dir')             % creation of result folder
 mkdir(Param.pathout)
 end

 %% Visualization of the parameters
GUI_deformation_ft(Param)  


%% Analysis alone if required

%full_analysis(Param)

%% Save results 

save([Param.pathout,'Param.mat'],'Param','-v7.3');


















