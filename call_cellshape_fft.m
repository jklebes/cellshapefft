addpath(genpath("..")); %folder utilities is expected to be sister directory

%standalong script to run cellshapefft on a directory of files.
%In general in our experiment, use call_analysis_pipeline.m instead to
% create an analysis_object rA and call its cellshapeft function

param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\Confocal\myosin';
param.contour = [];%'\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0011_mask\mask';
param.pathout = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\Confocal\';
%param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\expk0699\results_actin_', char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];

%insert your ffmpeg path here if it's not just command 'ffmpeg'
if ispc
    param.ffmpeg_path = 'C:\Users\jklebes001\Miniconda3\pkgs\ffmpeg-4.2.2-he774522_0\Library\bin\ffmpeg';
    % caution : has to be single quotes, or there will be "text scalar"
    % error in system command
else
    param.ffmpeg_path = [];
end

param.time_points = [];
param.time_avg = 1;
param.method = []; %choose between ellipse-fitting and inertia matrix
                   %analysis (see paper). "ellipse", "matrix", or "radon"
param.workers = 16;
param.chunk_size = 96; %ideally multiple of number of workers

param.tileSize = []; %size of tiles.
param.overlap= [];
param.cut = []; %masking a circle in the middle of Fourier spectrum
param.propor = []; 
param.sigma = [];
param.strel = [];

%scale, color of lines to draw in visualization:
param.scale = 300;
param.col='green';

param.writeSpectra = []; %to choose whether spectrum data is saved

obj = cellshapefft(param);
obj.full_analysis;
% obj.full_analysis;