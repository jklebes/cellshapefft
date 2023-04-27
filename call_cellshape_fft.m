addpath(genpath("..")); %folder utilities is expected to be sister directory

%standalong script to run cellshapefft on a directory of files.
%In general in our experiment, use call_analysis_pipeline.m instead to
% create an analysis_object rA and call its cellshapeft function

param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\cable detection\exp0186\focused_one';
param.contour = [];%'\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0011_mask\mask';
param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\cable detection\exp0186\results_actin_', ...
    char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];

%insert your ffmpeg path here if it's not just command 'ffmpeg'
if ispc
    param.ffmpeg_path = 'C:\Users\jklebes001\Miniconda3\pkgs\ffmpeg-4.2.2-he774522_0\Library\bin\ffmpeg';
    % caution : has to be single quotes, or there will be "text scalar"
    % error in system command
else
    param.ffmpeg_path = [];
end

param.timePoints = [1:48];
param.time_avg = 1;
param.method = 'matrix'; %choose between ellipse-fitting and inertia matrix
                  %analysis (see paper). "ellipse", "matrix"
param.workers = 16;
param.chunk_size = 96; %ideally multiple of number of workers
param.tileSize = 100; %size of tiles.
param.overlap= 0.5;
param.cut = 5; %masking a circle in the middle of Fourier spectrum
param.propor = 0.02; 
param.sigma = 0.8;
param.strel = 4;
%scale, color of lines to draw in visualization:
param.scale = [];
param.col='green'; 

param.writeSpectra = []; %to choose whether spectrum data is saved
param.stretchedNoise = true; %correction particular to interpolated DSLM data


obj = cellshapefft(param);
obj.full_analysis();