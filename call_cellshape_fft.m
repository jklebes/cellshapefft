addpath("..")

param = struct();
param.pathin = '\\cjw-t640.lifesci.dundee.ac.uk\md3400big\analysis\expk0699_ACTM1\0003_merging_actin\results\merge0000\';
param.contour = [];%'\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0011_mask\mask';
param.pathout = '\\cjw-t640.lifesci.dundee.ac.uk\md3400big\analysis\expk0699_ACTM1\0004_cellshapefft_actin\';
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
param.time_avg = 4;
param.method = []; %choose between ellipse-fitting and inertia matrix
                        %analysis (see paper). "ellipse", "matrix", or "radon"
param.workers = 16;
param.chunk_size = 48; %ideally multiple of number of workers

param.tile_size = []; %size of tiles.
param.overlap= [];
param.cut = []; %masking a circle in the middle of Fourier spectrum
param.propor = []; 
param.stripe_sigma=4; %blurring real space image for stripe removal 
param.sigma = [];
param.strel = [];

%scale, color of lines to draw in visualization:
param.scale = 300;
param.col='green';

param.register = []; %to choose whether spectrum data is saved
param.regsize = []; %??

obj = cellshapefft(param);
obj.full_analysis;
% obj.full_analysis;