addpath("..")

param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0636_ACTM1\proper\ch1_myosin\0003_merging\results\merge';
param.contour = [];%'\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0011_mask\mask';
param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\expk0636\results_myosin_', char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];
if ispc
    param.ffmpeg_path = 'C:\Users\jklebes001\Miniconda3\pkgs\ffmpeg-4.2.2-he774522_0\Library\bin\ffmpeg';
    % caution : has to be single quotes, or there will be "text scalar"
    % error in system command
else
    param.ffmpeg_path = [];
end

param.time_points = [1:220];
param.method = []; %choose between ellipse-fitting and inertia matrix
                        %analysis (see paper). "ellipse", "matrix", or "radon"
param.workers = 16;
param.chunk_size = 48; %ideally multiple of number of workers

param.tile_size = []; %size of tiles.
param.overlap= [];
param.cut = []; %masking a circle in the middle of Fourier spectrum
param.propor = []; 
param.threshold = []; %in addition to keeping top propor% of points, 
                      %to dealt with varying signal-noise ration ignore 
                      %intensities below this threshold
param.stripe =true ; %bool - mask x and y stripes on Fourier spectrum
param.stripewidth = [];
param.sigma = [];
param.stripe_sigma = 3; %blur real space tiles to remove stripe artfact
param.strel = [];
%scale, color of lines to draw in visualization:
param.scale = [];
param.col=[];

param.register = []; %to choose whether spectrum data is saved
param.regsize = []; %??

obj = cellshapefft(param);
obj.full_analysis;
% obj.full_analysis;