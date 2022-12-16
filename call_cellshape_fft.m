addpath("..")

param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\cable detection\exp0186\segmentation_frame 1-100\186_0057_segments\';
param.contour = [];
param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\exp0186_segmented\results_', char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];
if ispc
    param.ffmpeg_path = 'C:\Users\jklebes001\Miniconda3\pkgs\ffmpeg-4.2.2-he774522_0\Library\bin\ffmpeg';
    % caution : has to be single quotes, or there will be "text scalar"
    % error in system command
else
    param.ffmpeg_path = [];
end
param.siz = [];
param.time_points = [1:56];
param.chunk_size = [];
param.tleng = [];
param.timestep = [];
param.tile_size = 128; %size of tiles.
param.overlap=0;
param.fres = [];
param.cut = [];
param.propor = [];
param.threshold = []; %in addition to keeping top propor% of points, 
                      %to dealt with varying signal-noise ration ignore 
                      %intensities below this threshold
param.stripe = []; %bool - mask x and y stripe on Fourier spectrum
param.stripewidth = [];
param.sigma = [];
param.scale = [];
param.strel = [];
param.register = [];
param.regsize = []; %??
param.workers = 16;

obj = cellshapefft(param);
obj.full_analysis_chunks;
% obj.full_analysis;