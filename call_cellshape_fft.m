addpath("..")

param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0003_merge\merged_images\merge\';
param.contour = [];
param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\expk0245\results_', char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];
param.ffmpeg_path = 'C:\Users\jklebes001\Miniconda3\pkgs\ffmpeg-4.2.2-he774522_0\Library\bin\ffmpeg';
param.siz = [];
param.time_points = [1:208];
param.chunk_size = [];
param.tleng = [];
param.timestep = [];
param.tile_size = 256; %size of tiles.
param.fres = [];
param.cut = [];
param.propor = [];
param.threshold = []; %in addition to keeping top propor% of points, 
                      %to dealt with varying signal-noise ration ignore 
                      %intensities below this threshold
param.stripe = true; %bool - mask x and y stripe on Fourier spectrum
param.stripewidth = [];
param.sigma = 0; %try not doing Fourier space Gaussian blur
param.scale = [];
param.strel = [];
param.register = [];
param.regsize = []; %??
param.workers = 16;

obj = cellshapefft(param);
obj.full_analysis_chunks;
% obj.full_analysis;