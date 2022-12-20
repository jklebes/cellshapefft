addpath("..")

param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0003_merge\merged_images\merge\';
param.contour = []; %'\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0011_mask\mask';
param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\expk0245\results_', char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];
if ispc
    param.ffmpeg_path = 'C:\Users\jklebes001\Miniconda3\pkgs\ffmpeg-4.2.2-he774522_0\Library\bin\ffmpeg';
    % caution : has to be single quotes, or there will be "text scalar"
    % error in system command
else
    param.ffmpeg_path = [];
end
param.siz = [];
param.time_points = [1:484];
param.chunk_size = 48;
param.tleng = [];
param.timestep = [];
param.tile_size = []; %size of tiles.
param.overlap=[];
param.fres = [];
param.cut = [];
param.propor = [];
param.threshold = []; %in addition to keeping top propor% of points, 
                      %to dealt with varying signal-noise ration ignore 
                      %intensities below this threshold
param.stripe =true; %bool - mask x and y stripes on Fourier spectrum
param.stripewidth = [];
param.sigma = [];
param.scale = [];
param.strel = [];
param.register = []; %to choose whether spectrum data is saved
param.regsize = []; %??
param.workers = 16;

obj = cellshapefft(param);
obj.full_analysis_chunks;
% obj.full_analysis;