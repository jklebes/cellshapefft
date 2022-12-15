addpath("..")

param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0643_ACTM1\proper\ch1_myosin\0003_merging\results\merge\';
param.contour = [];
param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\expk0643\results_myosin_', char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];
param.ffmpeg_path = 'C:\Users\jklebes001\Miniconda3\pkgs\ffmpeg-4.2.2-he774522_0\Library\bin\ffmpeg';
param.siz = [];
param.time_points = [1:263];
param.chunk_size = [];
param.tleng = [];
param.timestep = [];
param.pas2 = 256; %size of tiles.
param.fres = [];
param.cut = [];
param.propor = [];
param.stripe_sigma = 3; %filter on real space tiles, attempt to get rid of 
                        %a high frequency scan stripe artefact causing
                        %apparent horizontal anisotropy
param.sigma = []; %blurring of Fourier space images
param.scale = [];
param.strel = [];
param.register = [];
param.regsize = [];
param.workers = 16;

obj = cellshapefft(param);
obj.full_analysis_chunks;
% obj.full_analysis;