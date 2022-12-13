param = struct();
param.pathin = pathConverter('Y:\mDrives\storage4\raid0_backup\expk\expk0245\0003_merge\merged_images\merge');
param.contour = [];%pathConverter('Y:\mDrives\storage4\raid0_backup\expk\expk0245\0011_mask\mask');
param.pathout = [pwd filesep 'results_parallel'];
param.siz = [];
param.time_points = [1:220];
param.chunk_size = [];
param.tleng = [];
param.timestep = [];
param.pas2 = [];
param.fres = [];
param.cut = [];
param.propor = [];
param.sigma = [];
param.scale = [];
param.strel = [];
param.register = [];
param.regsize = [];
param.workers = 20;

obj = cellshapefft(param);
obj.full_analysis_chunks;
% obj.full_analysis;

