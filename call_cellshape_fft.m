param = struct();
param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\DSLM_expk\expk0245\0003_merge\merged_images\';
param.contour = [];%pathConverter('Y:\mDrives\storage4\raid0_backup\expk\expk0245\0011_mask\mask');
param.pathout = ['//lfs.lifesci.dundee.ac.uk/lfs/cjw/Jason/expk0245/results_parallel/'];
param.siz = [];
param.time_points = [1:484];
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
param.workers = 12;

obj = cellshapefft(param);
obj.full_analysis_chunks;
% obj.full_analysis;

