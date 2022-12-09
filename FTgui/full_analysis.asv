function Results = full_analysis(Param)

    %% Initialization of the matrices we will use

    % a structure with info about each image - path, filename, date, etc
    % what does this variable name mean?
    % TODO what 's a better name?
    Anamesi = cleanfolder(Param.name{1}); 

    Param.siz = size_image(Param.name{1},Anamesi(1).name);

    Results = struct('regl',[],'Posi',[],'im_regav',[],'ci',[],'numX',0,'numY',0, 'tile_slopes',0); % result structure

    % Find the positions we will use
    [Results.Posi,Results.regl,Results.numX,Results.numY] = Position_gui(Param);
    
    % Length of all the matrices we will later use
    Len = numel(Param.tsart:Param.tleng-Param.timestep+1);
    Results.ci = Len;
    
    Results.im_regav = struct('spect',zeros(Len,Results.regl,Param.subwindow_size),'M',zeros(Len,Results.regl,2),'S',...
        zeros(Len,Results.numY,Results.numX),...
        'angS',zeros(Len,Results.numY,Results.numX),'a',zeros(Len,Results.numY,Results.numX),...
        'b',zeros(Len,Results.numY,Results.numX),'phi',zeros(Len,Results.numY,Results.numX));

    %% Computation of all the averaged (FT) spectra on time and sample
    
   Results.im_regav.spect = avfft_forall(Param.name,Results.regl,Param,Results.Posi);
   %% Analysis in size when required
   
    if Param.regsize == 1
        Results = size_analysis(Param,Results);
    else
    end

   %%
   if ~isempty(Param.heightmapfiles)
    Results.tile_slopes = basis_from_heightmap(Param.heightmapfiles, Param.voxel, Results.regl, Results.Posi, Param.subwindow_size, Results.ci);
   end
    
   %% Computation of cell deformation 
   
   if Param.def_ellipse == 1 && Param.regsize == 0
        Results = size_analysis(Param,Results);
        Results = def_analysis_ellips(Results);     % compute cell orientation subimages
   end
   if Param.def_ellipse == 1 && Param.regsize == 1
        Results = def_analysis_ellips(Results);     % compute cell orientation subimages
   else %no ellipse analysis, inertia matrix analysis
       Results = def_analysis(Param,Results);     % compute cell orientation subimages
   end
   
   %% Print result
   
   affich_result(Param,Results,'yellow');        % create maps of cell orientations

%%

    if Param.regsize == 1
          affich_size(Param,Results,'red');
    else
    end
    %% Register results 

    save([Param.pathout,'Results.mat'],'Results','-v7.3');% save results
end












