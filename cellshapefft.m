classdef cellshapefft < handle
    %{
    Original code by Durande (2017)
    
    Further modification by Guillermo Serrano Najera (2019)
     > Converted into a class for convenience
     > Plotting functions modified to be faster with very large images
     > Parallelization def_analysis and plotting
     > Array multiplication (.*) for masking (much faster than original)

    Modifications Jason Klebes (2022)
     > splitting off some functions (periodicDecomposition) for re-use.
     > Should work as installed from github with sister projects (+utilities).
     > paralellized differently, each worker should run through all steps
     independently.  
         > Had to undo some of the one-file, one-object structure and make
         function files, for parallelization reasons.
     > removed function variants, in particular made parallel variants the
     only and main variant
     > added correction for stretched noise, specific to our microscopy setup
    
    Code to analyse the cell deformation on large images with Fourier 
    transform. This code generates maps of cell deformation based on the 
    computation of the inertia matrix on Fourier transforms of subimages.
    
    Known issues:
     > When using mask for computations (not only for plotting) the x0, y0
     positions of the subplots can be different between images. At the 
     moment it is better to run the program without mask and apply the mask 
     only for plotting.
    
    To run the code:
        1. Create param struct. Empty values will use the defaults.
        2. Create cellshapefft object
        3. Run analysis by calling method "full_analysis"
    
    Example:
    
    %TODO check current

    param = struct();
    param.pathin = 'C:PATH\TO\IMAGES';
    param.contour = 'C:PATH\TO\MASK';
    param.pathout = 'C:PATH\TO\RESULTS';
    param.siz = [];
    param.timePoints = [];
    param.chunk_size = [];
    param.tleng = [];
    param.timestep = [];
    param.tileSize = [];
    param.cut = [];
    param.propor = [];
    param.sigma = [];
    param.scale = [];
    param.strel = [];
    param.register = [];
    param.regsize = [];
    param.workers = [];

    obj = cellshapefft(param);
    obj.full_analysis_chunks;
    %}
    properties
        param
        im_regav
        tileCoords
        nTiles
    end

    methods

        function obj = cellshapefft(param)
            %If no contour is needed, enter param.contour = [];
            param.reader = expReader(param.pathin);
            param.contour_reader = expReader(param.contour);

            if isempty(param.timePoints) % indices of images to analyze
                param.timePoints = param.reader.timePoints;
            end

            param.tleng = size(param.timePoints,1);
            
            if isempty(param.time_avg) %sliding average over how many time points
                param.time_avg = 1; %equivalent to no averaging
            end

            if isempty(param.chunk_size) % indices of images to analyze
                param.chunk_size = 96;
            end

            if isempty(param.tileSize)
                param.tileSize = 128; % size of subimages
            end

            if isempty(param.overlap) % indices of images to analyze
                param.overlap = 0.5;
            end

            if isempty(param.scale)
                param.scale = 150;  % scaling for cell deformation visualization
            end

            if isempty(param.col)
                param.col = 'yellow';  % scaling for cell deformation visualization
            end

            if isempty(param.propor)
                param.propor = 0.02;  % value of the percentile of points to keep
            end

            if isempty(param.cut)
                param.cut = 2;  % size of the cut around the center of the Fourier spectrum
            end

            if isempty(param.strel)
                param.strel = 4;  % strel value in the smoothing filter
            end

            if isempty(param.sigma)
                param.sigma = 0.8;  % sigma value in the gaussian smoothing kernel
            end

            if isempty(param.register)
                param.register = 0;  %
            end

            if isempty(param.regsize)
                param.regsize = 0;  %
            end

            if isempty(param.method)
                param.method = "matrix";
            end

            if isempty(param.ffmpeg_path)
                param.ffmpeg_path = 'ffmpeg';
            end

            if isempty(param.workers)
                param.workers = 4;  % cores for parallel computation
            end


            obj.param = param;
            mkdir(obj.param.pathout);

            %save a copy of input files
            copyfile(which('spectrum_analysis'),[obj.param.pathout filesep 'spectrum_analysis_' datestr(datetime,'yymmdd_HHMM') '.m']);
            if obj.param.method=="matrix"
                copyfile(which('inertia_matp_sigma'),[obj.param.pathout filesep 'inertia_matp_sigma_' datestr(datetime,'yymmdd_HHMM') '.m']);
                copyfile(which('S_angS'),[obj.param.pathout filesep 'S_angS_' datestr(datetime,'yymmdd_HHMM') '.m']);
            elseif obj.param.method=="ellipse"
                copyfile(which('size_def'),[obj.param.pathout filesep 'size_def_' datestr(datetime,'yymmdd_HHMM') '.m']);
            
            end

            info = imfinfo(obj.param.reader.fileName); % get information about the images
            obj.param.siz=[info.Width info.Height];        % extract size of the images
            %obj.results = struct('nTiles',[],'tileCoords',[],'im_regav',[],'ci',[]); % result structure
            obj.tileCoords=[];
            obj.im_regav =[];
        end

        function full_analysis(obj)
            % This function is the main function of the program, it uses the parameters
            % determined in the preliminary analysis to perform the functions
            %
            % # Position
            % # spectrum_analysis
            % # def_analysis
            % # affich_result
            %
            % and register the results in a matrix.
            %
            %--------------------------------------------------------------------------
            %
            % *INPUT*: + obj.param for the analysis in Fourier transform, obj.param is a structure
            %        initialized in FullGUI_analysis_chicken.m and that contains fields
            %        for 'siz': the size of the image in pixels
            %        'tleng': the number of images to treat
            %        'pathin': the path to the folder containing the folder with the
            %        images
            %        'pathout': path to the folder that will contain the results
            %        'names': name of the folder with the images
            %        'fres': name of the folder with the results images
            %        'timestep': number of images to perform the average on
            %        'tileSize': the size of the subimages
            %        'sigma': sigma value in the gaussian soothing kernel
            %        'cut': size of the cut around the center of the Fourier spectrum
            %        'proportion': value of the percentile of points to keep
            %        'strel': strel value in the soothing filter
            %
            % *OUTPUT*: modifies the properties obj.tileCoords and obj.im_regav:
            %          'tileCoords': a matrix of size (nTiles x 2) in case with no mask
            %          or (nTiles x 2 x timepoints) in case with contour mask, that
            %          contains the positions of the center of the subimages at each
            %          times.
            %          'im_regav: Structure array containing, for each timepoint, 
            %          a structure 
            %          containing fields 'S' for the amplitude of the cell deformation,
            %          'angS' for the angle of the cell deformation, 'spect' for the
            %          average spectrum, 'a' for the half major axis of the ellipse in
            %          real space, 'b' the half minor axis of the ellipse in real
            %          space, phi the angle of the ellipse 
            %          for each tile.
            %
            %--------------------------------------------------------------------------
            pool = gcp('nocreate');
            if isempty(pool)
                parpool(obj.param.workers);
            else
                if pool.NumWorkers ~= obj.param.workers
                    delete(gcp('nocreate'))
                    parpool(obj.param.workers);
                end
            end

          
            %once per run
            obj.Position();               % register positions of each subimages
            %TODO could make part of parallel loop, esp in case with mask

            spectsize=ceil(obj.param.tileSize/4);
            %obj.im_regav = repmat(struct(repmat(struct('S',[]...
            %    ,'angS',[],'a',[],'b',[],'phi',[]),max(obj.nTiles),1)),1,obj.param.tleng);
            % Initialization of the structure im_regav, containing spectrums
            % cell deformations and sizes for each averaged subimage
            %instead 2 or 3 arrays to hold results
            resultsdims = [max(obj.nTiles), obj.param.tleng];

            %from affich - prepare diretory to save images in
            out_folder = [obj.param.pathout filesep 'img'];
            mkdir(out_folder);

            time_vector = obj.param.timePoints;
            quality = zeros([resultsdims]);
            switch obj.param.method
                case "matrix"
                    %create variables
                    inertia_matrix = zeros([2 2 resultsdims]);
                    SangS = zeros([2 resultsdims]);
                    write_spectra=false;

                    %chunked
                    for t_ind = time_vector(1:obj.param.chunk_size:end)
                        time_lims = [t_ind min(t_ind+obj.param.chunk_size-1, time_vector(end))];
                        if write_spectra
                            spectra_chunk=zeros([spectsize, spectsize, resultsdims(1), obj.param.chunk_size]);
                        end
                        param=obj.param; %take a copy of the structs to send to each worker
                        tileCoords = obj.tileCoords; %intentional, ignore the yellow suggestions
                       
                        parfor c = time_lims(1):time_lims(2) %main parfor (substitute 'for' for debugging only)
                            %per timepoint
                            if ~isempty(param.contour)
                                Posi_c=tileCoords(:,:,c);
                            else
                                Posi_c=tileCoords;
                            end
                            nTiles = size(Posi_c,1);
                        %inertia matrix route
                        data=param.reader.stack(c);
                        spectra=spectrum_analysis(param, data, nTiles, Posi_c, spectsize, c);        % compute spectrums of subimages
                        %here save spectra in chunks if applicable
                        if write_spectra
                            spectra_chunk(:,:,:,c)=spectra;
                        end
                        [ quality(:, c),inertia_matrix(:,:,:,c)]= inertia_matp_sigma(param, spectra, nTiles);  % compute cell orientation subimages
                    end %first parfor
                    if write_spectra
                        save([obj.param.pathout filesep 'spectra_' num2str(time_lims(1),'%04d') '_' num2str(time_lims(2),'%04d') '.mat'], 'spectra_chunk', '-v7.3');
                    end
                    end
                    %time average
                    inertia_matrix=movmean(inertia_matrix, obj.param.time_avg, 4);
                    quality=movmean(quality, obj.param.time_avg, 2);
                    
                    %transform quality scores to range 0 to 1 for
                    %visualization
                    quality_rescaled=rescale(quality);

                    for t_ind = time_vector(1:obj.param.chunk_size:end) %another chunked process
                        time_lims = [t_ind min(t_ind+obj.param.chunk_size-1, time_vector(end))];
                        %parfor to process and visualize
                        param=obj.param; %take a copy of the structs to send to each worker
                        tileCoords = obj.tileCoords;
                        nTiles = obj.nTiles; %intentional, ignore the yellow suggestions
                       
                        parfor c = time_lims(1):time_lims(2) %parfor
                            if ~isempty(param.contour)
                                Posi_c=tileCoords(:,:,c);
                            else
                                Posi_c=tileCoords;
                            end
                            SangS_c = S_angS(inertia_matrix(:,:,:,c), nTiles);
                            visualization_strain(param, SangS_c, quality_rescaled(:, c), Posi_c, nTiles, out_folder, c); % create maps of cell orientations
                            SangS(:,:,c) = SangS_c;
                        end
                        obj.save_chunk(inertia_matrix(:,:,:,time_lims(1):time_lims(2)), 'inertia_matrix', time_lims)    % save results (time-averaged)
                        obj.save_chunk(SangS(:,:,time_lims(1):time_lims(2)),'SangS', time_lims) 
                        obj.save_chunk(quality(:,time_lims(1):time_lims(2)),'quality', time_lims) 
                    end
                    
                    obj.im_regav = struct('S', SangS(1,:,:), 'angS', SangS(2,:,:));
                    
                case "ellipse"
                    abphi = zeros([3 resultsdims]); %create variables
                    SangS = zeros([2 resultsdims]);

                    %parfor and save chunks
                    %chunk
                    for t_ind = time_vector(1:obj.param.chunk_size:end)
                        time_lims = [t_ind min(t_ind+obj.param.chunk_size-1, time_vector(end))];

                        param=obj.param; %take a copy of the structs to send to each worker
                        tileCoords = obj.tileCoords; %intentional, ignore the yellow suggestions
                        
                    parfor c = time_lims(1):time_lims(2) %main parfor (substitute 'for' for debugging only)
                        %per timepoint
                        if ~isempty(param.contour)
                            Posi_c=tileCoords(:,:,c);
                        else
                            Posi_c=tileCoords;
                        end
                        nTiles = size(Posi_c,1);
                        %inertia matrix route
                        data=param.reader.stack(c);
                        spectra=spectrum_analysis(param, data, nTiles, Posi_c, spectsize, c);        % compute spectrums of subimages
                        %here save spectra in chunks if applicable
                        [ quality(:, c),abphi(:,:,c)] = deformation_ellipse(param, spectra, nTiles);  % compute cell orientation subimages
                    end
                    end

                    %time average - TODO is this th right ellipse quantities
                    %to avg?
                    abphi=movmean(abphi, obj.param.time_avg, 4);
                    quality=movmean(quality, obj.param.time_avg, 2);
                    
                    %transform quality scores to range 0 to 1 for
                    %visualization
                    quality=rescale(quality);
                    for t_ind = time_vector(1:obj.param.chunk_size:end) %another chunked process
                        time_lims = [t_ind min(t_ind+obj.param.chunk_size-1, time_vector(end))];
                        %parfor to process and visualize
                        param=obj.param; %take a copy of the structs to send to each worker
                        tileCoords = obj.tileCoords;
                        nTiles = obj.nTiles; %intentional, ignore the yellow suggestions
                        parfor c = time_lims(1):time_lims(2) %parfor
                            if ~isempty(param.contour)
                                Posi_c=tileCoords(:,:,c);
                            else
                                Posi_c=tileCoords;
                            end
                            visualization_ellipse(param, abphi, quality(:,c), Posi_c, nTiles, out_folder, c); % create maps of cell orientations
                        end
                        obj.save_chunk(abphi, time_lims)    % save results (time-averaged)
                    end

                    obj.im_regav = struct('a', abphi(1,:,:), 'b', abphi(2,:,:),'phi', abphi(3,:,:));
                case "radon"
                    disp("radon method not yet implemented")
                    exit()
            end

            %once per experiment - make movie
            %eR = expReader([obj.param.pathout]);

            %eR.ffmpeg_path = obj.param.ffmpeg_path;
            %eR.ffmpeg;
        end

        function save_chunk(obj,chunk,name, time_lims)
            save([obj.param.pathout filesep 'cellshapefft_results_'  name '_' num2str(time_lims(1),'%04d') '_' num2str(time_lims(2),'%04d') '.mat'], 'chunk', '-v7.3');% save results
            if ~isempty(obj.param.contour)  
                tileCoords = obj.tileCoords(time_lims);
                save([obj.param.pathout filesep 'cellshapefft_results_' name '_'  num2str(time_lims(1),'%04d') '_' num2str(time_lims(2),'%04d') '.mat'], 'tileCoords', '-v7.3');
            end
        end

        function Position(obj)
            % results = Position(obj.param)
            %
            % This function registers the positions of each subimage at each
            % time, and the number of subimage at each time. The nTiles and tileCoords fields
            % of the Result structure are filled.
            disp('Registering...');

            overlap = floor((1-obj.param.overlap)*obj.param.tileSize);                      % number of shared pixels between
            % subimages

            if isempty(obj.param.contour)                  % if no mask is required
                    x = 1:overlap:(obj.param.siz(1)- obj.param.tileSize+1);
                    y = 1:overlap:(obj.param.siz(2)- obj.param.tileSize+1);
                    [X,Y] = meshgrid(x,y);
                    positions = [X(:) Y(:)];
                    obj.nTiles=size(positions,1) ;
                    obj.tileCoords(1:obj.nTiles,:) =  [X(:) Y(:)];
             % register number of subimages at each time
                %only one value - should be same number of subimages each
                %timepoint

            else                                        % when a mask is required

                obj.nTiles = zeros(obj.param.tleng,1);
                for t = 1:obj.param.tleng
                    b = obj.param.contour_reader.readSpecificImage(obj.param.timePoints(t)); % load mask images
                    b = im2double(b);

                    x = 1:overlap:(obj.param.siz(1)- obj.param.tileSize+1);
                    y = 1:overlap:(obj.param.siz(2)- obj.param.tileSize+1);
                    [X,Y] = meshgrid(x,y);

                    mask = imresize(b, size(X))>=0.5;
                    X = X.*mask;
                    Y = Y.*mask;

                    positions = [X(:) Y(:)];
                    positions(prod(positions,2)==0,:) = [];
                    reg = size(positions,1);
                    obj.tileCoords(1:reg,:,t) = positions;
                    obj.nTiles(t) = reg;

                end
            end
        end


    end
end