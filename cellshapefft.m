classdef cellshapefft < handle
    %{
    Original code by Durande (2017)
    
    Further modification by Guillermo Serrano Najera (2019)
     > Converted into a class for convenience
     > Plotting functions modified to be faster with very large images
     > Parallelization def_analysis and plotting
     > Array multiplication (.*) for masking (much faster than original)

    Modifications Jason Klebes (2022)
     > splitting off some functions (perdecomp) for re-use.
     > Should work as installed from github with sister projects (+utilities).
     > paralellized differently, each worker should run through all steps
     independently.  
         > Time averaging of spectrums taken out (for now?)
         > Had to undo some of the one-file, one-object structure and make
         function files, for parallelization reasons.
     > removed function variants, in particular made parallel variants the
     only and main variant
    
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
    
    param = struct();
    param.pathin = 'C:PATH\TO\IMAGES';
    param.contour = 'C:PATH\TO\MASK';
    param.pathout = 'C:PATH\TO\RESULTS';
    param.siz = [];
    param.time_points = [];
    param.chunk_size = [];
    param.tleng = [];
    param.timestep = [];
    param.tile_size = [];
    param.fres = [];
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
        Posi
        regl
    end

    methods

        function obj = cellshapefft(param)

            import utilities.expReader;
            %If no contour is needed, enter param.contour = [];
            param.reader = expReader(param.pathin);
            param.contour_reader = expReader(param.contour);

            if isempty(param.time_points) % indices of images to analyze
                param.time_points = param.reader.timePoints;
            end

            if isempty(param.chunk_size) % indices of images to analyze
                param.chunk_size = 48;
            end

            if isempty(param.tleng) % number of images to analyze
                param.tleng = length(param.time_points);
            end

            if isempty(param.tile_size)
                param.tile_size = 256; % size of subimages
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
                param.propor = 0.13;  % value of the percentile of points to keep
            end

            if isempty(param.cut)
                param.cut = 2;  % size of the cut around the center of the Fourier spectrum
            end

            if isempty(param.stripe)
                param.stripe = false;  % whether to mask for x- y- stripes
            end

            if isempty(param.stripewidth)
                param.stripewidth = 4;  % default number of fourier space pixels to mask
            end


            if isempty(param.stripe_sigma)
                param.stripe_sigma = [];  % default number of fourier space pixels to mask
            end

            if isempty(param.strel)
                param.strel = 4;  % strel value in the soothing filter
            end

            if isempty(param.threshold)
                param.threshold = 0;  % intensity threshold for FT pixels
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

            if isempty(param.ellipse_fit)
                param.ellipse_fit = false;
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
            try
                if ispc
                    system(['copy /b call_cellshape_fft.m ' obj.param.pathout filesep]);
                    system(['copy /b cellshapefft.m ' obj.param.pathout filesep]);
                else %unix
                    system(['cp ./call_cellshape_fft.m' obj.param.pathout filesep]);
                    system(['cp ./cellshapefft.m' obj.param.pathout filesep]);
                end
            catch
                %if pathing problems, skip
                disp(["Problems encountered while attempting to copy call_cellshape_fft.m and/or cellshapefft.m to" ...
                    obj.param.pathout ", input file and/or code was not saved"])
            end

            info = imfinfo(obj.param.reader.fileName); % get information about the images
            obj.param.siz=[info.Width info.Height];        % extract size of the images
            %obj.results = struct('regl',[],'Posi',[],'im_regav',[],'ci',[]); % result structure
            obj.Posi=[];
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
            %        'tile_size': the size of the subimages
            %        'sigma': sigma value in the gaussian soothing kernel
            %        'cut': size of the cut around the center of the Fourier spectrum
            %        'proportion': value of the percentile of points to keep
            %        'strel': strel value in the soothing filter
            %
            % *OUTPUT*: modifies the properties obj.Posi and obj.im_regav:
            %          'Posi': a matrix of size (regl x 2) in case with no mask
            %          or (regl x 2 x timepoints) in case with contour mask, that
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
            import utilities.expReader;
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

            spectsize=ceil(obj.param.tile_size/4);
            obj.im_regav = repmat(struct(repmat(struct('spect',[],'S',[]...
                ,'angS',[],'a',[],'b',[],'phi',[]),max(obj.regl),1)),1,obj.param.tleng);
            % Initialization of the structure im_regav, containing spectrums
            % cell deformations and sizes for each averaged subimage

            %from affich - prepare diretory to save images in
            out_folder = [obj.param.pathout filesep 'img'];
            mkdir(out_folder);

            time_vector = obj.param.time_points;
            %chunk
            for t_ind = time_vector(1:obj.param.chunk_size:end)
                time_lims = [t_ind min(t_ind+obj.param.chunk_size-1, time_vector(end))];
                
                param=obj.param; %take a copy of the structs to send to each worker
                Posi = obj.Posi; %intentional, ignore the yellow suggestions
                im_regav=obj.im_regav;
                parfor c = time_lims(1):time_lims(2) %parfor (substitute 'for' for debugging only)
                    %per timepoint
                    if ~isempty(param.contour)
                        Posi_c=Posi(:,:,c);
                    else
                        Posi_c=Posi;
                    end
                    regl = size(Posi_c,1);
                    spectra=spectrum_analysis(param, regl, Posi_c, spectsize, c);        % compute spectrums of subimages
                    if ~param.ellipse_fit
                        %inertia matrix route
                        results_tmp=def_analysis(param, spectra, regl, c);  % compute cell orientation subimages
                        affich_result_parallel(param, results_tmp, Posi_c, regl, out_folder, c); % create maps of cell orientations
                    else
                        %ellipse fitting route
                        results_tmp=size_analysis(param, spectra, regl, c);  % compute cell orientation subimages
                        affich_size(param, results_tmp, Posi_c, regl, out_folder, c); % create maps of cell orientations
                    end
                    im_regav(:,c)=results_tmp;
                end
                obj.im_regav = im_regav;
                %once per chunk
                obj.save_chunk(time_lims)    % save results
                %obj.results = struct('regl',[],'Posi',[],'im_regav',[],'ci',[]); % result structure
            end
            %once per experiment
            rR = expReader([obj.param.pathout filesep 'img']);
            rR.ffmpeg_path = obj.param.ffmpeg_path;
            rR.ffmpeg;
        end

        function save_chunk(obj, time_lims)
            im_regav = obj.im_regav(:, time_lims);
            save([obj.param.pathout filesep 'cellshapefft_results_' num2str(time_lims(1),'%04d') '_' num2str(time_lims(2),'%04d') '.mat'], 'im_regav', '-v7.3');% save results
            if ~isempty(obj.param.contour)  
                Posi = obj.Posi(time_lims);
                save([obj.param.pathout filesep 'cellshapefft_results_' num2str(time_lims(1),'%04d') '_' num2str(time_lims(2),'%04d') '.mat'], 'Posi', '-v7.3');
            end
        end

        function Position(obj)
            % results = Position(obj.param)
            %
            % This function registers the positions of each subimage at each
            % time, and the number of subimage at each time. The regl and Posi fields
            % of the Result structure are filled.
            disp('Registering...');

            overlap = floor((1-obj.param.overlap)*obj.param.tile_size);                      % number of shared pixels between
            % subimages

            if isempty(obj.param.contour)                  % if no mask is required
                    x = 1:overlap:(obj.param.siz(1)- obj.param.tile_size+1);
                    y = 1:overlap:(obj.param.siz(2)- obj.param.tile_size+1);
                    [X,Y] = meshgrid(x,y);
                    positions = [X(:) Y(:)];
                    obj.regl=size(positions,1) ;
                    obj.Posi(1:obj.regl,:) =  [X(:) Y(:)];
             % register number of subimages at each time
                %only one value - should be same number of subimages each
                %timepoint

            else                                        % when a mask is required

                obj.regl = zeros(obj.param.tleng,1);
                for t = 1:obj.param.tleng
                    b = obj.param.contour_reader.readSpecificImage(obj.param.time_points(t)); % load mask images
                    b = im2double(b);

                    x = 1:overlap:(obj.param.siz(1)- obj.param.tile_size+1);
                    y = 1:overlap:(obj.param.siz(2)- obj.param.tile_size+1);
                    [X,Y] = meshgrid(x,y);

                    mask = imresize(b, size(X))>=0.5;
                    X = X.*mask;
                    Y = Y.*mask;

                    positions = [X(:) Y(:)];
                    positions(prod(positions,2)==0,:) = [];
                    reg = size(positions,1);
                    obj.Posi(1:reg,:,t) = positions;
                    obj.regl(t) = reg;

                end
            end
        end


    end
end