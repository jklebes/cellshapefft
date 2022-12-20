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
         > Time averaging of spectrums taken out temporarily.
         > Had to undo some of the one-file, one-object structure and make
         function files, for parallelization reasons.
     > removed function variants
     > notebooks, new GUI app for calibrating
    
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
        3. Run analysis by calling method "full_analysis" or "full_analysis_chunks"
    
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

            if isempty(param.propor)
                param.propor = 0.13;  % value of the percentile of points to keep
            end

            if isempty(param.cut)
                param.cut = 2;  % size of the cut around the center of the Fourier spectrum
            end

            if isempty(param.stripe)
                param.stripe = false;  % whether to mask for x- y- stripes
            end

            if isempty(param.stripe)
                param.stripewidth = 4;  % default number of fourier space pixels to mask
            end

            if isempty(param.strel)
                param.strel = 4;  % strel value in the soothing filter
            end

            if isempty(param.threshold)
                param.threshold = 0;  % intensity threshold for FT pixels
            end

            if isempty(param.sigma)
                param.sigma = 0.8;  % sigma value in the gaussian soothing kernel
            end

            if isempty(param.register)
                param.register = 0;  %
            end

            if isempty(param.regsize)
                param.regsize = 0;  %
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

        function full_analysis_chunks(obj)
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
            % *OUTPUT*:'Posi': a matrix of size regl x 2 x total number of images that
            %          contains the positions of the center of the subimages at each
            %          times.
            %          'im_regav: a structure containing field 'c' for each slice of
            %          averaged spectrums. In c is found,an intricated structure
            %          containing fields 'S' for the amplitude of the cell deformation,
            %          'angS' for the angle of the cell deformation, 'spect' for the
            %          average spectrum, 'a' for the half major axis of the ellipse in
            %          real space, 'b' the half minor axis of the ellipse in real
            %          space, phi the angle of the ellipse.
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
                
                param=obj.param; %take a copy out of the structs to send to each worker
                Posi = obj.Posi; %intentional, ignore the yellow suggestions
                im_regav=obj.im_regav;
                parfor c = time_lims(1):time_lims(2) %parfor (substitute 'for' for debugging only)
                    %per timepoint
                    Posi_c=Posi(:,:,c);
                    regl = size(Posi_c,1);
                    spectra=spectrum_analysis(param, regl, Posi_c, spectsize, c);        % compute spectrums of subimages
                    results_tmp=def_analysis(param, spectra, regl, c);  % compute cell orientation subimages
                    affich_result_parallel(param, results_tmp, Posi_c, regl, out_folder, c); % create maps of cell orientations
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

            rec = floor((1-obj.param.overlap)*obj.param.tile_size);                      % number of shared pixels between
            % subimages

            if isempty(obj.param.contour)                  % if no mask is required

                obj.regl = zeros(obj.param.tleng,1);   % number of subimages at each time
                for t = 1:obj.param.tleng
                    reg = 0;

                    for x = 1:rec:(obj.param.siz(1)- obj.param.tile_size+1)    % register position
                        for y = 1:rec:(obj.param.siz(2)- obj.param.tile_size+1)
                            reg = reg+1;
                            obj.Posi(reg,:,t) = [x,y];
                        end
                    end
                    obj.regl(t) = reg;               % register number of subimages at each time
                end

            else                                        % when a mask is required

                obj.regl = zeros(obj.param.tleng,1);
                for t = 1:obj.param.tleng
                    reg = 1; % 0 in original
                    b = obj.param.contour_reader.readSpecificImage(obj.param.time_points(t)); % load mask images
                    b = im2double(b);

                    x = 1:rec:(obj.param.siz(1)- obj.param.tile_size+1);
                    y = 1:rec:(obj.param.siz(2)- obj.param.tile_size+1);
                    [X,Y] = meshgrid(x,y);

                    mask = imresize(b, size(X))>=0.5;
                    X = X.*mask;
                    Y = Y.*mask;

                    positions = [X(:) Y(:)];
                    positions(prod(positions,2)==0,:) = [];
                    obj.Posi(reg:size(positions,1),:,t) = positions;
                    reg = size(positions,1);

                    %                     for x = 1:rec:(obj.param.siz(1)- obj.param.tile_size+1)
                    %                         for y = 1:rec:(obj.param.siz(2)- obj.param.tile_size+1)
                    %                              if sum(sum(b(y:y+ obj.param.tile_size-1,x:x+ obj.param.tile_size-1)))<=0.70*obj.param.tile_size^2
                    %                              % Condition on the intensity of the subimages, ie: get rid of borders
                    %                              % if condition not held, don't register position
                    %                              else
                    %                                  reg = reg+1;
                    %                                  obj.results.Posi(reg,:,t) = [x,y];
                    %                              end
                    %                         end
                    %                     end
                    obj.regl(t) = reg;

                end
            end
        end

        function size_analysis(obj)
            % results = def size_analysis(obj.param,results)
            % Function that computes the ellipse for the siz representation
            % on each averaged subimage,  fills the im_regav structure with the fields
            % a, b, phi
            %
            % *OUTPUT*:+ results
            %--------------------------------------------------------------------------

            for c = 1:obj.results.ci                                   % for each averaged time
                display(['Computation of cell size, time: ',num2str(c)]); % for the user to keep track
                for re = 1:obj.results.regl(1+obj.param.timestep*(c-1))     % for each region
                    if sum(sum(to_fill.c(re).spect))==0 % if the spectrum is empty
                        to_fill.c(re).a = 0;            % don't register a, b, phi
                        to_fill.c(re).b = 0;
                        to_fill.c(re).phi = 0;

                    else
                        [to_fill.c(re).a,to_fill.c(re).b,...
                            to_fill.c(re).phi] =...
                            obj.size_def(to_fill.c(re).spect,obj.param);
                        % then compute
                        % the cell size
                        % parameters
                    end

                end
            end
        end

        function affich_size(obj, col)
            % affich_size(obj.param,results,col)
            % A function that plots maps of cell size, registers the figures
            % as .png and .fig
            %
            % + col the colour for the plot of the cell size

            for c = 1:obj.results.ci                            % for each averaged time
                display(['Representation of size map, time: ',num2str(c)]); % for the user to keep track
                Figure = figure(20);                            % Open figure
                imshow(ones(obj.param.siz(2),obj.param.siz(1)));        % Initialize figure with zeros, same size as the original
                hold on                                         % wait

                for win =1:obj.results.regl((c-1)*(obj.param.timestep)+1) % for each subimage
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.tile_size/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.tile_size/2);

                    a = obj.results.im_regav(c).c(win).a;    % extract half major axis
                    b = obj.results.im_regav(c).c(win).b; % extract half minor axis
                    phi = obj.results.im_regav(c).c(win).phi; % extract angle


                    obj.f_ellipseplotax(obj,a,b,xo,yo,phi,col,0.5,Figure.CurrentAxes);
                end
                set(gca,'XTickLabel','','YTickLabel','');
                X3 = [obj.param.siz(1)-250 obj.param.siz(1)-100];
                Y3 = [obj.param.siz(2)-250 obj.param.siz(2)-250];
                plot(X3,Y3,'.b-','LineWidth',1);                                   % scaling bar
                text(obj.param.siz(1)-250,obj.param.siz(2)-300,'0.5','Color','black','FontSize',12);
                img = getframe(gcf);
                imwrite(img.cdata, [obj.param.pathout  ['Size_map' num2str(obj.param.timestep) '_'] num2str(c) '.png']); % save image png
                savefig(Figure,[obj.param.pathout   ['Size_map' num2str(obj.param.timestep) '_'] num2str(c) '.fig'])     % save image fig
                close(Figure);

            end
        end

        function affich_result_original(obj)
            % affich_result(obj.param,results,col)
            % A function that plots maps of cell deformations, registers the figures
            % as .png and .fig
            %          + col the colour for the plot of the cell deformation
            %
            %--------------------------------------------------------------------------
            col = [0.9490    0.7294    0.0471];
            for c = 1:obj.results.ci                            % for each averaged time
                display(['Representation of deformation map, time: ',num2str(c)]); % for the user to keep track

                Figure = figure(20);                            % Open figure
                %                 imshow(ones(obj.param.siz(2),obj.param.siz(1)));        % Initialize figure with zeros, same size as the original
                imshow(obj.param.reader.readSpecificImage(c));
                hold on                                         % wait
                for win =1:obj.results.regl((c-1)*(obj.param.timestep)+1) % for each subimage
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.tile_size/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.tile_size/2);

                    amp = obj.results.im_regav(c).c(win).S;    % extract amplitude
                    ang = obj.results.im_regav(c).c(win).angS; % extract angle of deformation

                    x1 = xo+cos(ang+pi/2)*amp*obj.param.scale; % define the extremities for the plot of deformation
                    y1 = yo+sin(ang+pi/2)*amp*obj.param.scale;
                    x2 = xo-cos(ang+pi/2)*amp*obj.param.scale;
                    y2 = yo-sin(ang+pi/2)*amp*obj.param.scale;

                    X1 = [xo,x1];
                    Y1 = [yo,y1];
                    X2 = [xo,x2];
                    Y2 = [yo,y2];

                    plot(X1,Y1,'Marker','.','LineStyle','-','Color',col,'LineWidth',3); % map cell deformation top half
                    plot(X2,Y2,'Marker','.','LineStyle','-','Color',col,'LineWidth',3); % map cell deformation bottom half

                end

                set(gca,'XTickLabel','','YTickLabel','');
                X3 = [obj.param.siz(1)-250 obj.param.siz(1)-100];
                Y3 = [obj.param.siz(2) obj.param.siz(2)];
                plot(X3,Y3,'.b-','LineWidth',1);                                   % scaling bar
                text(obj.param.siz(1)-250,obj.param.siz(2),'0.5','Color','black','FontSize',12);

                img = getframe(gcf);
                imwrite(img.cdata, [obj.param.pathout  ['Deformation_map' num2str(obj.param.timestep) '_'] num2str(c) '.png']); % save image png
                savefig(Figure,[obj.param.pathout   ['Deformation_map' num2str(obj.param.timestep) '_'] num2str(c) '.fig'])     % save image fig
                close(Figure);

            end
        end

      
        function affich_result_alt_bg(obj, path_alt_bg, out_folder)
            % affich_result(obj.param,results,col)
            % A function that plots maps of cell deformations, registers the figures
            % as .png and .fig
            %          + col the colour for the plot of the cell deformation
            %
            %--------------------------------------------------------------------------
            eR = expReader(path_alt_bg);
            mkdir(out_folder);

            %             col = round(255*[0.2235    0.8588         0]);
            for c = 1:obj.results.ci                            % for each averaged time
                display(['Representation of deformation map, time: ',num2str(obj.param.time_points(c))]); % for the user to keep track
                I = eR.readSpecificImage(obj.param.time_points(c));

                line_array = zeros(obj.results.regl((c-1)*(obj.param.timestep)+1), 4);
                for win = 1:obj.results.regl((c-1)*(obj.param.timestep)+1) % for each subimage
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.tile_size/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.tile_size/2);

                    amp = obj.results.im_regav(c).c(win).S;    % extract amplitude
                    ang = obj.results.im_regav(c).c(win).angS; % extract angle of deformation

                    x1 = xo+cos(ang+pi/2)*amp*obj.param.scale; % define the extremities for the plot of deformation
                    y1 = yo+sin(ang+pi/2)*amp*obj.param.scale;
                    x2 = xo-cos(ang+pi/2)*amp*obj.param.scale;
                    y2 = yo-sin(ang+pi/2)*amp*obj.param.scale;

                    line_array(win,:) = [x1 y1 x2 y2];
                end
                RGB = insertShape(I, 'line', line_array, 'LineWidth', 10, 'Color', 'green');
                imwrite(RGB, [out_folder filesep 'img_' num2str(obj.param.time_points(c), '%04d') '.png']);

            end
        end

        %% Support functions

        function [a,b,phi]= size_def(obj,abs_im_fft_w, param)
            % [a,b,phi]= size_def(abs_im_fft_w,obj.param)
            % Function that computes the halm major and minor axis in the direct space
            % based on a search of local maxima and a fit of ellipse on those local
            % maxima.
            %
            % abs_im_fft_w is a spectrum on which to find the local maximas

            % *OUTPUT*:+ a the half major axis in the direct space
            %          + b the half minor axis in the direct space
            %          + phi the angle of the ellipse
            %--------------------------------------------------------------------------

            [yu,xu] = localMaximum_h(abs_im_fft_w,2,0,40);
            [~,~,ai,bi,phi,~]=ellipsefit(xu,yu);
            a = floor(obj.param.tile_size/2)*1/ai;
            b =floor(obj.param.tile_size/2)*1/bi;
        end

        function [hEllipse,hAxes] = f_ellipseplotax(obj,a,b,x0,y0,phi,col,ldw,axh)
            %  Isabelle Bonnet
            %  Begin: 2009/10/12
            %  Draw an ellipse with given parameters.
            %
            %   Input parameters:
            %       a           Value of the HALF-major axis
            %       b           Value of the HALF-minor axis
            %       x0          Abscissa of the center of the ellipse
            %       y0          Ordinate of the center of the ellipse
            %       phi         Angle between X-axis and the MAJOR axis - units: RADIANS
            %       lineStyle   Definition of the plotted line style
            %
            %   Output:
            %       h    Handle of the ellipse
            %
            %  Usage:
            %       f_ellipseplot(maxi,mini,a,b,phi);
            %       or
            %       hEllipse = f_ellipseplot(5,3,1,-2,pi/4);
            %       set(h,'LineWidth',2);
            %   A slight change added by Melina Durande regarding the handling of the
            %   axes and the gestion of the colors and line width
            %


            if (nargin < 2)|(nargin > 8),
                error('Please see help for INPUT DATA.');
                return;
            end

            angle = [-0.003:0.01:2*pi];

            % parametric equation of the ellipse
            %----------------------------------------
            x = a*cos(angle);
            y = b*sin(angle);


            % Coordinate transform
            %----------------------------------------
            X = cos(phi)*x - sin(phi)*y;
            Y = sin(phi)*x + cos(phi)*y;
            X = X + x0;
            Y = Y + y0;


            % Plot the ellipse
            %----------------------------------------
            hEllipse = plot(X,Y,'Color',col,'parent',axh);
            set(hEllipse,'LineWidth',ldw);
            axis equal;

            %%% Add lines representing the major and minor axes of an
            % ellipse on the current plot.

            % INPUT
            %   axes      [major minor] axes of the ellipse

            % Major Axis coordinates
            x_maj = a*cos(phi);
            y_maj = a*sin(phi);

            % Draw the line
            Major = line(x_maj*[-1 1]+x0,y_maj*[-1 1]+ y0,'Color',col,'parent',axh);
            set(Major,'LineWidth',ldw);

            % Minor Axis coordinates
            x_min = b*cos(phi+pi/2);
            y_min = b*sin(phi+pi/2);

            % Draw the line
            Minor = line(x_min*[-1 1]+ x0,y_min*[-1 1]+ y0,'Color',col,'parent',axh);
            set(Minor,'LineWidth',ldw);
            % Return the line element handles so it can be customized.
            hAxes = [Major Minor];
        end

    end
end