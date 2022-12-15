classdef cellshapefft < handle
%{
    Original code by Durande (2017)
    
    Further modification by Guillermo Serrano Najera (2019)
     > Converted into a class for convenience
     > Plotting functions modified to be faster with very large images
     > Parallelization def_analysis and plotting
     > Array multiplication (.*) for masking (much faster than original)
    
    Code to analyse the cell deformation on large images with Fourier 
    transform. This code generates maps of cell deformation based on the 
    computation of the inertia matrix on Fourier transforms of subimages.
    
    Known issues:
     > When using mask for computations (not only for plotting) the x0, y0
     positions of the subplots can be different between images. At the 
     moment it is better to run the program without mask and apply the mask 
     only for plotting.
    
    To run the code:
        1. Create param struct. Empty values to use default.
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
    param.pas2 = [];
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
        results
    end
    
    methods
        
        function obj = cellshapefft(param)

            %If no contour is needed, enter param.contour = [];
            param.reader = expReader(param.pathin);
            param.contour_reader = expReader(param.contour);
            
            if isempty(param.time_points) % indices of images to analyze
                param.time_points = param.reader.timePoints;
            end
            
            if isempty(param.chunk_size) % indices of images to analyze
                param.chunk_size = 16;
            end
            
            if isempty(param.tleng) % number of images to analyze
                param.tleng = length(param.time_points);
            end
            
            if isempty(param.timestep)
                param.timestep = 1; % time stack to perform averages
            end
            
            if isempty(param.pas2)
                 param.pas2 = 128; % size of subimages
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
            
            if isempty(param.strel)
                 param.strel = 4;  % strel value in the soothing filter
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
                 param.ffmpeg_path = "ffmpeg";  % 
            end
            
            if isempty(param.workers)
                 param.workers = 4;  % cores for parallel computation
            end
            
            obj.param = param;
            info = imfinfo(obj.param.reader.fileName); % get information about the images 
            obj.param.siz=[info.Width info.Height];        % extract size of the images
            obj.results = struct('regl',[],'Posi',[],'im_regav',[],'ci',[]); % result structure
            mkdir(obj.param.pathout);

            %save a copy of input file
            try
            if ispc
                system(['copy /b callcellshapefft.m ' obj.param.pathout filesep]);
            else %unix
                system(['cp ./callcellshapefft.m' obj.param.pathout filesep]);
            end
            catch
                %if pathing problems, skip
                disp(["Problems encountered while attempting to copy callcellshapefft.m to" ...
                    obj.param.pathout ", input parramters were not saved"])
            end
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
        %        'pas2': the size of the subimages
        %        'sigma': sigma value in the gaussian soothing kernel
        %        'cut': size of the cut around the center of the Fourier spectrum
        %        'proportion': value of the percentile of points to keep
        %        'strel': strel value in the soothing filter
        %
        % *OUTPUT*: + results, a structure containing fields for
        %          'regl': a list of the number of subimages in each image
        %          'Posi': a matrix of size regl x 2 x total number of images that
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

            obj.Position();                 % register positions of each subimages
            obj.spectrum_analysis; % compute averages spectrums of subimages
            obj.save_object()
            obj.def_analysis_parallel();     % compute cell orientation subimages
            obj.save_object()
            if obj.param.regsize == 1 %GUILLERMO> WHAT IS THIS?
               obj.size_analysis();
               obj.affich_size('red');
            else
            end
            obj.affich_result_parallel();        % create maps of cell orientations

        end
        
        function full_analysis_chunks(obj)
            pool = gcp('nocreate');
            if isempty(pool)
                parpool(obj.param.workers);
            else
                if pool.NumWorkers ~= obj.param.workers
                    delete(gcp('nocreate'))
                    parpool(obj.param.workers); 
                end
            end

            time_vector = obj.param.time_points;
            for t_ind = time_vector(1:obj.param.chunk_size:end)
                time_lims = [t_ind min(t_ind+obj.param.chunk_size-1, time_vector(end))];
                obj.param.time_points = time_lims(1):time_lims(2);
                obj.param.tleng = length(obj.param.time_points);
                
                obj.Position();               % register positions of each subimages
                obj.spectrum_analysis;        % compute averages spectrums of subimages
                obj.def_analysis_parallel();  % compute cell orientation subimages
                obj.save_object(time_lims)    % save results
                obj.affich_result_parallel(); % create maps of cell orientations
                
                obj.results = struct('regl',[],'Posi',[],'im_regav',[],'ci',[]); % result structure
            end
                rR = expReader([obj.param.pathout filesep 'img']);
                rR.ffmpeg_path = obj.param.ffmpeg_path;
                rR.ffmpeg;
        end
        
        function save_object(obj, time_lims)
            
            save([obj.param.pathout filesep 'cellshapefft_results_' num2str(time_lims(1),'%04d') '_' num2str(time_lims(2),'%04d') '.mat'], 'obj', '-v7.3');% save results
        end
        
        function Position(obj)  
            % results = Position(obj.param)
            %
            % This function registers the positions of each subimage at each
            % time, and the number of subimage at each time. The regl and Posi fields 
            % of the Result structure are filled.
            disp('Registering...');
                        
            rec = floor(1/2*obj.param.pas2);                      % number of shared pixels between
                                                       % subimages

            if isempty(obj.param.contour)                  % if no mask is required

                obj.results.regl = zeros(obj.param.tleng,1);   % number of subimages at each time
                for t = 1:obj.param.tleng
                        reg = 0;
                        % GUILLERMO. WHY LOADING IMAGE?
%                         a = obj.param.reader.readSpecificImage(t); % read image 
%                         a = im2double(a);

                    for x = 1:rec:(obj.param.siz(1)- obj.param.pas2+1)    % register position
                        for y = 1:rec:(obj.param.siz(2)- obj.param.pas2+1)
                                 reg = reg+1;
                                 obj.results.Posi(reg,:,t) = [x,y];        
                        end
                    end
                   obj.results.regl(t) = reg;               % register number of subimages at each time
                end

            else                                        % when a mask is required

                obj.results.regl = zeros(obj.param.tleng,1);
                for t = 1:obj.param.tleng
                        reg = 1; % 0 in original
                        b = obj.param.contour_reader.readSpecificImage(obj.param.time_points(t)); % load mask images
                        b = im2double(b);

                        x = 1:rec:(obj.param.siz(1)- obj.param.pas2+1);
                        y = 1:rec:(obj.param.siz(2)- obj.param.pas2+1);
                        [X,Y] = meshgrid(x,y);
                        
                        mask = imresize(b, size(X))>=0.5;
                        X = X.*mask;
                        Y = Y.*mask;
                        
                        positions = [X(:) Y(:)];
                        positions(prod(positions,2)==0,:) = [];
                        obj.results.Posi(reg:size(positions,1),:,t) = positions; 
                        reg = size(positions,1);

%                     for x = 1:rec:(obj.param.siz(1)- obj.param.pas2+1)
%                         for y = 1:rec:(obj.param.siz(2)- obj.param.pas2+1)
%                              if sum(sum(b(y:y+ obj.param.pas2-1,x:x+ obj.param.pas2-1)))<=0.70*obj.param.pas2^2  
%                              % Condition on the intensity of the subimages, ie: get rid of borders
%                              % if condition not held, don't register position
%                              else
%                                  reg = reg+1;
%                                  obj.results.Posi(reg,:,t) = [x,y];        
%                              end
%                         end
%                     end
                   obj.results.regl(t) = reg;
                   
                end
            end
        end
        
        function spectrum_analysis(obj)
        % results = spectrum_analysis(obj.param,results)
        %expRe
        % This function performs Fourier Transform on subimages and averages the 
        % spectrum on time.
        % The im_regav structure is filled with averaged spectrums.
        %
        % *OUTPUT*: results
        %--------------------------------------------------------------------------


            obj.results.im_regav = struct(repmat(struct('c',struct(repmat(struct('spect',zeros(ceil(obj.param.pas2/4)),'S',[]...
                ,'angS',[],'a',[],'b',[],'phi',[]),max(obj.results.regl),1))),obj.param.tleng/obj.param.timestep,1));
            % Initialization of the structure im_regav, containing spectrums
            % cell deformations and sizes for each averaged subimage
            %TODO why divide tile size by 4? seems to only use middle half
            %of Fourier transforms

            obj.results.ci =0;                             % ci will be the final number of slices
            if isempty(obj.param.contour)                  % portion of code when a mask is not required
                for to =1:obj.param.timestep:obj.param.tleng-obj.param.timestep+1  % index of first images for each slice
                    display(['Computation of FT, time: ',num2str(obj.param.time_points(to))]);% for the user to keep track
                    obj.results.ci = obj.results.ci+1;         % increment the counter of ci 
                    for t = to:to+obj.param.timestep-1     % for the images index contained between to and
                                                       % to+thickness of the average
                        a = obj.param.reader.readSpecificImage(obj.param.time_points(t));% read the image
                        a = im2double(a);

                        for win = 1:obj.results.regl(t)    % for each subimage of the image
                           x = obj.results.Posi(win,1,t);  % get positionin x and y
                           y = obj.results.Posi(win,2,t);
                           FT = obj.fft_adir(a(y:y+obj.param.pas2-1,x:x+obj.param.pas2-1));
                           start= floor(3/8*obj.param.pas2);
                           len = size(obj.results.im_regav(obj.results.ci).c(win).spect);
                           obj.results.im_regav(obj.results.ci).c(win).spect=obj.results.im_regav(obj.results.ci).c(win).spect...
                               +FT(start+1:start+len(1),start+1:start+len(1)); % perform FT analysis 
                        end
                    end
                    obj.results.im_regav(obj.results.ci).c(win).spect =...
                        obj.results.im_regav(obj.results.ci).c(win).spect/obj.param.timestep; % divide intensity by the number of spectrums averaged
                end

            else                                           % when a mask is required

                for to = 1:obj.param.timestep:obj.param.tleng-obj.param.timestep+1
                    Para = obj.param.timestep;                 % counter for the number of images to average
                    display(['Computation of FT, time: ',num2str(obj.param.time_points(to))]);
                    obj.results.ci = obj.results.ci+1;
                    for t = to:to+obj.param.timestep-1
                        a = obj.param.reader.readSpecificImage(obj.param.time_points(t));
                        a = im2double(a);

                        b = obj.param.contour_reader.readSpecificImage(obj.param.time_points(t));% load mask images
                        b = im2double(b);

                        for win = 1:obj.results.regl(t)
                           x = obj.results.Posi(win,1,t);
                           y = obj.results.Posi(win,2,t);
                            if x==0 || y==0 || sum(sum(b(y:y+ obj.param.pas2-1,x:x+ obj.param.pas2-1)))...
                                    <=0.70*obj.param.pas2^2; % Condition on the intensity of the subimages, ie: get rid of borders
                                obj.results.im_regav(obj.results.ci).c(win).spect =...
                                    obj.results.im_regav(obj.results.ci).c(win).spect+zeros(obj.param.pas2/4); % if condition not held, add zeros
                                Para = obj.param.timestep-1;                                                                                 % decrement the countour on which to average afterwards
                            else
                                   FT = obj.fft_adir(a(y:y+obj.param.pas2-1,x:x+obj.param.pas2-1));

                                obj.results.im_regav(obj.results.ci).c(win).spect=...
                                    obj.results.im_regav(obj.results.ci).c(win).spect+...
                                    +FT(3/8*obj.param.pas2+1:5/8*obj.param.pas2,3/8*obj.param.pas2+1:5/8*obj.param.pas2);    
                            end

                        end
                        
                    end
                    obj.results.im_regav(obj.results.ci).c(win).spect = obj.results.im_regav(obj.results.ci).c(win).spect/Para;
                end
            end
        end
        
        function def_analysis(obj)
            % results = def analysis(obj.param,results)
            % Function that computes the inertia matrix and cell deformation on 
            % each averaged subimage, fills the im_regav structure with the fields S
            % and angS.
            %
            % *OUTPUT*: results
            %--------------------------------------------------------------------------

             for c = 1:obj.results.ci                                   % for each averaged time
                display(['Computation of cell deformation, time: ',num2str(c)]); % for the user to keep track
                for re = 1:obj.results.regl(1+obj.param.timestep*(c-1))     % for each region
                    if sum(sum(obj.results.im_regav(c).c(re).spect))==0 % if the spectrum is empty
                        obj.results.im_regav(c).c(re).angS = 0;         % don't register cell deformation
                        obj.results.im_regav(c).c(re).S = 0;

                    else
                        M = inertia_matp_sigma(obj.results.im_regav(c).c(re).spect,obj.param);      
                                                                    % else compute
                                                                    % the inertia matrix
                        [obj.results.im_regav(c).c(re).S,obj.results.im_regav(c).c(re).angS] =S_angS(M);
                                                                    % then compute the cell
                                                                    % deformation amplitude and angle
                    end
                    if obj.param.register == 0 && obj.param.regsize == 0
                       obj.results.im_regav(c).c(re).spect = [];       % delete the spectrum to keep space
                    else
                    end
                end
             end
        end
        
        function def_analysis_parallel(obj)
            % results = def analysis(obj.param,results)
            % Function that computes the inertia matrix and cell deformation on 
            % each averaged subimage, fills the im_regav structure with the fields S
            % and angS.
            %
            % *OUTPUT*: results
            %--------------------------------------------------------------------------

            im_regav = cell(obj.param.tleng,1);
            
             parfor c = 1:obj.results.ci                                   % for each averaged time
                display(['Computation of cell deformation, time: ',num2str(obj.param.time_points(c))]); % for the user to keep track
                to_fill = obj.results.im_regav(c);
                for re = 1:obj.results.regl(1+obj.param.timestep*(c-1))     % for each region
                    if sum(sum(obj.results.im_regav(c).c(re).spect))==0 % if the spectrum is empty
                        to_fill.c(re).angS = 0;         % don't register cell deformation
                        to_fill.c(re).S = 0;

                    else  % else compute the inertia matrix
                        M = obj.inertia_matp_sigma(to_fill.c(re).spect,obj.param);      
                        % then compute the cell deformation amplitude and angle
                        [to_fill.c(re).S,to_fill.c(re).angS] = obj.S_angS(M);                                                                    
                    end
                    
                    if obj.param.register == 0 && obj.param.regsize == 0
                       to_fill.c(re).spect = [];       % delete the spectrum to keep space
                    else
                    end
                end
                im_regav{c} = to_fill;
             end
             
            % return to the original organization
            for c = 1:obj.results.ci
                obj.results.im_regav(c) = im_regav{c};
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
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2);

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
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2);

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
        
        function affich_result(obj)
            % affich_result(obj.param,results,col)
            % A function that plots maps of cell deformations, registers the figures
            % as .png and .fig
            %          + col the colour for the plot of the cell deformation 
            %
            %--------------------------------------------------------------------------
            out_folder = [obj.param.pathout filesep 'img'];
            mkdir(out_folder);
%             col = round(255*[0.2235    0.8588         0]);           
            for c = 1:obj.results.ci                            % for each averaged time
                
                display(['Representation of deformation map, time: ',num2str(obj.param.time_points(c))]); % for the user to keep track
                I = imadjust(imcomplement(obj.param.reader.readSpecificImage(obj.param.time_points(c))));
                line_array = zeros(obj.results.regl((c-1)*(obj.param.timestep)+1), 4);
                for win = 1:obj.results.regl((c-1)*(obj.param.timestep)+1) % for each subimage
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2);

                    amp = 2*obj.results.im_regav(c).c(win).S;    % extract amplitude 
                    ang = obj.results.im_regav(c).c(win).angS; % extract angle of deformation

                    x1 = xo+cos(ang+pi/2)*amp*obj.param.scale; % define the extremities for the plot of deformation
                    y1 = yo+sin(ang+pi/2)*amp*obj.param.scale;
                    x2 = xo-cos(ang+pi/2)*amp*obj.param.scale;
                    y2 = yo-sin(ang+pi/2)*amp*obj.param.scale;
                    
                    line_array(win,:) = [x1 y1 x2 y2];
                end
                        
                RGB = insertShape(I, 'line', line_array, 'LineWidth', 30, 'Color', 255*[ 0.9490    0.7294    0.0471]);
                
                X3 = [obj.param.siz(1)-1000 obj.param.siz(1)-700];
                Y3 = [0.9*obj.param.siz(2) 0.9*obj.param.siz(2)];
                line_array = [X3(1) Y3(1) X3(2) Y3(2)]; %scale bar
                RGB = insertShape(RGB, 'line', line_array, 'LineWidth', 30, 'Color', 'black');
                RGB = insertText(RGB, [X3(1)+diff(X3)/2 Y3(1)+diff(Y3)/2; X3(1)+diff(X3)/2 Y3(1)+diff(Y3)/2+150], {'0.5', '0.195 mm'}, 'AnchorPoint', 'CenterTop', 'FontSize', 120, 'BoxOpacity', 0);
                        
                imwrite(RGB, [out_folder filesep 'img_' num2str(obj.param.time_points(c), '%04d') '.png']);

            end
        end
        
        function affich_result_parallel(obj)
            % affich_result(obj.param,results,col)
            % A function that plots maps of cell deformations, registers the figures
            % as .png and .fig
            %          + col the colour for the plot of the cell deformation 
            %
            %--------------------------------------------------------------------------
            out_folder = [obj.param.pathout filesep 'img'];
            
            mkdir(out_folder);
            
%             col = round(255*[0.2235    0.8588         0]);           
            parfor c = 1:obj.results.ci                            % for each averaged time
                display(['Representation of deformation map, time: ',num2str(obj.param.time_points(c))]); % for the user to keep track
%                 I = imadjust(imcomplement(obj.param.reader.readSpecificImage(obj.param.time_points(c))));
                I = imadjust((obj.param.reader.readSpecificImage(obj.param.time_points(c))));
                line_array = zeros(obj.results.regl((c-1)*(obj.param.timestep)+1), 4);
                for win = 1:obj.results.regl((c-1)*(obj.param.timestep)+1) % for each subimage
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2);

                    amp = obj.results.im_regav(c).c(win).S;    % extract amplitude 
                    ang = obj.results.im_regav(c).c(win).angS; % extract angle of deformation

                    x1 = xo+cos(ang+pi/2)*amp*obj.param.scale; % define the extremities for the plot of deformation
                    y1 = yo+sin(ang+pi/2)*amp*obj.param.scale;
                    x2 = xo-cos(ang+pi/2)*amp*obj.param.scale;
                    y2 = yo-sin(ang+pi/2)*amp*obj.param.scale;
                    
                    line_array(win,:) = [x1 y1 x2 y2];
                end
                RGB = insertShape(I, 'line', line_array, 'LineWidth', 5, 'Color', 'yellow');
                imwrite(RGB, [out_folder filesep 'img_' num2str(obj.param.time_points(c), '%04d') '.png']);

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
                    xo =  obj.results.Posi(win,1,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2); % find positions of the center of the subimages
                    yo =  obj.results.Posi(win,2,(c-1)*(obj.param.timestep)+1)+floor(obj.param.pas2/2);

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
        function [abs_im_fft] = fft_adir(obj,im)
            % abs_im_fft = fft_adir(im)
            % function that computes the FT of an image

            %--------------------------------------------------------------------------
            %
            % *INPUT*: + im the image on which to perform the fourier transform
            %
            % *OUTPUT*: + abs_im_fft the spectrum normalized
            %
            %--------------------------------------------------------------------------

                if(isnan(im)==1)                % if image has NaN
                    abs_im_fft = NaN;
                else
                    [p,~] = obj.perdecomp(im);      % reduce image size effects by periodizing borders
                    p=imgaussfilt(p,obj.param.stripe_sigma);    
                    im_fft = fftn(p);               % 2D fast fourier transform
                    im_fft_shift = fftshift(im_fft);% shifts the FT for representation purpose
                    abs_im_fft = abs(im_fft_shift); % takes the module
                    [~,ind]=max(abs_im_fft(:));     % finds index of the max
                    sumabs = sum(sum(abs_im_fft));  % computes total value 
                    abs_im_fft(ind) = 0;            % puts center to zero

                    abs_im_fft = abs_im_fft/sumabs; % normalize the rest of the image
                end
        end
        
        function [p,s] = perdecomp(obj,u)
            %       Periodic plus Smooth Image Decomposition
            %
            %               author: Lionel Moisan
            %
            %   This program is freely available on the web page
            %
            %   http://www.mi.parisdescartes.fr/~moisan/p+s
            %
            %   I hope that you will find it useful.
            %   If you use it for a publication, please mention 
            %   this web page and the paper it makes reference to.
            %   If you modify this program, please indicate in the
            %   code that you did so and leave this message.
            %   You can also report bugs or suggestions to 
            %   lionel.moisan [AT] parisdescartes.fr
            %
            % This function computes the periodic (p) and smooth (s) components
            % of an image (2D array) u
            %
            % usage:    p = perdecomp(u)    or    [p,s] = perdecomp(u)
            %
            % note: this function also works for 1D signals (line or column vectors)
            %
            % v1.0 (01/2014): initial Matlab version from perdecomp.sci v1.2
            [ny,nx] = size(u); 
            u = double(u);
            X = 1:nx; Y = 1:ny;
            v = zeros(ny,nx);
            v(1,X)  = u(1,X)-u(ny,X);
            v(ny,X) = -v(1,X);
            v(Y,1 ) = v(Y,1 )+u(Y,1)-u(Y,nx);
            v(Y,nx) = v(Y,nx)-u(Y,1)+u(Y,nx);
            fx = repmat(cos(2.*pi*(X -1)/nx),ny,1);
            fy = repmat(cos(2.*pi*(Y'-1)/ny),1,nx);
            fx(1,1)=0.;   % avoid division by 0 in the line below
            s = real(ifft2(fft2(v)*0.5./(2.-fx-fy)));
            p = u-s;

        end
        
        function [S,angS] = S_angS(obj, M)
            % [S,angS] = S_angS(M)
            % This function defines the angle and amplitude of cell deformation
            % from inertia matrix

            %--------------------------------------------------------------------------
            %
            % *INPUT*: M inertia matrix
            %
            % *OUTPUT*: + S the amplitude of the cell deformation defined by the 
            %           inertia matrix
            %           + angS angle of the cell deformation defined by the inertia
            %           matrix
            %
            %--------------------------------------------------------------------------

            [~,E] = eig(M);                % extract the eigenvalues of the inertia matrix
            L1 = 2*sqrt(E(2,2)/2);         % definition of the minor and major axis of the inertia ellipse     
            L2 = 2*sqrt(E(1,1)/2);
            lgE = -1/2*(log(L1/L2));       % amplitude of the cell deformation
            Dev = M-1/2*trace(M)*eye(2);   % deviation matrix
            [V,~] = eig(Dev);              % eigenvectors of deviation matrix to extract the angle
            ang=atan(V(1,2)/V(2,2));       % definition of angle
            angS = ang;                    % register angle
            S = lgE;                       % register amplitude
        end
        
        function [M] = inertia_matp_sigma(obj, al, param)
            % M = inertia_matp_sigma(al,obj.param)
            % function that computes the inertia matrix based on a spectrum averaged
            % smoothed with a gaussian filter and with a percentile of points to keep.
            %
            % *OUTPUT*: + M the inertia matrix
            %--------------------------------------------------------------------------
            
            filtre = obj.gaussianFilter(ceil(2*obj.param.sigma),obj.param.sigma);
            th_spec = filter2(filtre,al);      % th_spec gaussian filter of the spectrum

            ic = obj.param.pas2*1/8+1;
            jc = obj.param.pas2*1/8+1;              % get rid of the pixels in the center in a radius
                                                % of obj.param.cut
                                                % because the proportion was defined this way
            [I,J] = ndgrid(1:obj.param.pas2/4,1:obj.param.pas2/4);
            M = double((I-ic).^2+(J-jc).^2<=(obj.param.cut)^2); 
            th_spec(M>0) = 0;


            reshaped = reshape(th_spec,[obj.param.pas2*obj.param.pas2*1/16 ,1]);  % reshape the spectrum  
            [~,index]= sort(reshaped,'descend');         % sort in intensity and keep indexes
            nbel = obj.param.propor * obj.param.pas2*obj.param.pas2*1/16;         % nbel number of points determined by obj.param.propor

            if nbel > obj.param.pas2*obj.param.pas2/16             % if nbel is to large, take the whole picture
                nbel = obj.param.pas2*obj.param.pas2/16;
            else
            end

            th_spec(index(1:floor(nbel))) = 1;           % put to one the points kept by proportion
            th_spec(index(floor(nbel)+1:end)) = 0;       % put to zero the others

            th_spec(M>0) = 1;                            % get back the pixels in the center in a 
                                                         % radius of obj.param.cut 

            B = false(obj.param.pas2/4);                       % create a logical matrix of falses of size pas2         
            B(th_spec>0)=1;                              % put rights where th_spec = 1
            se=strel('disk',obj.param.strel);                % smooth with a strel
            th_spec=imclose(B,se);                       % close image

            norm = sum(sum((th_spec).^2));               % computes norm 
            [Yg,Xg]=meshgrid(1:(obj.param.pas2/4),1:(obj.param.pas2/4)); % create meshgrid of correct size
            xc = obj.param.pas2/8;                               % position of the center 
            yc = obj.param.pas2/8;
            Xo = Xg-xc;                           
            Yo = Yg-yc;
            Y = (Yg-yc).^2; % Lignes
            X = (Xg-xc).^2; % Colonnes
            YY = sum(sum((th_spec).^2.*Y))/norm;         % compute the 3 different values in the inertia matrix
            XX = sum(sum((th_spec).^2.*X))/norm;
            XY = sum(sum((th_spec).^2.*Xo.*Yo))/norm;
            M = [XX,XY;XY,YY];                           % keep matrix of inertia
        end
        
        function filter = gaussianFilter(obj, r,s)
            % Make a gaussian filter of radius 'r' and standard deviation 's'
            %   filter = gaussianFilter(r,s)

            [x,y]  = meshgrid(-r:r,-r:r);
            filter = exp(-(x.^2 + y.^2)/(2*s.^2));
            filter = filter./sum(filter(:));
        end
        
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
            a = floor(obj.param.pas2/2)*1/ai;
            b =floor(obj.param.pas2/2)*1/bi;
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