function Results = Position(Param)  
% Results = Position(Param)
%
% This function registers the positions of each subimage at each
% time, and the number of subimage at each time. The regl and Posi fields 
% of the Result structure are filled.
%
%--------------------------------------------------------------------------
%
% *INPUT*: Param for the analysis in Fourier transform, Param is a structure
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
%        'register': is 1 whe the spectrums are registered
%        'regsize': 1 if you want to register the size analysis
%        'nbpoints': the number of points to keep in the fit for the size
%            analysis
%        'sizeview': 1 if you want to visualise the analysis in size
%        'rec': the value of the overlap

% *OUTPUT*: Results, a structure containing fields for
%          'regl': a list of the number of subimages in each image
%          'Posi': a matrix of size regl x 2 x total number of images that
%          contains the positions of the center of the subimages at each
%          times.
%          'im_regav: a structure containing field 'c' for each slice of
%          averaged spectrums. In c is found,an intricated structure 
%          containing fields 'S' for the amplitude of the cell deformation,
%          'angS' for the angle of the cell deformation, 'spect' for the 
%          average spectrum,'a' for the half major axis of the ellipse in
%          real space, 'b' the half minor axis of the ellipse in real
%          space, phi the angle of the ellipse.
%
%--------------------------------------------------------------------------



    rec = (1-Param.rec)*Param.pas2;                      % number of shared pixels between
                                               % subimages
        
    if isempty(Param.contour)                  % if no mask is required
        Anames = ([Param.name]); % path of the images
        Anamesi = dir(Anames);                 % get all the images names
        Anamesi([Anamesi.isdir]) = [];         % get rif of hidden files
        A = strfind({Anamesi.name},'.DS_Store');
        Anamesi = Anamesi(cellfun('isempty',A));
    A = strfind({Anamesi.name},'.mat');
    Anamesi = Anamesi(cellfun('isempty',A));
   
        
        Results.regl = zeros(Param.tleng-Param.tsart,1);   % number of subimages at each time
        for t = Param.tsart:Param.tleng
                reg = 0;

            if numel(extractfield(Anamesi,'name'))==1
                ibl = extractfield(Anamesi,'name');
                
                if numel(strfind([Anames,  ibl{1}],'.png'))>0
                    a = imread([Anames,  ibl{1}]);
                    a = im2double(a);
                   Param.siz= size(a);

                else
                    a = imread([Anames,  ibl{1}],t);
                    a = im2double(a);
                    Param.siz= size(a);
                end

            else
                 a = imread([Anames,  Anamesi(t).name]);
                 a = im2double(a);         
                 Param.siz= size(a);

            end
            
            for x = 1:rec:(Param.siz(2)- Param.pas2+1)    % register position
                for y = 1:rec:(Param.siz(1)-  Param.pas2+1)
                         reg = reg+1;
                         Results.Posi(reg,:,t-Param.tsart+1) = [floor(x),floor(y)];        
                end
            end
           Results.regl(t-Param.tsart+1) = reg;               % register number of subimages at each time
        end
        
    else                                        % when a mask is required
        Bnames = ([Param.pathin, Param.contour]);  % path of the mask images
        Bnamesi = dir(Bnames);                     % get mask images names
        Bnamesi([Bnamesi.isdir]) = [];             % get rif of hidden files
        B = strfind({Bnamesi.name},'.DS_Store');
        Bnamesi = Bnamesi(cellfun('isempty',B));


        Results.regl = zeros(Param.tleng,1);
        for t = Param.tsart:Param.tleng
                reg = 0;           
            if numel(extractfield(Bnamesi,'name'))==1
                ibl = extractfield(Bnamesi,'name');
                 b = imread([Bnames,  ibl{1}],t);
                 b= im2double(b);
            else
                 b = imread([Bnames, Bnamesi(t).name]);
                 b = im2double(b);         
            end

            for x = 1:rec:(Param.siz(1)- Param.pas2+1)
                for y = 1:rec:(Param.siz(2)- Param.pas2+1)
                     if sum(sum(b(y:y+ Param.pas2-1,x:x+ Param.pas2-1)))<=0.70*Param.pas2^2  
                     % Condition on the intensity of the subimages, ie: get rid of borders
                     % if condition not held, don't register position
                     else
                         reg = reg+1;
                         Results.Posi(reg,:,t-Param.tsart+1) =[floor(x),floor(y)];       
                     end
                end
            end
           Results.regl(t-Param.tsart+1) = reg;
        end
    end
end