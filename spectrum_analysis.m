function spectra=spectrum_analysis(param, regl, Posi, spectsize, t_index)
% results = spectrum_analysis(obj.param,results)
% expRe
% This function performs Fourier Transform on subimages and averages the
% spectrum on time.
% The im_regav structure is filled with averaged spectrums.
%
% *OUTPUT*: results
%--------------------------------------------------------------------------

if isempty(param.contour)                  % portion of code when a mask is not required
    display(['Computation of FT, time: ',num2str(t_index)]);% for the user to keep track

    a = param.reader.readSpecificImage(param.time_points(t_index));% read the image
    a = im2double(a);
    
    spectra = repmat(zeros(spectsize), [1 1 regl]);
    for win = 1:regl    % for each subimage of the image
        x = Posi(win,1);  % get positionin x and y
        y = Posi(win,2);
        FT = fft_adir(a(y:y+param.tile_size-1,x:x+param.tile_size-1)); % perform FT analysis
        start= floor(3/8*param.tile_size);  % return the middle 1/2 only
        spectra(:,:, win)= FT(start+1:start+spectsize,start+1:start+spectsize);
    end

else                                           % when a mask is required
    display(['Computation of FT, time: ',num2str(t_index)]);
    a = param.reader.readSpecificImage(param.time_points(t_index));
    a = im2double(a);

    b = param.contour_reader.readSpecificImage(param.time_points(t_index));% load mask images
    b = im2double(b);

    for win = 1:results.regl(t_index)
        x = results.Posi(win,1,t_index);
        y = results.Posi(win,2,t_index);
        if x==0 || y==0 || sum(sum(b(y:y+ param.tile_size-1,x:x+ param.tile_size-1)))...
                <=0.70*param.tile_size^2; % Condition on the intensity of the subimages, ie: get rid of borders
            spectrum =...
                results.im_regav(t_index).c(win).spect+zeros(param.tile_size/4); % if condition not held, add zeros
        else
            FT = fft_adir(a(y:y+param.tile_size-1,x:x+param.tile_size-1));

            spectrum=...
                results.im_regav(t_index).c(win).spect+...
                +FT(3/8*param.tile_size+1:5/8*param.tile_size,3/8*param.tile_size+1:5/8*param.tile_size);
        end
    end
end
end