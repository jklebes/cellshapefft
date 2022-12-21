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
        FT = fft_adir(a(y:y+param.tile_size-1,x:x+param.tile_size-1), param.stripe_sigma); % perform FT analysis
        start= floor(3/8*param.tile_size);  % return the middle 1/2 only
        spectra(:,:, win)= FT(start+1:start+spectsize,start+1:start+spectsize);
    end

else                                           % when a mask is required
    display(['Computation of FT, time: ',num2str(t_index)]);
    a = param.reader.readSpecificImage(param.time_points(t_index));
    a = im2double(a);

    b = param.contour_reader.readSpecificImage(param.time_points(t_index));% load mask images
    b = im2double(b);

    spectra = repmat(zeros(spectsize), [1 1 regl]);
    for win = 1:regl
        x = Posi(win,1);
        y = Posi(win,2);
        if x==0 || y==0 || sum(sum(b(y:y+ param.tile_size-1,x:x+ param.tile_size-1)))...
                <=0.70*param.tile_size^2 % Condition on the intensity of the subimages, ie: get rid of borders
            % if condition not held, spectrum remains a zeros matrix
        else
            FT = fft_adir(a(y:y+param.tile_size-1,x:x+param.tile_size-1), param.stripe_sigma);
            start= floor(3/8*param.tile_size);
            spectra(:,:, win)=FT(start+1:start+spectsize,start+1:start+spectsize);
        end
    end
end
end

function [abs_im_fft] = fft_adir(im, sigma)
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
                if ~isempty(sigma)
                    im = imgaussfilt(im, sigma);
                end
                [p,~] = perdecomp(im);      % reduce image size effects by periodizing borders
                im_fft = fftn(p);               % 2D fast fourier transform
                im_fft_shift = fftshift(im_fft);% shifts the FT for representation purpose
                abs_im_fft = abs(im_fft_shift); % takes the module
                [~,ind]=max(abs_im_fft(:));     % finds index of the max
                sumabs = sum(sum(abs_im_fft));  % computes total value
                abs_im_fft(ind) = 0;            % puts center to zero

                abs_im_fft = abs_im_fft/sumabs; % normalize the rest of the image
            end
        end