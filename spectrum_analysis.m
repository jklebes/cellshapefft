function spectra=spectrum_analysis(param, data, nTiles, tileCoords, spectsize, t_index)
% 
% This function performs Fourier Transform on subimages and averages the
% spectrum on time.
% The spectra structure is filled with middle area of Fourier transofrmed
% tiles
%
% *OUTPUT*: spectra: array of dimensions [spectsize spectsize nTiles] 
%--------------------------------------------------------------------------
if isempty(param.contour)                  % portion of code when a mask is not required
    display(['Computation of FT, time: ',num2str(t_index)]);% for the user to keep track
    if isempty(data)
    a = param.reader.readSpecificImage(param.time_points(t_index));% read the image
    a = im2double(a);
    else
        a=data;
    end
    spectra = repmat(zeros(spectsize), [1 1 nTiles]);

    %quality=zeros([1 nTiles]);
    for win = 1:nTiles    % for each subimage of the image
        x = tileCoords(win,1);  % get positionin x and y
        y = tileCoords(win,2);
        FT = fftTile(a(y:y+param.tileSize-1,x:x+param.tileSize-1), param.sigma); % perform FT analysis
        %quality(win)=signalToNoise(a(y:y+param.tileSize-1,x:x+param.tileSize-1), 4);
        start= floor(3/8*param.tileSize);  % return the middle 1/2 only
        spectra(:,:, win)= FT(start+1:start+spectsize,start+1:start+spectsize);
    end

else                                           % when a mask is required
    display(['Computation of FT, time: ',num2str(t_index)]);
    a = param.reader.readSpecificImage(param.time_points(t_index));
    a = im2double(a);

    b = param.contour_reader.readSpecificImage(param.time_points(t_index));% load mask images
    b = im2double(b);

    spectra = repmat(zeros(spectsize), [1 1 nTiles]);
    quality=zeros([nTiles 1]);
    for win = 1:nTiles
        x = tileCoords(win,1);
        y = tileCoords(win,2);
        if x==0 || y==0 || sum(sum(b(y:y+ param.tileSize-1,x:x+ param.tileSize-1)))...
                <=0.70*param.tileSize^2 % Condition on the intensity of the subimages, ie: get rid of borders
            % if condition not held, spectrum remains a zeros matrix
        else
            FT = fftTile(a(y:y+param.tileSize-1,x:x+param.tileSize-1), param.stripe_sigma);
            quality(:,win)=mean_intensity(a(y:y+param.tileSize-1,x:x+param.tileSize-1));
            start= floor(3/8*param.tileSize);
            spectra(:,:, win)=FT(start+1:start+spectsize,start+1:start+spectsize);
        end
    end
end
end

function [abs_im_fft] = fftTile(im, sigma)
            % function that computes the FT of an image

            %--------------------------------------------------------------------------
            %
            % *INPUT*: + im the image on which to perform the fourier transform
            %
            % *OUTPUT*: + abs_im_fft the spectrum normalized
            %
            %--------------------------------------------------------------------------
            if isnan(im)              % if image has NaN
                abs_im_fft = NaN;
            else
                [p,~] = periodicDecomposition(im); % reduce image size effects by periodizing borders
                im_fft = fftn(p);               % 2D fast fourier transform
                im_fft_shift = fftshift(im_fft);% shifts the FT for representation purpose
                abs_im_fft = abs(im_fft_shift); % takes the modulus
                [~,ind]=max(abs_im_fft(:));     % finds index of the max = center
                %sumabs = sum(sum(abs_im_fft));  % computes total value
                if ~isempty(sigma)
                    abs_im_fft=imgaussfilt(abs_im_fft, sigma); %blurring of Fourier space image
                end
                abs_im_fft(ind) = 0;            % puts center to zero
                siz=size(abs_im_fft,1);

                %% the profile-subtracting correction!
                %need a tile-size-independent way to figure out where to
                %cut for top and bottom regions that are definitely just
                %noise and no signal
                top_bottom=[abs_im_fft(1:round(siz/8)+1,:);abs_im_fft(end-round(siz/8):end,:)];
                col_means = mean(top_bottom, 1);
                abs_im_fft =abs_im_fft-col_means;

                %abs_im_fft = abs_im_fft/sumabs; %normalize by tile area ?
            end
end
