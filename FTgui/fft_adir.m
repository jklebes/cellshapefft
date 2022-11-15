function [abs_im_fft] = fft_adir(im)
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

    [p,~] = perdecomp(im);          % reduce image size effects by periodizing borders
   
   im_fft = fftn(p);               % 2D fast fourier transform
    im_fft_shift = fftshift(im_fft);% shifts the FT for representation purpose
    abs_im_fft = abs(im_fft_shift).^2; % takes the module
    [~,ind]=max(abs_im_fft(:));     % finds index of the max
    sumabs = sum(sum(abs_im_fft));  % computes total value 
    abs_im_fft(ind) = 0;            % puts center to zero

    abs_im_fft = abs_im_fft/size(im,1)^2; % normalize the rest of the image
    end
end