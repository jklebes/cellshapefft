function [quality,Ms] = inertia_matp_sigma(param, spectra, regl)
% M = inertia_matp_sigma(al,obj.param)
% function that computes the inertia matrix based on a spectrum averaged
% smoothed with a gaussian filter and with a percentile of points to keep.
%
% *OUTPUT*: + M the inertia matrix
%--------------------------------------------------------------------------
Ms = zeros(2,2,regl);
quality=zeros([1 regl]);
for re=1:regl
al=spectra(:,:,re);
th_spec = al;


%% binarization procedure
ic = ceil(param.tile_size*1/8+1);
jc = ceil(param.tile_size*1/8+1);      % position of the center 
            % x/2+1 in case of even dimensions seems to be where center of 
                           %of spectrum after fftshift is!  x/2 ->
                           %non-centered spectrum -> diagonal orientation
                           %artefact    

% get rid of the pixels in the center in a radius
% of obj.param.cut
% because the proportion was defined this way
[I,J] = ndgrid(1:param.tile_size/4,1:param.tile_size/4);
filt = double((I-ic).^2+(J-jc).^2<=(param.cut)^2);
th_spec(filt>0) = 0;


reshaped = reshape(th_spec,[param.tile_size*param.tile_size*1/16 ,1]);  % reshape the spectrum
[values,index]= sort(reshaped,'descend');         % sort in intensity and keep indexes
nbel = param.propor * param.tile_size*param.tile_size*1/16;         % nbel number of points determined by obj.param.propor

if nbel > param.tile_size*param.tile_size/16             % if nbel is to large, take the whole picture
    nbel = param.tile_size*param.tile_size/16;
end
keep = index(1:floor(nbel)); %indices of top propor% birghtest pixels
nbel = size(keep);

%remember avg intensity of kept pixels, as quality score 

%get the normalized intensities of fourier space pixels to keep
%param.noise_threshold;
%filter for those above user-defined noise intensity threshold
%if less than 10% were kept, mark the whole tile as noise
if false
    MS(:,:,re) = [nan, nan; nan, nan]
else

%% mean intensity of region used as quality indicator
quality(re)=mean(th_spec(keep));

th_spec(keep) =1;           % binarize- put to one the points kept by proportion
th_spec(index(nbel+1:end)) = 0;       % put to zero the others

th_spec(filt>0) = 1;                            % get back the pixels in the center in a
% radius of param.cut

%% morphological closing
B = false(param.tile_size/4);                       % create a logical matrix of falses of size tile_size
B(th_spec>0)=1;                              % put rights where th_spec = 1
se=strel('disk',param.strel);                % smooth with a strel
th_spec=imclose(B,se);                       % close image

%% correlations between the binary points as inertia matrix

norm = sum(sum((th_spec).^2));               % computes norm
[Yg,Xg]=meshgrid(1:(param.tile_size/4),1:(param.tile_size/4)); % create meshgrid of correct size

Xo = Xg-ic;
Yo = Yg-jc;
Y = (Yg-ic).^2; % Lignes
X = (Xg-jc).^2; % Colonnes
YY = sum(sum((th_spec).^2.*Y))/norm;         % compute the 3 different values in the inertia matrix
XX = sum(sum((th_spec).^2.*X))/norm;
XY = sum(sum((th_spec).^2.*Xo.*Yo))/norm;
Ms(:,:,re) = [XX,XY;XY,YY];                           % keep matrix of inertia
end
end
end
