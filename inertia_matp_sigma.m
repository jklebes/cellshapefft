function [Ms, quality] = inertia_matp_sigma(param, spectra, regl)
% M = inertia_matp_sigma(al,obj.param)
% function that computes the inertia matrix based on a spectrum averaged
% smoothed with a gaussian filter and with a percentile of points to keep.
%
% *OUTPUT*: + M the inertia matrix
%--------------------------------------------------------------------------
Ms = zeros(2,2,regl);
quality=zeros([regl 1]);
for re=1:regl
al=spectra(:,:,re);
%Blurring the Fourier space image
if ~isempty(param.sigma) && param.sigma>0
    th_spec = imgaussfilt(al,param.sigma);
else
    th_spec = al;
end
% it's equivalent to real space pointwise
% multipication with a Gaussian, darkens central real space
% area

ic = param.tile_size*1/8+1;
jc = param.tile_size*1/8+1;          % get rid of the pixels in the center in a radius
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
else
end
keep = index(1:floor(nbel)); %indices of top propor% birghtest pixels
nbel = size(keep);

%remember avg intensity of kept pixels, as quality score 
quality(re) = mean(values(keep));

th_spec(keep) = 1;           % binarize- put to one the points kept by proportion
th_spec(index(nbel+1:end)) = 0;       % put to zero the others

th_spec(filt>0) = 1;                            % get back the pixels in the center in a
% radius of param.cut

B = false(param.tile_size/4);                       % create a logical matrix of falses of size tile_size
B(th_spec>0)=1;                              % put rights where th_spec = 1
se=strel('disk',param.strel);                % smooth with a strel
th_spec=imclose(B,se);                       % close image

norm = sum(sum((th_spec).^2));               % computes norm
[Yg,Xg]=meshgrid(1:(param.tile_size/4),1:(param.tile_size/4)); % create meshgrid of correct size
xc = param.tile_size/8;                               % position of the center
yc = param.tile_size/8;
Xo = Xg-xc;
Yo = Yg-yc;
Y = (Yg-yc).^2; % Lignes
X = (Xg-xc).^2; % Colonnes
YY = sum(sum((th_spec).^2.*Y))/norm;         % compute the 3 different values in the inertia matrix
XX = sum(sum((th_spec).^2.*X))/norm;
XY = sum(sum((th_spec).^2.*Xo.*Yo))/norm;
Ms(:,:,re) = [XX,XY;XY,YY];                           % keep matrix of inertia
end
end