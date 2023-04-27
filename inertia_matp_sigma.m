function [quality,Ms] = inertia_matp_sigma(param, spectra, nTiles)
% M = inertia_matp_sigma(al,obj.param)
% function that computes the inertia matrix based on a spectrum averaged
% smoothed with a gaussian filter and with a percentile of points to keep.
%
% *OUTPUT*: + M the inertia matrix
%--------------------------------------------------------------------------
Ms = zeros(2,2,nTiles);
quality=zeros([1 nTiles]);

for re=1:nTiles
th_spec=spectra(:,:,re);

%% binarization procedure
    ic = param.tileSize*1/2+1;
    jc = param.tileSize*1/2+1;          
            % x/2+1 in case of even dimensions seems to be where center of 
                           %of spectrum after fftshift is!  x/2 ->
                           %non-centered spectrum -> diagonal orientation
                           %artefact    

% get rid of the pixels in the center in a radius
% of obj.param.cut
% because the proportion was defined this way
[I,J] = ndgrid(1:param.tileSize,1:param.tileSize);
filt = double((I-ic).^2+(J-jc).^2<=(param.cut)^2);
%th_spec(filt>0) = 0;

%% Boolean matrix used after and strel
     se=strel('disk',param.strel);                % smooth with a strel
th_spec(filt>0) = 0;

reshaped = reshape(th_spec,[param.tileSize*param.tileSize ,1]);  % reshape the spectrum
[~,index]= sort(reshaped,'descend');                            % sort in intensity and keep indexes
nbel = param.propor * param.tileSize*param.tileSize;         % nbel number of points determined by obj.param.propor

%remember avg intensity of kept pixels, as quality score 
%% mean intensity of region used as quality indicator
quality(re)=mean(index(1:floor(nbel)));

        th_spec(index(1:floor(nbel))) = 1;           % put to one the points kept by proportion
        th_spec(index(floor(nbel)+1:end)) = 0;       % put to zero the others
        th_spec(filt>0) = 1;                            % get back the pixels in the center
                           % put rights where th_spec = 1
        th_spec=imclose(th_spec,se);                       % close image
        
        %% Inertia matrix computation
        
        norm = sum(sum((th_spec).^2));               % computes norm 
        [Yg,Xg]=meshgrid(1:(param.tileSize),1:(param.tileSize)); % create meshgrid of correct size
        xc = param.tileSize/2;                               % position of the center 
        yc = param.tileSize/2;
        Xo = Xg-xc;                           
        Yo = Yg-yc;
        Y = (Yg-yc).^2; % Lignes
        X = (Xg-xc).^2; % Colonnes
        YY = sum(sum((th_spec).^2.*Y))/norm;         % compute the 3 different values in the inertia matrix
        XX = sum(sum((th_spec).^2.*X))/norm;
        XY = sum(sum((th_spec).^2.*Xo.*Yo))/norm;
        
        %% Registering of inertia matrix

Ms(:,:,re) = [XX,XY;XY,YY];                           % keep matrix of inertia
end
end
