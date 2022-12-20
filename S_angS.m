function [S,angS] = S_angS(param, spectrum)
            % [S,angS] = S_angS(spectrum)
            % This function defines the angle and amplitude of cell deformation
            % from inertia matrix

            %--------------------------------------------------------------------------
            %
            % *INPUT*: spectrum Fourier space spectrum
            %
            % *OUTPUT*: + S the amplitude of the cell deformation defined by the
            %           inertia matrix
            %           + angS angle of the cell deformation defined by the inertia
            %           matrix
            %
            %--------------------------------------------------------------------------
            M=inertia_matp_sigma(param, spectrum);

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

function [M] = inertia_matp_sigma(param, al)
% M = inertia_matp_sigma(al,obj.param)
% function that computes the inertia matrix based on a spectrum averaged
% smoothed with a gaussian filter and with a percentile of points to keep.
%
% *OUTPUT*: + M the inertia matrix
%--------------------------------------------------------------------------

%Blurring the Fourier space image
if ~isempty(param.sigma) && param.sigma>0
    th_spec = imgaussfilt(al,param.sigma);
else
    th_spec = al;
end
% try not doing this - it's equivalent to real space pointwise
% multipication with a Gaussian, darkens central real space
% area

ic = param.tile_size*1/8+1;
jc = param.tile_size*1/8+1;          % get rid of the pixels in the center in a radius
% of obj.param.cut
% because the proportion was defined this way
[I,J] = ndgrid(1:param.tile_size/4,1:param.tile_size/4);
M = double((I-ic).^2+(J-jc).^2<=(param.cut)^2);
th_spec(M>0) = 0;

%also get rid of x and y vertical lines if applicable
if param.stripe
    dims=size(th_spec);
    Mask = zeros(dims(1), dims(2), "logical");
    Mask(round((dims(1)-param.stripewidth)/2):round((dims(1)+param.stripewidth)/2),:) = true;
    Mask(:,round((dims(2)-param.stripewidth)/2):round((dims(2)+param.stripewidth)/2)) = true;
    th_spec(Mask>0) = 0;
end


reshaped = reshape(th_spec,[param.tile_size*param.tile_size*1/16 ,1]);  % reshape the spectrum
[values,index]= sort(reshaped,'descend');         % sort in intensity and keep indexes
nbel = param.propor * param.tile_size*param.tile_size*1/16;         % nbel number of points determined by obj.param.propor

if nbel > param.tile_size*param.tile_size/16             % if nbel is to large, take the whole picture
    nbel = param.tile_size*param.tile_size/16;
else
end
keep = index(1:floor(nbel)); %indices of top propor% birghtest pixels
keep = keep(values(1:floor(nbel))>param.threshold); %keep only those with intensity above noise threshold
nbel = size(keep);
th_spec(keep) = 1;           % put to one the points kept by proportion
th_spec(index(nbel+1:end)) = 0;       % put to zero the others

th_spec(M>0) = 1;                            % get back the pixels in the center in a
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
M = [XX,XY;XY,YY];                           % keep matrix of inertia
end