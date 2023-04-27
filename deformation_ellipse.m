function [quality,abphi] =deformation_ellipse(param, spectra, nTiles)
% results = def size_analysis(obj.param,results)
% Function that computes the ellipse for the siz representation
% on each averaged subimage,  fills the im_regav structure with the fields
% a, b, phi
%
% *OUTPUT*:+ results
%--------------------------------------------------------------------------

abphi = zeros(3,nTiles);
quality=zeros([1 nTiles]);
for re = 1:nTiles     % for each region
    if sum(sum(spectra(:,:,re)))==0 % if the spectrum is empty
        abphi(:,re)= [0;0;0];
        quality(re)=0;
    else
        [a, b, phi] = size_def(spectra(:,:,re),param);
        quality(re)=mean(spectra(:,:,re), 'all');
        abphi(:,re) = [a b phi];
        % then compute
        % the cell size
        % parameters
    end

end
end

%% Support functions

function [a,b,phi]= size_def(abs_im_fft_w, param)
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
a = floor(param.tileSize/2)*1/ai;
b =floor(param.tileSize/2)*1/bi;
end

function [varargout]=ellipsefit(x,y)
%ELLIPSEFIT Stable Direct Least Squares Ellipse Fit to Data.
% [Xc,Yc,A,B,Phi,P]=ELLIPSEFIT(X,Y) finds the least squares ellipse that
% best fits the data in X and Y. X and Y must have at least 5 data points.
% Xc and Yc are the x- and y-axis center of the ellipse respectively.
% A and B are the major and minor axis of the ellipse respectively.
% Phi is the radian angle of the major axis with respect to the x-axis.
% P is a vector containing the general conic parameters of the ellipse.
% The conic representation of the ellipse is given by:
%
% P(1)*x^2 + P(2)*x*y + P(3)*y^2 + P(4)*x + P(5)*y + P(6) = 0
%
% S=ELLIPSEFIT(X,Y) returns the output data in a structure with field names
% equal to the variable names given above, e.g., S.Xc, S.Yc, S.A, S.B,
% S.Phi and S.P
%
% Reference: R. Halif and J. Flusser, "Numerically Stable Direct Least
% Squares FItting of Ellipses," Department of Software Engineering, Charles
% University, Czech Republic, 2000.

% Conversion from conic to conventional ellipse equation inspired by
% fit_ellipse.m on MATLAB Central

% D.C. Hanselman, University of Maine, Orono, ME 04469
% Mastering MATLAB 7
% 2005-02-28
% Rotation angle fixed 2005-08-09

%--------------------------------------------------------------------------
x=x(:); % convert data to column vectors
y=y(:);
if numel(x)~=numel(y) || numel(x)<5
    error('X and Y Must be the Same Length and Contain at Least 5 Values.')
end
D1=[x.*x x.*y y.*y]; % quadratic terms
D2=[x y ones(size(x))]; % linear terms
S1=D1'*D1;
S2=D1'*D2;

[Q2,R2]=qr(D2,0);
if condest(R2)>1.0e10
    warning('ellipsefit',...
        'Data is Poorly Conditioned and May Not Represent an Ellipse.')
end
T=-R2\(R2'\S2'); % -inv(S3) * S2'

M=S1+S2*T;
CinvM=[M(3,:)/2; -M(2,:); M(1,:)/2];
[V,na]=eig(CinvM);
c=4*V(1,:).*V(3,:) - V(2,:).^2;
A1=V(:,c>0);
P=[A1; T*A1];

% correct signs if needed
P=sign(P(1))*P;

Phi=atan(P(2)/(P(3)-P(1)))/2;
c=cos(Phi);
s=sin(Phi);

% rotate the ellipse parallel to x-axis
Pr=zeros(6,1);
Pr(1)=P(1)*c*c - P(2)*c*s + P(3)*s*s;
Pr(2)=2*(P(1)-P(3))*c*s + (c^2-s^2)*P(2);
Pr(3)=P(1)*s*s + P(2)*s*c + P(3)*c*c;
Pr(4)=P(4)*c - P(5)*s;
Pr(5)=P(4)*s + P(5)*c;
Pr(6)=P(6);

% extract other data
XcYc=[c s;-s c]*[-Pr(4)/(2*Pr(1));-Pr(5)/(2*Pr(3))];
Xc=XcYc(1);
Yc=XcYc(2);
F=-Pr(6) + Pr(4)^2/(4*Pr(1)) + Pr(5)^2/(4*Pr(3));
AB=sqrt(F./Pr(1:2:3));
A=AB(1);
B=AB(2);
Phi=-Phi;
if A<B % x-axis not major axis, so rotate it pi/2
    Phi=Phi-sign(Phi)*pi/2;
    A=AB(2);
    B=AB(1);
end
S.Xc=Xc;
S.Yc=Yc;
S.A=A;
S.B=B;
S.Phi=Phi;
S.P=P;
if nargout==1
    varargout{1}=S;
else
    outcell=struct2cell(S);
    varargout=outcell(1:nargout);
end
end

function varargout = localMaximum_h(x,minDist, exculdeEqualPoints,number)
% function varargout = localMaximum(x,minDist, exculdeEqualPoints)
%
% This function returns the indexes\subscripts of local maximum in the data x.
% x can be a vector or a matrix of any dimension
%
% minDist is the minimum distance between two peaks (local maxima)
% minDist should be a vector in which each argument corresponds to it's
% relevant dimension OR a number which is the minimum distance for all
% dimensions
%
% exculdeEqualPoints - is a boolean definning either to recognize points with the same value as peaks or not
% x = [1     2     3     4     4     4     4     4     4     3     3     3     2     1];
%  will the program return all the '4' as peaks or not -  defined by the 'exculdeEqualPoints'
% localMaximum(x,3)
% ans =
%      4     5     6     7     8     9    11    12
%
%  localMaximum(x,3,true)
% ans =
%      4     7    12
%
%
% Example:
% a = randn(100,30,10);
% minDist = [10 3 5];
% peaks = localMaximum(a,minDist);
%
% To recieve the subscript instead of the index use:
% [xIn yIn zIn] = localMaximum(a,minDist);
%
% To find local minimum call the function with minus the variable:
% valleys = localMaximum(-a,minDist);

if nargin < 3
    exculdeEqualPoints = false;
    if nargin < 2
        minDist = size(x)/10;
    end
end

if isempty(minDist)
    minDist = size(x)/10;
end

dimX = length ( size(x) );
if length(minDist) ~= dimX
    % In case minimum distance isn't defined for all of x dimensions
    % I use the first value as the default for all of the dimensions
    minDist = minDist( ones(dimX,1) );
end

% validity checks
minDist = ceil(minDist);
minDist = max( [minDist(:)' ; ones(1,length(minDist))] );
minDist = min( [minDist ; size(x)] );

% ---------------------------------------------------------------------
if exculdeEqualPoints

    % this section comes to solve the problem of a plato
    % without this code, points with the same hight will be recognized as peaks
    y = sort(x(:));
    dY = diff(y);
    % finding the minimum step in the data
    minimumDiff = min( dY(dY ~= 0) );
    %adding noise which won't affect the peaks
    x = x + rand(size(x))*minimumDiff;
end
% ---------------------------------------------------------------------


se = ones(minDist);
X = imdilate(x,se);
% figure
% imagesc(X(:,:,20))
% figure
% imagesc(x)
f = find(x == X);

[~ ,order] = sort(x(f),'descend');
ind      = f(order(1:min(number,end)));


if nargout
    [varargout{1:nargout}] = ind2sub( size(x), ind);

else
    varargout{1} = f;
end
end