function SangS = S_angS(Ms, regl)
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
SangS=zeros(2,regl);
for re=1:regl
    M=Ms(:,:,re);
    [~,E] = eig(M);                % extract the eigenvalues of the inertia matrix
    L1 = 2*sqrt(E(2,2)/2);         % definition of the minor and major axis of the inertia ellipse
    L2 = 2*sqrt(E(1,1)/2);
    SangS(1,re) = 1/2*(log(L1/L2));       % amplitude of the cell deformation
    Dev = M-1/2*trace(M)*eye(2);   % deviation matrix
    [V,~] = eig(Dev);              % eigenvectors of deviation matrix to extract the angle
    SangS(2,re)=atan(V(1,2)/V(2,2));       % definition of angle
end
end