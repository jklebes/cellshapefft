function [Si,angSi] = S_angS(M)
% [S,angS] = S_angS(M)
% This function defines the angle and amplitude of cell deformation
% from inertia matrix

%--------------------------------------------------------------------------
%
% *INPUT*: M inertia matrix
%
% *OUTPUT*: + S the amplitude of the cell deformation defined by the 
%           inertia matrix
%           + angS angle of the cell deformation defined by the inertia
%           matrix
%
%--------------------------------------------------------------------------

    Si = zeros(size(M,1),1);
    angSi = zeros(size(M,1),1);
    for kk = 1:size(M,1)
        Mo=squeeze(M(kk,:,:));
        [~,E] = eig(Mo);                % extract the eigenvalues of the inertia matrix
        L1 = sqrt(E(2,2)/2);         % definition of the minor and major axis of the inertia ellipse     
        L2 = sqrt(E(1,1)/2);
        lgE = -1/2*(log(L1/L2));       % amplitude of the cell deformation
        Dev = Mo-1/2*trace(Mo)*eye(2);   % deviation matrix
        [V,~] = eig(Dev);              % eigenvectors of deviation matrix to extract the angle
        ang=atan(V(1,2)/V(2,2));       % definition of angle
        angSi(kk) = ang;                    % register angle
        Si(kk) = lgE;                       % register amplitude
    end
end