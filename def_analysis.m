function im_regav_c = def_analysis(param, spectra, regl,  c)
% im_regav_c = def analysis(obj.param,results)
% deformation analysis, inertia matrix route.
% Function that computes the inertia matrix and cell deformation on
% each averaged subimage, fills the im_regav structure with the fields S
% and angS.
%
% *OUTPUT*: a single time slice c of im_regav final results structure, with
% fields S and S_ang filled in
%--------------------------------------------------------------------------


display(['Computation of cell deformation, time: ',num2str(param.time_points(c))]); % for the user to keep track
im_regav_c = repmat(struct('spect', [], 'S', [], 'angS', [], 'a', [], 'b', [], 'phi', []), regl, 1);
for re = 1:regl    % for each region
    if all(spectra(:,:,re)==0)  % if the spectrum is empty
        im_regav_c(re).angS = 0;         % don't register cell deformation
        im_regav_c(re).S = 0;

    else 
        % compute the cell deformation amplitude and angle via inertia
        % matrix
        spectrum=spectra(:,:,re);
        [im_regav_c(re).S,im_regav_c(re).angS] = S_angS(param, spectrum);
        if ~(param.register == 0 || param.regsize == 0)
            im_regav_c(re).spect = spectrum;
        end
    end

    if param.register == 0 && param.regsize == 0 %TODO ??
        im_regav_c(re).spect = [];       % delete the spectrum to keep space
    else
    end
end

end