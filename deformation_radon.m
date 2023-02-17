function [outputArg1,outputArg2] = deformation_radon(inputArg1,inputArg2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% after Streichan et al 2018

%do radon transform
r=radon(im)
%identify N peaks ( extended Maxima transform)
BW = imextendedmax(r,n_lines)
%real space lines
%determine relevant segment of real space line, length
%merge sufficiently similar lines

%total/average instensity along lines

%return - tile avg orientation and intensity? all the lines?
end