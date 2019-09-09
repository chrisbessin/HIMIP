function [refl_cr, refl_hull, seg_idx] = ContinuumRemovalSeg (wavelength, refl, NB)

% This function performs hull corrections for the input reflectance, i.e., remove the background, and 
% find indices of each extracted segment (start - end)
% Input:
%     wavelength (n x 1): wavelengths of the reflectance in nm
%     refl (n x 1): the relative reflectance of hyperspectral data
%     NB (integer): minimum number of bands for a segment to be extracted
% Output:
%     refl_cr: the continuum removal of the input reflectance
%     refl_hull: the background (hull) curve
%     seg_idx: the indices of extracted segments (one segment for each row: start index - end index)
% Author: Thanh Bui (thanh.bui@erametgroup.com)
% Example: 
% Given a hyperspectral reflectance data Spectra (n,m, number of bands) and
% the wavelength
% figure, plot(wavelength, squeeze(Spectra(100,100,:)))
% figure, plot(squeeze(Spectra(100,100,:)))
% [cr, hull, seg] = ContinuumRemovalSeg(wavelength, squeeze(Spectra(100,100,:)), 6);
% figure, plot(wavelength, cr)

refl = refl(:);
% Remove noise 
order = 2;
framelen = 5;
refl = sgolayfilt(refl, order, framelen);
% Compute hull curve and hull removal 
[hull,~] = convhull(wavelength, refl);
start_idx = find(hull == length(refl));
up_hull = hull(end:-1:start_idx);
refl_hull = interp1(wavelength(up_hull), refl(up_hull), wavelength);
refl_cr = refl./(refl_hull +eps);
% Compute seg_idx
temp = [100; diff(up_hull)];
idx = find(temp>NB);  % for keeping the segments with the number of bands > NB
idx1 = idx-1;
seg_idx = [up_hull(idx1(2:end)) up_hull(idx(2:end))];


