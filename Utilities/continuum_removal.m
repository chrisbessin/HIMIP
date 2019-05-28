function [refl_cr, refl_hull] = continuum_removal (wavelength, refl)

% This function performs hull corrections for the input reflectance, i.e., remove the background
% Input:
%     wavelength (n x 1): wavelengths of the reflectance in nm
%     refl (n x 1): the relative reflectance of hyperspectral data
% Output:
%     refl_cr: the continuum removal of the input reflectance
%     refl_hull: the background (hull) curve
% @Author: Thanh Bui
% Example: 
% Given a hyperspectral reflectance data Spectra (n,m, number of bands) and
% the wavelengthSWIR
% figure, plot(wavelengthSWIR, squeeze(Spectra(100,100,:)))
% figure, plot(squeeze(Spectra(100,100,:)))
% [cr, hull] = continuum_removal(wavelengthSWIR, squeeze(Spectra(100,100,:)));
% figure, plot(wavelengthSWIR, cr)


refl = refl(:);
% Remove noise 
order = 2;
framelen = 11;
refl = sgolayfilt(refl, order, framelen);
% Compute hull curve and hull removal 
[hull,~] = convhull(wavelength, refl);
start_idx = find(hull == length(refl));
upper_hull = hull(end:-1:start_idx);
refl_hull = interp1(wavelength(upper_hull),refl(upper_hull), wavelength);
refl_cr = refl./(refl_hull +eps);


