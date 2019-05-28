
function [rgb_img] = create_rgb_img(refl, wavelength, R_WL, G_WL, B_WL)

resolution = (max(wavelength)-min(wavelength))/length(wavelength);
R_index = round((R_WL-min(wavelength))/resolution);
G_index = round((G_WL-min(wavelength))/resolution);
B_index = round((B_WL-min(wavelength))/resolution);
R = refl(:,:,R_index); G = refl(:,:,G_index); B = refl(:,:,B_index);

rgb_img = cat(3,R,G,B);