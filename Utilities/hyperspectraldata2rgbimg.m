
function rgb_img = hyperspectraldata2rgbimg(data, wavelength, R_WL, G_WL, B_WL)
% Create an rgb image from hyperspectral data with specified wavelengths

% Author: Thanh Bui (thanh.bui@erametgroup.com)

resolution = (max(wavelength)-min(wavelength))/length(wavelength);
r_index = round((R_WL-min(wavelength))/resolution);
g_index = round((G_WL-min(wavelength))/resolution);
b_index = round((B_WL-min(wavelength))/resolution);

if r_index < 1, r_index = 1; end
if g_index < 1, g_index = 1; end
if b_index < 1, b_index = 1; end

R = data(:,:,r_index); G = data(:,:,g_index); B = data(:,:,b_index);
rgb_img = cat(3,R,G,B);
% rgb_img_norm = rgb_img;
% rgb_img_norm(:,:,1) = normalize(rgb_img_norm(:,:,1));
% rgb_img_norm(:,:,2) = normalize(rgb_img_norm(:,:,2));
% rgb_img_norm(:,:,3) = normalize(rgb_img_norm(:,:,3));
rgb_img = double(rgb_img/max(rgb_img(:)));