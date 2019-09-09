function [data, info, rgb_img] = access_spectra_data(dataFile)

% Read hyperspectral data and return data, info and rgb image of the data

% Author: Thanh Bui (thanh.bui@erametgroup.com)

hdrFile = strcat(dataFile(1:end-4), '.hdr');
info = read_envihdr(hdrFile);
wavelength = info.Wavelength;
lines = info.lines;
samples = info.samples;
bands = info.bands;
precision = datatype2precision(info.data_type);
data = multibandread(dataFile, [lines, samples, bands], precision,0 , 'bil','ieee-le'); 

if(mean(wavelength > 1100))
    R_WL = 2000; G_WL = 2200; B_WL = 2350; 
else
    R_WL = 700; G_WL = 600; B_WL = 500; 
end
rgb_img = hyperspectraldata2rgbimg(data, wavelength, R_WL, G_WL, B_WL);
