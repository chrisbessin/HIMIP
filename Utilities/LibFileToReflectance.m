function [reflsm, refl, wavelength] = LibFileToReflectance(filePath, SkipNoLines)
order = 2;
framelen = 5;

data = ReadTextFile(filePath, SkipNoLines);
wavelength_temp = data(:,1);
refl_temp = data(:,2);
if(mean(wavelength_temp) < 10)
    wavelength_temp = wavelength_temp*1000;
end
index = find(wavelength_temp >= 1000 & wavelength_temp <=2500);

wavelength = wavelength_temp(index);
refl = refl_temp(index);
length(wavelength);
% Reduce noise
reflsm = sgolayfilt(refl, order, framelen);