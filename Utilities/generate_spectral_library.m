function generate_spectral_library(spectraPath, spectralLibFile)

% Generate the spectral library from files containing in a specified path

% Author: Thanh Bui (thanh.bui@erametgroup.com)

addpath('D:\Matlab\Utilities\')
fid = fopen(spectralLibFile, 'w');
%% Generate spectral library A (Lxm): L bands, m endmembers
listdir = dir(spectraPath);
m = length(listdir) - 2;
SkipNoLines = 3;
filePath = fullfile(spectraPath, listdir(5).name);
[reflsm, refl, wavelength_temp] = LibFileToReflectance(filePath, SkipNoLines);
WaveBottom = 1050; WaveTop = 2450;
idx = find(wavelength_temp >=WaveBottom & wavelength_temp <=WaveTop);
A = zeros(length(idx), m);
min_names = {m};

wavelength_ref = wavelength_temp(idx);


for i = 1: m% length(listdir)
    filePath = fullfile(spectraPath, listdir(i+2).name)
    [reflsm_temp, refl_temp, wavelength_temp] = LibFileToReflectance(filePath, SkipNoLines);
    idx = find(wavelength_temp >=WaveBottom & wavelength_temp <=WaveTop);
    if(length(idx) < length(wavelength_temp))
        wavelength = wavelength_temp(idx);
        reflsm = reflsm_temp(idx);
        refl = refl_temp(idx);
    else
        wavelength = wavelength_temp;
        reflsm = reflsm_temp;
        refl = refl_temp;
    end
    
    % Interpolation with csapi
    if (length(wavelength) ~= length(wavelength_ref))
        reflsm =  csapi (wavelength, reflsm, wavelength_ref);
        refl =  csapi (wavelength, refl, wavelength_ref);
        wavelength = wavelength_ref;
    end
    if(max(reflsm) > 10)
        reflsm = reflsm/100;
        refl = refl/100;
    end
    
    fprintf(fid,'%s\n', listdir(i+2).name);
    A(:,i) = refl;
    min_names{i} = listdir(i+2).name;

end
fclose(fid);
save (spectralLibFile, 'A', 'min_names', 'wavelength')

figure, imagesc(A), xlabel('# endmembers'), ylabel('# bands')
