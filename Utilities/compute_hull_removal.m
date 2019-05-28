function [refl_cr, wavelength, rgb_img_cr] = compute_hull_removal(reflFile, reflCRFile)

% Read ENVI reflectance file and compute hull correction data
if ~exist('reflCRFile', 'var') % if reflFile is not an argument   
    reflCR_file = strcat(reflFile(1:end-4), 'CR', sampleFile(end-3:end));
    reflCR_hdr_file = strcat(reflFile(1:end-4), 'CR.hdr');
else  
    reflCR_file = reflCRFile;
    reflCR_hdr_file = strcat(reflCRFile(1:end-4), '.hdr');
end
% Read header file
refl_hdr = strcat(reflFile(1:end-4), '.hdr');
if (exist(reflCR_file, 'file'))   %the reflectance was already computed
    fprintf('The continuum removal existed \n')
    info = read_envihdr(reflCR_hdr_file);
    wavelength = info.Wavelength;
    lines = info.lines;
    samples = info.samples;
    bands = info.bands;
    refl = multibandread(reflFile, [lines, samples, bands],'single',0, 'bil','ieee-le' );
    refl_cr = multibandread(reflCR_file, [lines, samples, bands],'single',0, 'bil','ieee-le' ); 
else        
    info = read_envihdr(refl_hdr);
    lines = info.lines;
    samples = info.samples;
    bands = info.bands;
    wavelength = info.Wavelength;
    % Read the white, dark and sample data
    refl = multibandread(reflFile, [lines, samples, bands],'single',0, 'bil','ieee-le' );
    
    % Create an waitbar
    f = waitbar(0,'1','Name','Computing continuum removal ...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(f,'canceling',0);
    
    NN = size(refl);
    refl_cr = zeros(NN);
    for i = 1: NN(1)
        for ii = 1:NN(2)
            [refl_cr(i,ii,:), temp] = ContinuumRemovalSeg(wavelength, refl(i,ii,:), 10);
        end
        % Escape when user clicks Cancel
        if getappdata(f,'canceling')
            break
        end
        % Update waitbar and message
        waitbar(i/NN(1),f,sprintf('%d%%',int16(i*100/NN(1))))
    end
    
    multibandwrite(refl_cr, reflCR_file, 'bil', 'precision', 'single'); % specify the precision is very important
    write_envihdr(info, reflCR_hdr_file);
    delete(f)
    
end
% Display the false-color image of the sample
%R_WL = 1250; G_WL = 1500; B_WL = 1750; 
if(mean(wavelength > 1100)) % SWIR data
    R_WL = 2000; G_WL = 2200; B_WL = 2350; 
else                        % VNIR data
    R_WL = 700; G_WL = 600; B_WL = 500; 
end
Resolution = (max(wavelength)-min(wavelength))/length(wavelength);
R_index = round((R_WL-min(wavelength))/Resolution);
G_index = round((G_WL-min(wavelength))/Resolution);
B_index = round((B_WL-min(wavelength))/Resolution);

R = refl(:,:,R_index); G = refl(:,:,G_index); B = refl(:,:,B_index);
rgb_img = cat(3,R,G,B);

R = refl_cr(:,:,R_index); G = refl_cr(:,:,G_index); B = refl_cr(:,:,B_index);
rgb_img_cr = cat(3,R,G,B);



figure, subplot 121, imshow(rgb_img), title('False color image')
subplot 122, imshow(rgb_img_cr), title('Continuum removal image')
%imwrite(rgb_img_norm, refl_rgb_file);
   
end
    



