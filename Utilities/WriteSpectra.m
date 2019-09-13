function WriteSpectra (filename, imgFile, spectra, spectraIndices, saveAll)
% Write spectra to a .txt file complied with ENVI spectra file format

% Author: Thanh Bui (thanh.bui@erametgroup.com)

fid = fopen(filename, 'wt');
fprintf(fid, 'File %s, %s\n', imgFile , datestr(datetime('now')));
fprintf(fid, '%s\n', 'Column 1: Wavelength');

% Write the header
if saveAll % Save all spectra in the file
    startIndex = 2;
    for i = startIndex: length(spectra)
        fprintf(fid, '%s %d: %s\n', 'Column', i, spectra(i).pos );
    end
else % save the selected spectra specified in spectraIndices
    for i = 1:length(spectraIndices)
        fprintf(fid, '%s %d: %s\n', 'Column', spectraIndices(i), spectra(spectraIndices(i)).pos );
    end
end
       
% Write the wavelength and the spectra
for i = 1: length(spectra(2).wavelength)   
    fprintf(fid, '  %5.3f\t', spectra(2).wavelength(i));
    if saveAll
        for ii = startIndex: length(spectra)
            fprintf(fid, '  %5.6f\t', spectra(ii).spectrum(i));
        end
        fprintf(fid, '\n');
    else
        for ii = 1:length(spectraIndices)
            fprintf(fid, '  %5.6f\t', spectra(spectraIndices(ii)).spectrum(i));
        end
        fprintf(fid, '\n');
    end
end 
fclose(fid);