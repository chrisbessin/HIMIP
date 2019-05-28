
function write_envihdr(info, filereflhdr)
% Write header file for ENVI file
% Example: 
%    write_envihdr(info, filename.hdr);
% where info has the same format as the one obtained from the read_envihdr
% function

% Author: Thanh Bui (thanh.bui@erametgroup.com)

fid = fopen(filereflhdr, 'w');
fprintf(fid,'%s\n', 'ENVI');
fprintf(fid,'%s\n','description = {');
fprintf(fid,'%s\n','Exported from MATLAB}');

fprintf(fid,'\n%s%d', 'samples = ', info.samples);
fprintf(fid,'\n%s%d', 'bands = ', info.bands);
fprintf(fid,'\n%s%d', 'lines = ', info.lines);

if(exist('info.header_offset'))
    fprintf(fid,'\n%s%d', 'header offset = ', info.header_offset);
end
if(exist('info.file_type'))
    fprintf(fid,'\n%s%s', 'file type = ', info.file_type);
end
fprintf(fid,'\n%s%d', 'data type = ', 4);
fprintf(fid,'\n%s%s', 'interleave = ', info.interleave);
if(exist('info.sensor_type'))
    fprintf(fid,'\n%s%s', 'sensor type = ', info.sensor_type);
end
fprintf(fid,'\n%s%d', 'byte order = ', info.byte_order);
fprintf(fid,'\n%s%s', 'default bands = ', '{70, 165, 237}');
fprintf(fid,'\n%s%s', 'wavelength units = ', 'Nanometers');
fprintf(fid,'\n%s%s', 'z plot title = ', '{nm, Reflectance}');
fprintf(fid,'\n\n%s', 'Wavelength = {');
for i = 1:length(info.Wavelength)
    if(i ==1)
        fprintf(fid,'\n%5.2f', info.Wavelength(i));
    else
        fprintf(fid,'%s\n%5.2f', ',', info.Wavelength(i));
    end
end
fprintf(fid,'\n%s', '}');

if(exist('info.fwhm'))
    fprintf(fid,'\n\n%s', 'fwhm = {');
    for i = 1:length(info.fwhm)
        if (i == 1)
            fprintf(fid,'\n%5.2f', info.fwhm(i));
        else
            fprintf(fid,'%s\n%5.2f', ',', info.fwhm(i));
        end
    end
end
fprintf(fid,'\n%s', '}');

fclose(fid);




