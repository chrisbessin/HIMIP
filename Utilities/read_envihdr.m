function info = read_envihdr(hdrfile)
% read_envihdr reads and return ENVI image file header information.
%   info = read_envihdr('hdr_file.hdr') reads the ASCII ENVI-generated image
%   header file and returns all the information in a structure of
%   parameters.
% 
%   Example:
%   >> info = read_envihdr('my_envi_image.hdr')
% info = 
% description: '{File Imported into ENVI}'
%                    file_type: 'ENVI'
%                  sensor_type: 'SWIR , Lumo - Recorder v2016-424'
%             acquisition_date: 'DATE(yyyy-mm-dd): 2017-10-31'
%                   Start_Time: 'UTC TIME: 10:16:35'
%                    Stop_Time: 'UTC TIME: 10:17:19'
%                      samples: 384
%                        bands: 288
%                        lines: 880
%                       errors: '{none}'
%                   interleave: 'bil'
%                    data_type: 12
%                header_offset: 0
%                   byte_order: 0
%                      x_start: 0
%                      y_start: 0
%                default_bands: '{72, 160, 127}'
%                         himg: '{1, 384}'
%                         vimg: '{1, 288}'
%                         hroi: '{1, 384}'
%                         vroi: '{1, 288}'
%                          fps: 20
%                      fps_qpf: 20
%                         tint: 2.9998
%                      binning: '{1, 1}'
%                 trigger_mode: 'Internal'
%                 trigger_sync: 1
%                        fodis: '{0, 0}'
%                     sensorid: 431048
%       acquisitionwindow_left: 0
%        acquisitionwindow_top: 0
%             calibration_pack: 'C:/Users/Public/Documents/Specim/431048_20170223_OLES56.scp'
%             SWIR_temperature: 151
%     Scb_temperature_channel1: 28.6800
%     Scb_temperature_channel2: 26.2100
%     Scb_temperature_channel3: 26.9200
%     Scb_temperature_channel4: 26.4000
%                  temperature: [5×1 double]
%                   Wavelength: [288×1 double]
%                         fwhm: [288×1 double]
%
%   Author: Thanh Bui (thanh.bui@erametgroup.com), refer to the one written by Ian M. Howat (University
%   of Washington)

fid = fopen(hdrfile);
while fid
    line = fgetl(fid);
    if line == -1
        break
    else
        eqsn = findstr(line,'=');
        if ~isempty(eqsn)
            param = strtrim(line(1:eqsn-1));
            param(findstr(param,' ')) = '_';
            value = strtrim(line(eqsn+1:end));
            if isempty(str2num(value))
                if ~isempty(findstr(value,'{')) && isempty(findstr(value,'}'))
                    while isempty(findstr(value,'}'))
                        line = fgetl(fid);
                        value = [value,strtrim(line)];
                    end
                end
                eval(['info.',param,' = ''',value,''';'])
            else
                eval(['info.',param,' = ',value,';'])
            end
        end
    end
end
fclose all;

if isfield(info,'map_info')
    line = info.map_info;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    info.map_info = [];
    info.map_info.projection = line{1};
    info.map_info.image_coords = [str2num(line{2}),str2num(line{3})];
    info.map_info.mapx = str2num(line{4});
    info.map_info.mapy = str2num(line{5});
    info.map_info.dx  = str2num(line{6});
    info.map_info.dy  = str2num(line{7});
    if length(line) == 9
        info.map_info.datum  = line{8};
        info.map_info.units  = line{9}(7:end);
    elseif length(line) == 11
        info.map_info.zone  = str2num(line{8});
        info.map_info.hemi  = line{9};
        info.map_info.datum  = line{10};
        info.map_info.units  = line{11}(7:end);
    end
end
if isfield(info, 'Wavelength')
    line = info.Wavelength;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    wavelength = zeros(2,1);
    for i = 1:length(line)
        wavelength(i) = str2double(line{i});
    end
    info.Wavelength = wavelength;       
end
if isfield(info, 'fwhm')
    line = info.fwhm;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    fwhm = zeros(2,1);
    for i = 1:length(line)
        fwhm(i) = str2double(line{i});
    end
    info.fwhm = fwhm;
end
if isfield (info, 'temperature')
    line = info.temperature;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    temperature = zeros(2,1);
    for i = 1:length(line)
        temperature(i) = str2double(line{i});
    end
    info.temperature = temperature;
end
if isfield(info,'pixel_size')
    line = info.pixel_size;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    info.pixel_size = [];
    info.pixel_size.x = str2num(line{1});
    info.pixel_size.y = str2num(line{2});
    info.pixel_size.units = line{3}(7:end);
end

function A = split(s,d)
%This function by Gerald Dalley (dalleyg@mit.edu), 2004
A = {};
while (~isempty(s))
    [t,s] = strtok(s,d);
    A = {A{:}, t};
end




