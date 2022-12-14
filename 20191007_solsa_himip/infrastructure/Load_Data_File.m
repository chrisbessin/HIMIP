classdef Load_Data_File

    properties 
        data_names
        data_loaded
    end
        
    methods
        function obj = Load_Data_File(dir_path_acquisition,dir_path_fusion,type_file)
            % Load the data of any type (RGB, VNIR, SWIR, XRF)
            %
            % Parameters
            % ----------
            % dir_path_acquisition: string
            %       path of the acquisition directory
            % dir_path_fusion: string
            %       path of the fusion directory
            % type_file: string
            %       type of data to load (rgb, vnir, swir, xrf)
            %
            % Output
            % ------
            % data_names: array of string (1xr)
            %       names of the fields loaded
            % data_loaded: array of numbers (pxqxr)
            %       data loaded of the image (pxq) containing r fields
            %
            
            global cfg

            home_dir = pwd;

            if strcmp(type_file,'rgb')
                % Find file
                cd(fullfile(dir_path_fusion,'S2P','RGB'))
                files = dir('*.tif*');
                info_img = imfinfo(files(1).name); 

                % Load data
                data_loaded = imread(files(1).name,'tiff','Info',info_img);
                data_loaded = cast(data_loaded,'double')./max(info_img.MaxSampleValue);

                % Remove the 4th channel if it exists
                if size(data_loaded,3) > 3
                    data_loaded = data_loaded(:,:,1:3);
                end

                % 90° rotation to obtain the appropriate orientation
                data_loaded = rot90(data_loaded,-1);

                data_names = [{'RGB_N_R'},{'RGB_N_G'},{'RGB_N_B'}];

            elseif strcmp(type_file,'vnir') || strcmp(type_file,'swir')
                tag = upper(type_file);

                % Find file
                cd(fullfile(dir_path_acquisition,tag))
                files = dir('*.raw');
                i = 1;
                while i <= length(files)
                   if ~contains(files(i).name,cfg.darkref_tag) && ~contains(files(i).name, cfg.whiteref_tag)
                       file = files(i).name(1:end-4);
                       i = Inf;
                   end
                   i = i + 1;
                end

                st_t = tic;
                
                % Load data
                [info, data_loaded] = ProcessEnviDarkWhiteref(fullfile(dir_path_acquisition,tag),file);
                
                % data names 
                data_names = strcat(tag,'_C_',strrep(strsplit(num2str(info.Wavelength')),'.','_'));
                
                if cfg.print_times
                    fprintf('\tLoad_file ')
                    toc(st_t)
                end

            elseif strcmp(type_file,'xrf')
                % TODO when the xrf files will be available
                error('Format fo XRF not known or coded')
            else
                error('wrong file type')
            end
            cd(home_dir)
            
            % Properties 
            obj.data_names = data_names;
            obj.data_loaded = data_loaded;
        end
    end
end

function [darkref_path, whiteref_path] = find_ref_path(dir_sample,darkref_tag,whiteref_tag)
    % Finds the path of a reference file using a tag
    %
    % Parameters
    % ----------
    % dir_sample: string
    %       path of the data file
    % darkref_tag: string
    %       tag of the dark reference file
    % whiteref_tag: string
    %       tag of the white reference file
    %
    % Output
    % ------
    % darkref_path: string
    %       path of the dark reference file
    % whiteref_path: string
    %       path of the white reference file
    %
    home_dir = pwd;
    cd(dir_sample);
    
    % Dark ref
    darkref_file = dir(strcat(darkref_tag,'*.raw'));
    if ~isempty(darkref_file)   
        % Takes the first file found if more reference files are found
        darkref_path = fullfile(dir_sample,darkref_file(1).name);
    end
    
    % White ref
    whiteref_file = dir(strcat(whiteref_tag,'*.raw'));
    if ~isempty(whiteref_file)   
        % Takes the first file found if more reference files are found
        whiteref_path = fullfile(dir_sample,whiteref_file(1).name);
    else
        % Case where the whiteref file is outside the data directory
        cd ..
        % Search the folder of the whiteref
        curr_path = pwd;
        parts = strsplit(curr_path, '\');
        current_fold = parts{end};
        split = strsplit(current_fold,'20');
        tag_data = split{1};
        cd ..
        
        dir_ref  = dir(strcat(tag_data,whiteref_tag,'*'));
        if isempty(dir_ref)
            dir_ref  = dir(strcat(whiteref_tag,tag_data,'*'));
        end
        if isempty(dir_ref)
            dir_ref  = dir(strcat(tag_data,'*',whiteref_tag(1:end-1)));
        end
        
        % Set the folder of the whiteref (if more than one, takes the first one)
        try
            whiteref_dir = fullfile(pwd,dir_ref(1).name,'capture');
            cd(whiteref_dir);
        catch
            error("Couldn't find the whiteref folder")
        end
        
        % Search of the whiteref file
        files_ref = dir(strcat('*.raw'));
        for f = 1:length(files_ref)
           if ~contains(files_ref(f).name,darkref_tag)
               whiteref_file = files_ref(f).name;
           end
        end
        
        % Set of the whiteref file
        if ~isempty(whiteref_file)
            whiteref_path = fullfile(pwd,whiteref_file);
        end
    end
    
    cd(home_dir);
end

function [info, data] = ProcessEnviDarkWhiteref(dir_sample,sample_tag)
global cfg

% Input files
data_path = fullfile(dir_sample, strcat(sample_tag, '.raw'));
data_path_header = fullfile(dir_sample, strcat(sample_tag, '.hdr'));

% get info from Header
info = ReadEnviHdr(data_path_header);
if info.data_type == 12
    data_type = 'uint16';
elseif info.data_type == 4
    data_type = 'single';
end
if info.byte_order == 0
    byte_order = 'ieee-le';
end

% Read the sample data - lines-5 : to avoid issues with the
% multibandread check 'file is too small to contain the
% specified data', it would require further analysis
info.lines = info.lines - 5;
data = multibandread(data_path, [info.lines, info.samples, info.bands],data_type,0,info.interleave,byte_order);

% Get dark and white references [data_dark_ave_rep,data_white_ave_rep]
[darkref_path, whiteref_path] = find_ref_path(dir_sample,cfg.darkref_tag,cfg.whiteref_tag);

% Load the dark and white files
data_white = multibandread(whiteref_path,[70, info.samples, info.bands],data_type,0,info.interleave,byte_order); 
data_dark = multibandread(darkref_path,[70, info.samples, info.bands],data_type,0,info.interleave,byte_order); 

% Resample if sample less than 5nm
crit = (max(info.Wavelength) - min(info.Wavelength))/length(info.Wavelength);
if crit < 3.75
    data = data(:,:,1:2:end);
    data_white = data_white(:,:,1:2:end);
    data_dark = data_dark(:,:,1:2:end);
    info.Wavelength = info.Wavelength(1:2:end);
end

% Average of the references
data_dark_ave = mean(data_dark,1);
data_white_ave = mean(data_white,1);

% Reflectance calculation
data = (data - repmat(data_dark_ave,info.lines,1,1))./(repmat(data_white_ave,info.lines,1,1) - repmat(data_dark_ave,info.lines,1,1) + eps);

% Set impossible values to NaN, correction is made further on to accelerate
% the process
data(data < 0) = 0;
data(data > 1) = 1;
end

function info = ReadEnviHdr(hdrfile)
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
end

function A = split(s,d)
%This function by Gerald Dalley (dalleyg@mit.edu), 2004
A = {};
while (~isempty(s))
    [t,s] = strtok(s,d);
    A = {A{:}, t};
end
end