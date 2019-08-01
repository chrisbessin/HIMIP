function [sampleFile, whiterefFile, darkrefFile] = search_data_files (data_path)

% Search corresponding file names for hyperspectral data from a data path

% Author: Thanh Bui (thanh.bui@erametgroup.com)

file_name = dir(data_path);
file_name_raw = {};
k = 1;
for i = 3: length(file_name)
    if((length(file_name(i).name) > 2) && strcmp(file_name(i).name(end-2:end), 'raw'))
        file_name_raw{k} = file_name(i).name;
        k = k + 1;
    end
end
for i = 1: length(file_name_raw)
    if(strcmp(file_name_raw{i}(1:4), 'DARK'))
        darkrefFile  = fullfile(data_path, file_name_raw{i});
    elseif(strcmp(file_name_raw{i}(1:4), 'WHIT'))
        whiterefFile = fullfile(data_path, file_name_raw{i});
    elseif (~strcmp(file_name_raw{i}(end-7:end-4),'refl') && ~strcmp(file_name_raw{i}(end-9:end-4),'reflCR'))
        sampleFile = fullfile(data_path, file_name_raw{i});
    end
    
end
fprintf('===========Reading and preprocessing the data===================\n')
fprintf('Sample file: %s \n', sampleFile)
fprintf('White ref: %s \n', whiterefFile)
fprintf('Dark ref: %s \n', darkrefFile)