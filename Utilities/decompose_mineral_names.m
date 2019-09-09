function [mineralDict, mineralList] = decompose_mineral_names(min_names)
% Extract mineral names

% Author: Thanh Bui (thanh.bui@erametgroup.com)

mineralList = min_names';
% Remove file extension in the file name
for i = 1: length(mineralList)
    mineralList{i} = mineralList{i}(1:end-4);
end
keySet = cell(1);
valueSet = cell(1);

j = 1;
for i = 1: length(mineralList)
    % Extract a mineral name
    m = mineralList{i};
    underInd = strfind(m, '_');
    if isempty(underInd)
        mineralName = m;
    else
        mineralName = m(1:underInd(1)-1);
    end
    if (i == 1) % Add the first mineral to cells
        keySet{j} = mineralName;
        valueSet{j} = i;
    else
        if strcmp(mineralName, keySet{j})   % if the mineral already existed, append the value
            valueSet{j} = [valueSet{j}, i];
        else    % if the mineral does not exist, add it to key and values cells
            j = j + 1;
            keySet{j} = mineralName;
            valueSet{j} = i;           
        end 
    end  
end
% Create a dictionary
mineralDict = containers.Map(keySet, valueSet);