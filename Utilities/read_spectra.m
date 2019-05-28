function [data, legendTxt] = read_spectra(filename)
% Read spectra from a .txt file

% Author: Thanh Bui (thanh.bui@erametgroup.com)

fid = fopen(filename, 'r');
skipLines = 0;
i = 1;
legendTxt = cell(1,1);
while(1)
    tline = fgets(fid);   
    if strcmp(tline(1:2), '  ')
        break
    end
    skipLines = skipLines + 1;
    if skipLines >= 3
        if length(tline) > 25
            legendTxt{i} = tline(11:25);
        else
            legendTxt{i} = tline(11:end);
        end
        i = i + 1;
    end    
    
end

fid = fopen(filename, 'r');
for i = 1: skipLines
    tline = fgets(fid);
end
nCol = skipLines - 1;
formatSpec = repmat('%f',[1, nCol]); % Generate format spec automatically
sizeA = [nCol Inf];
data = fscanf(fid,formatSpec,sizeA);
data = data';
fclose(fid);