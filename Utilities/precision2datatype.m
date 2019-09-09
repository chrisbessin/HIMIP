function data_type = precision2datatype(precision)
% Convert precision to data type complied with ENVI data format

% Author: Thanh Bui (thanh.bui@erametgroup.com)

switch precision
    case 'uint16'
        data_type = 12;
    case 'uint32'
        data_type = 13;
    case 'single'
        data_type = 4;
    case 'double'
        data_type = 5;
    otherwise
        fprintf('Unknown precision, please check \n')
end
        