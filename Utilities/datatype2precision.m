function precision = datatype2precision(data_type)

% Convert data type (number) to a corresponding precision

% Author: Thanh Bui (thanh.bui@erametgroup.com)

switch data_type
    case 12
        precision = 'uint16';
    case 13
        precision = 'uint32';
    case 4
        precision = 'single';
    case 5
        precision = 'double';
    otherwise
        fprintf('Unknown data type, please check \n')
end