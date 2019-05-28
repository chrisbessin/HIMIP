function precision = data2precision(data_type)

switch data_type
    case 12
        precision = 'uint16';
    case 13
        precision = 'uinit32';
    case 4
        precision = 'single';
    case 5
        precision = 'double';
    otherwise
        fprintf('Do not know the precision, please check \n')
end