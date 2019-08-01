function y = normalize(x)

% Normalize the input matrix to range of [0 1]
% Author: Thanh Bui (thanh.bui@erametgroup.com)

n = size(x);
minX = repmat(min(x), [n(1), 1]);
maxX = repmat(max(x), [n(1), 1]);
y = (x - minX)./(maxX - minX);