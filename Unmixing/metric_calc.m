function [sam, mae] = metric_calc(x, y)

% Calculate spectral angle mapper and mean absolute error between hyperspetral images 
% Author: Thanh Bui (thanh.bui@erametgroup.com)

%x = randn(10, 4);
%y = randn(10, 4);
n = size(x);

mae = mean(abs(x-y),1);
%num = diag(x'*y);
num = sum(x.*y,1);
norm_x = zeros(1,n(2));
norm_y = zeros(1,n(2));
for i = 1: n(2)
    norm_x(i) = norm(x(:,i));
    norm_y(i) = norm(y(:,i));
end
sam = acos(num./(norm_x.*norm_y));
