function f = sre(X, X_hat)
% Compute signal to reconstruction error (RSNR); it gives more information regarding
% the power of the signal in relation with the power of the error

% Author: Thanh Bui (thanh.bui@erametgroup.com)
temp = norm(X)/norm(X-X_hat);
f = 10*log10(temp);