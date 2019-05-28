function [xreduced, ureduced, xapp] = BandReduction(X, VarThresh)
% Dimensionality reduction using PCA, and 
% Input:
%     X (n x m x number of bands): the input reflectance
%     VarThresh: the percentage of variance to be retained
% Output:
%     xreduced: (n x m x k): dimensionality
%     ureduced
%     xapp: reconstructed data
    
n= size(X);
X_2D = reshape(X, [n(1)*n(2),n(3)]);
sigma = 1/n(1) * (X_2D') *X_2D; % Compute covariance
[u,s,~] = svd(sigma);
if nargin == 1
    VarThresh =  0.999; % 99.9% of variance is retained
end
sDiag = diag(s);
for k = 1:n(1)
    if(sum(sDiag(1:k))/sum(sDiag) >= VarThresh ) 
        break
    end
end
ureduced = u(:,1:k);
xreduced_temp = X_2D*ureduced;
xreduced = reshape(xreduced_temp, [n(1), n(2), k]);
xapp_temp = (X_2D*ureduced)*ureduced';  % Reconstructed X
xapp = reshape(xapp_temp, [n(1), n(2), n(3)]);
    

