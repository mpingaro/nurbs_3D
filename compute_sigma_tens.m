function SG = compute_sigma_tens(D)
% Compute the tensor S = second Piola Kirchhoff


%sigma = zeros(3);
sigma = ones(3);
SG = blkdiag(sigma,sigma,sigma);

return