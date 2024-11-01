function [Z] = SpatialSignal(X, Y, n, N_std)

if nargin == 3
    N_std = 0;
end
name = sprintf('NewWave_%d',n);
load(name);
Z = zeros(size(X));

% Generation of the large scale Gaussians
for k = 1: N
    Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)).^2 + (Y-P(k,2)).^2)/2/Sigma^2);
end

% Genmeration of small scale Gaussians
for k = 1 : Ns
    Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
end

% Snapshot normalization
Z = Peak / 2 * Z / MaxMin + 50;

% Inclusion of additive white noise
Z = Z + N_std * randn(size(Z)); 

