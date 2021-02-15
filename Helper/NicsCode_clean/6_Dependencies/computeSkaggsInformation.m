function [pfInfo, rMean] = computeSkaggsInformation(pX,r)
% pfInfo = computeSkaggsInformation(pX,r)
%
%   computes information about position conveyed by the firing rate of a neuron
%
%
% pX - X-length vector, the probability of being at each position
%          pX is normalized to have sum = 1
% r  - firing rate at each location
%
%
% from Skaggs et al 1993


% ensure same size
if numel(pX) ~= numel(r)
    error('pX and r are different sizes')
end

% ensure all non-negative
if any(pX<0), error('some elements of pX are less than 0'), end
if any(r<0), error('some elements of r are less than 0'), end

% ensure both columns
pX = pX(:);
r = r(:);

% normalize probability
pX = pX/sum(pX);

% compute overall rate
rMean = sum(pX.*r);

% for information computation, ignore points with firing rate zero
% otherwise the results is a NaN
rZero = r==0;
pX = pX(~rZero);
r = r(~rZero);

% compute information
pfInfo = sum(pX.*r.*log2(r./rMean));