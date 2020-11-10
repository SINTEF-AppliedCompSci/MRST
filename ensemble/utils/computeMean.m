function mean = computeMean(meanA, meanB, nA, nB)
% Compute mean of two existing means that are based on different
% sampling sizes.
% 
% DESCRIPTION:
%   Computes the weighted mean between the two means.
%
% SYNOPSIS:
%   mean = computeMean(meanA, meanB, nA, nB)
%
% PARAMETERS:
%   meanA - mean of sample set A
%   meanB - mean of sample set B
%   nA    - sample size of set A
%   nB    - sample size of set B
%
% RETURNS:
%   mean - the combined mean

%{
#COPYRIGHT#
%}
    mean = (nA.*meanA + nB.*meanB)./(nA + nB);
end