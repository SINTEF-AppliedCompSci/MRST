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
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    mean = (nA.*meanA + nB.*meanB)./(nA + nB);
end