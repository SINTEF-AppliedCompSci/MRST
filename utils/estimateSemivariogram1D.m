function [curve, sill, range] = ...
       estimateSemivariogram1D(x_cell, fx_cell, num_curvepts, max_reldist, rel_radius)
%
% Estimate semivariogram from one or more sets of sampled 1D data
% 
% SYNOPSIS:
%   function [curve, sill, range] = ...
%      estimateSemivariogram1D(x_cell, fx_cell, num_curvepts, ...
%                              max_reldist, rel_radius)
%
% DESCRIPTION:
%
% This function estimates the semivariogram associated with a set of 1D
% correllated noise, which is assumed to be of mean 0.  It does so by computing
% all the mutual distances and value differences between pairs of data
% points, sorting into bins and computing means.
% 
% If more than one set of sampled 1D data is available, the function can work
% with multiple sets, as datasets are provided as cell arrays (see parameter
% description below).
% 
% PARAMETERS:
%   x_cell      - cell array with one or more sets containing the x-positions of
%                 the 1D data
%   fx_cell     - cell array with one or more sets containing the sample
%                 values of the 1D data
%   num_curvepts - The number of points on the variogram curve to compute
%   max_reldist - The maximum extent of the semivariogram curve to compute
%                 (as relative propotion of maximum distance in line data).
%                 This extent may be shortened within the algorithm if a definite
%                 sill is encountered at a shorted distance.
%   rel_radius  - Width of averaging window used when computing each point on
%                 the semivariogram curve.  The with is given in units of bin
%                 width (where each bin has the width of the max length of
%                 the semivariogram curve, divided by number of curvepoints
%                 (num_curvepts). 
%
% RETURNS:
%   curve - Variogram curve presented as a vector of regularly sampeld values
%   sill  - maximum value encountered (where curve flattens)
%   range - the distance from origin to the distance value where the sill is
%           encountered 
%
% EXAMPLE:
%   For an example of use, refer to the script 'test_estimate_variogram'

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
   
   [c, d] = deal(cell(size(x_cell)));
   
   for i = 1:numel(x_cell)
      x = x_cell{i};
      fx = fx_cell{i};
      
      % compute matrix of mutual spatial distances
      d{i} = repmat(x(:), 1, numel(x));
      d{i} = abs(d{i}-d{i}');

      % compute matrix of squared random variable differences
      c{i} = repmat(fx(:), 1, numel(fx));
      c{i} = (c{i}-c{i}').^2;
   end
   
   % define samples for which we want to evaluate the variogram, and radius
   % of region for which to apply weighted averages when estimating the
   % variogram value
   [c_all, d_all] = deal([]);
   tovec = @(x) x(:);
   
   for i = 1:numel(x_cell)
      d_all = [d_all; tovec(d{i})];
      c_all = [c_all; tovec(c{i})];
   end
   
   max_len = max(d_all);
   vg_max = max_reldist * max_len; % maximal range of variogram
   aver_rad = rel_radius * vg_max/num_curvepts; % width of averaging window
   
   curve = nan(num_curvepts, 1);
   for i = 1:num_curvepts
      curve(i) = estimate_variogram_value(d_all(:), ...
                                          c_all(:), ...
                                          (i/num_curvepts) * vg_max, ...
                                          aver_rad);
   end
   
   % identify sill
   sill = max(curve);
   
   % identify range
   fac = 1; % to limit the extent where we have basically reached the
            % sill, but curve is very slowly rising
   tmp = find(curve >= fac * sill, 1, 'first');
   range = (tmp / numel(curve)) * vg_max;
   curve = [0; curve(1:tmp)];
   
end

% ----------------------------------------------------------------------------
function val = estimate_variogram_value(dists, vals, pos, rad)
   w = max(1 - (dists - pos)/rad, 0); % value of weighting function
   val = sum(vals .* w) ./ (2 * sum(w));
end
