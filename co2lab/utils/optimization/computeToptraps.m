function dh = computeToptraps(trapfun, Gt, recenter)
%
% Compute subscale trapping potential for all cells in a top surface grid,
% based on a precomputed trap function.
%
% SYNOPSIS:
%   function dh = computeToptraps(trapfun, Gt, recenter)
%
% DESCRIPTION:
%
% PARAMETERS:
%   trapfun  - precomputed trap function (As computed by e.g. 'resTiltUtsira.m')
%   Gt       - top surface grid
%   recenter - whether or not to 'recenter', i.e. adjust angle so that maximum
%              subscale trapping is obtained when surface is horizontal
%
% RETURNS:
%   dh - vector with one entry per cell of 'Gt', giving the subscale trapping
%        potential for each cell.
%
% SEE ALSO:
%   `resTiltUtsira`

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

    TVol = @(theta_x, theta_y) interp2(trapfun.theta_x_vec, ...
                                       trapfun.theta_y_vec, ...
                                       trapfun.tv_per_area', ...
                                       theta_x, theta_y, 'spline');

    % Compute inclination of each gridcell in the 'Utsira' grid
    z_norm = bsxfun(@rdivide, Gt.cells.normals, Gt.cells.normals(:,3));
    tan_xy = -z_norm(:, 1:2);
    
    % If recentering is on, compute the adjustment needed to have maximum
    % trapping when surface is horizontal
    if recenter
        samples = 80;
        xsamples = linspace(trapfun.theta_x_vec(1), trapfun.theta_x_vec(end), samples);
        ysamples = linspace(trapfun.theta_y_vec(1), trapfun.theta_y_vec(end), samples);
        fsampled = computeSamples(TVol, xsamples, ysamples);

        [max_vol, max_ix] = max(fsampled(:));%#ok
        [max_x, max_y] = ind2sub(size(fsampled), max_ix);
        [adj_x, adj_y] = deal(xsamples(max_x), ysamples(max_y));
    else
        [adj_x, adj_y] = deal(0,0);
    end
    
    % We consider gridcells with steeper inclination than the precomputed
    % table to have practically zero subscale trapping.
    min_tilt_x = trapfun.theta_x_vec(1)   + adj_x;
    max_tilt_x = trapfun.theta_x_vec(end) + adj_x;
    min_tilt_y = trapfun.theta_y_vec(1)   + adj_y; 
    max_tilt_y = trapfun.theta_y_vec(end) + adj_y;
    
    % We only consider a subset of the gridcells to have nonzero subscale
    % trapping: 
    % (1) Those whose inclination doesn't surpass the precomputed angular domain
    % (2) Those who are not in a coarse trap
    ind = ( (tan_xy(:,1) < max_tilt_x) & (tan_xy(:,1) > min_tilt_x) & ...
            (tan_xy(:,2) < max_tilt_y) & (tan_xy(:,2) > min_tilt_y) );
    
    % computing estimated subscale trapping for each relevant cell, based on
    % local inclination 
    dh = zeros(Gt.cells.num, 1);
    dh(ind) = TVol(tan_xy(ind,1) - adj_x, tan_xy(ind,2) - adj_y);

end

% ----------------------------------------------------------------------------

function fun_sampled = computeSamples(fun, samples_x, samples_y)
% Sample 2D function 'fun' at a cartesian grid of sample points, defined by
% 'samples_x' and 'samples_y'
    [x, y] = meshgrid(samples_x, samples_y);
    fun_sampled = fun(x, y);
end

