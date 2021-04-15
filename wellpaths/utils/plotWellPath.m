function h = plotWellPath(wellpaths, varargin)
%Plot a well path
%
% SYNOPSIS:
%   plotWellPath(wellpath)
%   h = plotWellPath(wellpath, 'color', 'r')
%
% DESCRIPTION:
%   Plots a given well path, using colors and showing control points along
%   the curve.
%
% REQUIRED PARAMETERS:
%
%   wellpath - Well path to be plotted. See makeSingleWellpath.
%
% OPTIONAL PARAMETERS:
%
%   interpType  - Interpolation type used. Same possible values as for
%               MATLAB builtin interp1. Default: Spline. 
%
%   LineWidth   - Line width of curve used to draw wellpath segments.
%
%   MarkerColor - Used to colorize the control points.
%
%   Color       - Color of the line segments themselves.
% 
%   Refinement  - The path is refined using interpolation to produce nice
%                 curves. Entering a positive number here will refine the
%                 curve by a number of points. Interpreted as a
%                 multiplicative factor.
%
% RETURNS:
%   h           - Two-column array of handles. The first column contains
%               handles for the line segments and the second for the
%               control point markers.
%
% SEE ALSO:
%   `makeSingleWellpath`, `findWellPathCells`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('interpType',   'spline', ...
             'LineWidth',    2, ...
             'Color',        [], ...
             'MarkerColor', 'r', ...
             'refinement',   25);
opt = merge_options(opt, varargin{:});

washeld = ishold();
hold on

if ~iscell(wellpaths)
    wellpaths = {wellpaths};
end
h = [];
for wp = 1:numel(wellpaths)
    wellpath = wellpaths{wp};
    
    npts = numel(wellpath.points);

    % Preallocate storage for handles to well trajectories & markers
    h_loc = nan(npts, 2);

    for i = 1:npts
        pts = wellpath.points{i};

        % Create a interpolated trajectory from the control points with
        % refinement to make a nice, smooth curve.
        [vals, v] = refineSpline(pts, opt.refinement, opt.interpType);

        % Mask away any inactive segments from plotting
        act = wellpath.active{i};
        segInd = floor(v);
        segInd = min(segInd, numel(act));
        isActive = act(segInd);
        vals(~isActive, :) = nan;
        
        if size(pts, 2) == 3
            plotter = @(x, varargin) plot3(x(:, 1), x(:, 2), x(:, 3), varargin{:});
        else
            plotter = @(x, varargin) plot(x(:, 1), x(:, 2), varargin{:});
        end
        % Make the curves for current segment
        if ~isempty(opt.Color)
            h_loc(i, 1) = plotter(vals,...
                'LineWidth', opt.LineWidth, 'Color', opt.Color);
        else
            h_loc(i, 1) = plotter(vals,...
                'LineWidth', opt.LineWidth);
        end
        % Plot control points (note that pts are the prescribed values, and
        % not interpolated).
        h_loc(i, 2) = plotter(pts, 'ro', 'MarkerFaceColor', opt.MarkerColor);

    end
    h = [h; h_loc];
end
% Return hold state to whatever it was before we started.
if ~washeld
    hold off
end
% Reservoir convention.
set(gca, 'ZDir', 'reverse')
end
