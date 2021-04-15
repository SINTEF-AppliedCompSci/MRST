function colorizedHistogram(data, n_sample)
%Undocumented Utility Function

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

    if nargin < 2
        n_sample = 10;
    end
    assert(size(data, 1) == 1 || size(data, 2) == 1, ...
        'Matrix input not supported. Please provide a vector as input.')
    data = double(data);
    data = data(all(isfinite(data), 2), :);
    % Find data
    [n, X] = hist(data, n_sample);

    h = X(2) - X(1);
    for i = 1:n_sample
        plotRectangle(X(i) - h, n(i), h, X(i));
    end
end


function h = plotRectangle(x, height, width, value, varargin)
    h = patch([x; x+width; x+width; x], [0; 0; height; height], value, varargin{:});
end
