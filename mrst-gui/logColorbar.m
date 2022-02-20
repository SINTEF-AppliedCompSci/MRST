function h = logColorbar(varargin)
%Create colorbar for log dataset
%
% SYNOPSIS:
%   h = logColorbar();
%
% NOTE:
%   This function is identical to Matlab builtin colorbar, but with tick
%   labels indicating that the data is log10-transformed.
%
% SEE ALSO:
%   `colorbar`

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


    h = colorbar(varargin{:});
    
    ticks = get(h, 'YTick');
    newticks = arrayfun(@(x) ['1e', num2str(x)], ticks, 'UniformOutput', false);
    set(h, 'YTickLabel', newticks)
end
