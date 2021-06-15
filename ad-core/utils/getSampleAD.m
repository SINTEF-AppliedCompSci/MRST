function [x, isAD] = getSampleAD(varargin)
% Utility for getting a AD value if it exists from a list of possible
% AD-values

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

    x = varargin{1};
    isAD = false;
    for i = 1:numel(varargin)
        if isa(varargin{i}, 'ADI')
            x = varargin{i};
            isAD = true;
            return
        end
    end
end
