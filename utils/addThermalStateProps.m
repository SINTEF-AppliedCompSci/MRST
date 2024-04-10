function state = addThermalStateProps(state, varargin)
%Add thermal properties to an existing state structure

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('T', convertFromCelcius(20));
    opt = merge_options(opt, varargin{:});
    nc = length(state.pressure);
    if size(opt.T, 1) < nc
        opt.T = repmat(opt.T, nc, 1);
    end
    assert(numel(opt.T) == nc);
    state.T = opt.T;
    
end