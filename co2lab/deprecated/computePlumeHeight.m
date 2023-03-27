function [h, h_max] = computePlumeHeight(Gt, state, sw, sr)
% computing height corresponding to the current saturation state

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

    if isfield(state, 'sGmax')
        % We operate with dissolution, implying a gradually receding 'max' saturation
        smax = state.sGmax;
    else
        % There is no dissolution, so the historical maximum is the current
        % maximum
        smax = state.smax(:,2);
    end

    [h, h_max] = upscaledSat2height(state.s(:,2), smax, Gt, 'resSat', [sw sr]);
end
