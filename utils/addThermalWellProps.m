function W = addThermalWellProps(W, varargin)
%Add thermal properties to an existing well structure.
% 
% SYNOPSIS: 
%   W = addThermalWellProps(W, 'pn1', pv1)
% 
% PARAMETERS: 
%   W - Well structure created with e.g. addWell.
% 
%   T - Injection temperature at the well. 
% 
% RETURNS:
%   W - valid well structure with thermal properties.
% 
% SEE ALSO:
% 'addWell'.

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

    opt = struct('T', 400*Kelvin);
    opt = merge_options(opt, varargin{:});
    
    fNames = fieldnames(opt);
    nW = numel(W);
    for fNo = 1:numel(fNames)
        fn = fNames{fNo};
        v = opt.(fn);
        if numel(v) < nW
            v = repmat(v, nW, 1);
        end
        for wNo = 1:nW
            W(wNo).(fn) = v(wNo);
        end
    end
    
end