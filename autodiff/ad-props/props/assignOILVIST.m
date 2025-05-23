function f = assignOILVIST(f, table, reg)
   % Compute tables (static data)
   tab = cellfun(@(x)x(:,[1,2]), table, 'UniformOutput', false);
   f.muOT =@(T,varargin) func(T,tab,reg,varargin{:});
end
function v = func(x, tab, reg, varargin)
inx = getRegMap(x, reg.SATNUM, reg.SATKINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), tab, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, x, inx);
end

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
