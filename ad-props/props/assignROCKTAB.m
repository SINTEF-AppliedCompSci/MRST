function f = assignROCKTAB(f, rocktab, reg)
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

   if ~iscell(rocktab), rocktab = { rocktab }; end

   Tpv = cellfun(@(x) x(:, [1, 2]), rocktab, 'UniformOutput', false);
   Ttr = cellfun(@(x) x(:, [1, 3]), rocktab, 'UniformOutput', false);

   regmap = @(p, varargin) ...
      getRegMap(p, reg.ROCKNUM, reg.ROCKINX, varargin{:});

   f.pvMultR   = @(p, varargin) interpReg(Tpv, p, regmap(p, varargin{:}));
   f.tranMultR = @(p, varargin) interpReg(Ttr, p, regmap(p, varargin{:}));
end
