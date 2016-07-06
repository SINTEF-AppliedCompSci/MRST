function f = assignPVDO(f, pvdo, reg)
   cfun    = @(f) cellfun(f, pvdo, 'UniformOutput', false);

   % Compute tables (static data)
   TBO     = cfun(@(x) x(:, [1, 2]));
   TbO     = cfun(@(x) [x(:,1), 1 ./ x(:,2)]);
   TmuO    = cfun(@(x) x(:, [1, 3]));
   TBOxmuO = cfun(@(x) [x(:,1), prod(x(:, [2, 3]), 2)]);

   % Region mapping
   regmap = @(po, varargin) ...
      getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});

   % Region interpolator
   ireg = @(T, po, varargin) interpReg(T, po, regmap(po, varargin{:}));

   f.BO     = @(po, varargin) ireg(TBO,     po, varargin{:});
   f.bO     = @(po, varargin) ireg(TbO,     po, varargin{:});
   f.muO    = @(po, varargin) ireg(TmuO,    po, varargin{:});
   f.BoxmuO = @(po, varargin) ireg(TBOxmuO, po, varargin{:});
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
