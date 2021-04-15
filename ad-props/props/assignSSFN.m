function f = assignSSFN(f, ssfn, reg)
   cfun = @(f) cellfun(f, ssfn, 'UniformOutput', false);

   % Compute tables (static data)
   Tkrfg  = extendTab( cfun(@(x) x(:, [1, 2])) );
   Tkrfs = extendTab( cfun(@(x) x(:, [1, 3])) );

   % Region mapping
   regmap = @(sg, varargin) ...
      getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});

   % Region interpolator
   ireg = @(T, sg, varargin) interpReg(T, sg, regmap(sg, varargin{:}));

   f.krFG  = @(sg, varargin) ireg(Tkrfg , sg, varargin{:});
   f.krFS = @(sg, varargin) ireg(Tkrfs, sg, varargin{:});
end

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