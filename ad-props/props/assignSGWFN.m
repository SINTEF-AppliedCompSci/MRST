function f = assignSGWFN(f, sgfn, reg)
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

   cfun = @(f) cellfun(f, sgfn, 'UniformOutput', false);

   % Compute tables (static data)
   Tkrg  = extendTab( cfun(@(x) x(:, [1, 2])) );
   Tkrw  = extendTab( cfun(@(x) flipud([1 - x(:, 1), x(:, 3)])));
   Tpcwg = extendTab( cfun(@(x) x(:, [1, 4])) );

   % Region mapping
   regmap = @(sg, varargin) ...
      getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});

   % Region interpolator
   ireg = @(T, sg, varargin) interpReg(T, sg, regmap(sg, varargin{:}));

   f.krG  = @(sg, varargin) ireg(Tkrg , sg, varargin{:});
   f.krW  = @(sw, varargin) ireg(Tkrw, sw, varargin{:});
   f.pcWG = @(sg, varargin) ireg(Tpcwg, sg, varargin{:});
end
