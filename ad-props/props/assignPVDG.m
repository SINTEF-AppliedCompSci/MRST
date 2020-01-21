function f = assignPVDG(f, pvdg, reg)
    [f.bG, f.muG] = getFunctions(f, pvdg, reg);
end

function [bG, muG] = getFunctions(f, PVDG, reg)
    [bG, muG] = deal(cell(1, reg.pvt));
    
    for i = 1:reg.pvt
        pvdg = PVDG{i};
        
        pG = pvdg(:, 1);
        BG = pvdg(:, 2);
        mug = pvdg(:, 3);
        bG{i}  = @(pg) reg.interp1d(pG, 1./BG, pg);
        muG{i} = @(pg) reg.interp1d(pG, mug, pg);
    end
end



function f = assignPVDGOld(f, pvdg, reg)
   cfun = @(f) cellfun(f, pvdg, 'UniformOutput', false);

   % Compute tables (static data)
   TbG  = cfun(@(x) [x(:,1), 1 ./ x(:,2)]);
   TmuG = cfun(@(x) x(:, [1, 3]));

   % Region mapping
   regmap = @(pg, varargin) ...
      getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});

   % Region interpolator
   ireg = @(T, pg, varargin) interpReg(T, pg, regmap(pg, varargin{:}));

   f.bG  = @(pg, varargin) ireg(TbG , pg, varargin{:});
   f.muG = @(pg, varargin) ireg(TmuG, pg, varargin{:});
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
