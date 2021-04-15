function f = assignPVDG(f, pvdg, reg)
    [f.bG, f.muG] = getFunctions(pvdg, reg);
end

function [bG, muG] = getFunctions(PVDG, reg)
    [bG, muG] = deal(cell(1, reg.pvt));
    
    for i = 1:reg.pvt
        pvdg = PVDG{i};
        
        bg = 1./pvdg(:, 2);
        mug = pvdg(:, 3);
        % Extend table at lower end
        pG = [0; pvdg(:, 1)];
        bg = bg([1, 1:end]);
        mug = mug([1, 1:end]);
        bG{i}  = @(pg) reg.interp1d(pG, bg, pg);
        muG{i} = @(pg) reg.interp1d(pG, mug, pg);
    end
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
