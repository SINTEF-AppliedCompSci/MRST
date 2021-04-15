function [L, x, y, Z_L, Z_V, p, T] = getFlashTable(p_range, T_range, z, EOS)
% Utility for getting meshgrid plots of flash results (e.g. for phase diagram)

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

    [p, T] = meshgrid(p_range, T_range);
    np = numel(p_range);
    nT = numel(T_range);
    
    n = numel(p);
    ncomp = numel(z);
    Z = cell(1, ncomp);
    for i = 1:ncomp
        if iscell(z)
            zi = z{i};
        else
            zi = z(i);
        end
        Z{i} = repmat(zi, n, 1);
    end
    
    [L, x, y, Z_L, Z_V] = standaloneFlash(p(:), T(:), Z, EOS);
    
    fix = @(x) reshape(x, nT, np)';
    L = fix(L);
    Z_L = fix(Z_L);
    Z_V = fix(Z_V);
    for i = 1:numel(x)
        x{i} = fix(x{i});
        y{i} = fix(y{i});
    end
    T = fix(T);
    p = fix(p);
end
