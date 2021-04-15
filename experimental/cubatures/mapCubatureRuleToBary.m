function [lambda, w] = mapCubatureRuleToBary(rule, type)
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

    if strcmpi(type, 'triangle')
        t = '%f%f';
    elseif strcmpi(type, 'tetrahedron')
        t = '%f%f%f';
    end

    fid = fopen([rule, '_x.txt']);
    x = textscan(fid, t, 'HeaderLines', 0, 'CollectOutput', 1);
    x = x{:};
    fclose(fid);
    
    fid = fopen([rule, '_w.txt']);
    w = textscan(fid, '%f', 'HeaderLines', 0, 'CollectOutput', 1);
    w = w{:};
    fclose(fid);
    
    fid = fopen([rule, '_r.txt']);
    v = textscan(fid, t, 'HeaderLines', 0, 'CollectOutput', 1);
    v = v{:};
    fclose(fid);
    
    n = numel(w);
    
    x = [x, ones(size(x,1),1)];
    R = [v, ones(size(v,1),1)];
    lambda = x/R;

end
