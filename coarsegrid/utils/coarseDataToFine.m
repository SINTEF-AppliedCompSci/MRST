function data = coarseDataToFine(CG, data)
% Convert coarse grid dataset into fine grid representation

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

    if isempty(data)
        return
    end
    if isnumeric(data)
        data = getData(CG, data);
        return
    end
    
    ic = iscell(data);
    
    for i = 1:numel(data)
        if ic
            d = data{i};
        else
            d = data(i);
        end
        
        if isnumeric(d)
            d = getData(CG, d);
        else
            fn = fieldnames(d);
            for j = 1:numel(fn)
                f = fn{j};
                d.(f) = getData(CG, d.(f));
            end
        end
        if ic
            data{i} = d;
        else
            data(i) = d;
        end
    end
end

function data = getData(CG, data)
    sz = size(data);
    if sz(1) == CG.cells.num
        data = data(CG.partition, :);
    elseif all(sz == [1, CG.cells.num])
        data = data(CG.partition);
    end
end
