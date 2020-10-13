function state = dgLimiter2(disc, state, bad)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    nDofMax = disc.basis.nDof;
    G       = disc.G;
    
    varNames = {'s', 'rs', 'rv', 'c'};
    for vNo = 1:numel(varNames)
        v    = varNames{vNo};
        if isfield(state, v)
            vdof = [v, 'dof'];
            isActive = size(state.(v),1) == disc.G.cells.num;
            if isActive
                ix = disc.getDofIx(state, 1, bad);
                state.(vdof)(ix,:) = state.(v)(bad,:);
            end
            if disc.degree > 0 && isActive
                ix = disc.getDofIx(state, 2:nDofMax, bad);
                state.(vdof)(ix,:) = [];
            end
        end
    end
    
    if disc.degree > 0
        state.degree(bad) = 0;
    end
    state = disc.updateDofPos(state);
    for vNo = 1:numel(varNames)
        v    = varNames{vNo};
        if isfield(state, v) && size(state.(v),1) == disc.G.cells.num
            vdof = [v, 'dof'];
            state.(v) = zeros(G.cells.num, size(state.(v),2));
            for cNo = 1:size(state.(v),2)
                state.(v)(:,cNo) = disc.getCellMean(state, state.(vdof)(:,cNo));
            end
        end
    end
end
