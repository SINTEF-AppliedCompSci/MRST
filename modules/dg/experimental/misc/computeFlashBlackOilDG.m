function state = computeFlashBlackOilDG(state, state0, model, status)
%Undocumented Utility Function

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

    statePrev = state;
    % Flash
    state = computeFlashBlackOil(state, state0, model, status);
    names = {'s', 'rs', 'rv'};
    for vNo = 1:numel(names)
        v = names{vNo};
        if isfield(state, v) && size(state.(v),1) == model.G.cells.num
            state = scaleDof(model.disc, state, statePrev, v);
        end
    end

end

function state = scaleDof(disc, state, statePrev, v)

    vdof = [v, 'dof'];
    % Add to change to constant part
    d  = state.(v) - statePrev.(v);
    ix = disc.getDofIx(state, 1, Inf);
    state.(vdof)(ix,:) = state.(vdof)(ix,:) + d;
    % Adjust mean
    for cNo = 1:size(state.(vdof),2)
        val = disc.getCellMean(state, state.(vdof)(:,cNo));
        f   = state.(v)(:,cNo)./val;
        f(~isfinite(f)) = 1;
        state.(vdof)(:,cNo) = state.(vdof)(:,cNo).*rldecode(f, state.nDof, 1);
    end
    
    if 1
        for cNo = 1:size(state.(vdof),2)
            val = disc.getCellMean(state, state.(vdof)(:,cNo));
            d   = val - state.(v)(:,cNo);
        end
    end

end
