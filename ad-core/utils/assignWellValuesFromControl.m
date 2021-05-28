function wellSol = assignWellValuesFromControl(model, wellSol, W, wi, oi, gi)
% Assign wellSol values when values are set as controls
%
% SYNOPSIS:
%   wellSol = assignWellValuesFromControl(model, wellSol, W, wi, oi, gi)
%
% DESCRIPTION:
%   Well rates and pressures can be both controls and solution variables,
%   depending on the problem. For a subset of possible well controls, this
%   function explicitly assigns the values to the wellSol.
%
% REQUIRED PARAMETERS:
%   model   - ReservoirModel-derived subclass.
%
%   wellSol - wellSol to be updated.
%
%   W       - Well struct used to create wellSol.
%
%   wi, oi, gi - Indices for water, oil and gas respectively in the .compi
%                field of the well.
%
% RETURNS:
%   wellSol - Updated wellSol where fields corresponding to assigned
%             controls have been modified.
%
% SEE ALSO:
%   WellModel

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
    for w = 1:numel(wellSol)
        ws = wellSol(w);
        if ~ws.status
            continue
        end
        tp = ws.type;
        v  = ws.val;
        switch tp
            case 'bhp'
                ws.bhp = v;
            case 'rate'
                if ws.sign < 1
                    continue;
                end
                if model.water
                    ws.qWs = v*W(w).compi(wi);
                end
                if model.oil
                    ws.qOs = v*W(w).compi(oi);
                end
                if model.gas
                    ws.qGs = v*W(w).compi(gi);
                end
                if isprop(model, 'polymer') && model.polymer
                    ws.cWPoly = W(w).cp;
                end
                if isprop(model, 'surfactant') && model.surfactant
                    ws.qWSft = ws.qWs*W(w).cs;
                end
            case 'orat'
                ws.qOs = v;
            case 'wrat'
                ws.qWs = v;
            case 'grat'
                ws.qGs = v;
            case {'lrat', 'thp', 'resv'}
                % Do nothing
            otherwise
                error('Unknown well control mode %s', tp);
        end
        wellSol(w) = ws;
    end
end
