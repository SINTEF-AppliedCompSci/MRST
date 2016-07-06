function wellSol = assignWellValuesFromControl(model, wellSol, W, wi, oi, gi)
% Utility function to assign well values that are the resut of the controls
% being enabled.
    for w = 1:numel(wellSol)
        ws = wellSol(w);
        tp = ws.type;
        v  = ws.val;
        switch tp
            case 'bhp'
                ws.bhp = v;
            case 'rate'
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
                    ws.qWPoly = ws.qWs*W(w).poly;
                end
                if isprop(model, 'surfactant') && model.surfactant
                    ws.qWSft = ws.qWs*W(w).surfact;
                end
            case 'orat'
                ws.qOs = v;
            case 'wrat'
                ws.qWs = v;
            case 'grat'
                ws.qGs = v;
            case 'lrat'
                % Do nothing
            otherwise
                error('Unknown well control mode');
        end
        wellSol(w) = ws;
    end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
