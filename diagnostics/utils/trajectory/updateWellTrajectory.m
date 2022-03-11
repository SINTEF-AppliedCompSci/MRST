function [w_o, ws] = updateWellTrajectory(model, w, ws, traj)
% Undocumented utility function used for trajectory optimization 

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

opts = {'Name', w.name, 'Type', w.type, 'Val', w.val, 'Sign', w.sign, ...
        'lims', w.lims, 'Radius', w.r(1), 'compi', w.compi};
w_o  = addWellFromTrajectory([], model.G, model.rock, traj, ...
                    'exteriorFaceCorrection', true, ...
                    'errorOnEmptyCellSet', false, opts{:});
if isempty(w_o.cells)
    nc = 0;
else
    nc = numel(w_o.cells);
end
w_o.cell_origin = ones(nc, 1);
if isfield(w, 'posControl')
    w_o.posControl  = w.posControl;
end

if nc == 0
    warning('Empty set of well-cells for well: %s\n', w.name);
    w_o.status  = false;
    ws.status = false;
else
    if nargout > 1 
        if isempty(ws)
            warning('Can''t edit empty wellsol...')
        else
            ws.cstatus = true(nc,1);
            ws.cdp     = zeros(nc,1);
            flds = {'cqs', 'flux', 'ComponentTotalFlux'};
            for k = 1:numel(flds)
                ws.(flds{k}) = nan(nc, size(ws.(flds{k}),2));
            end
        end
    end
    
end
end


