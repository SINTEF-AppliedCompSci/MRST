function [s_min, s_max] = getMinMaxPhaseSaturations(model, satnum, cells)
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

    if nargin < 3
       cells = [];
    end
    if nargin < 2
        satnum = 1;
    end
    nc = numel(satnum);
    hasEndpoints = isfield(model.rock, 'krscale');
    if hasEndpoints
        % Endpoint scaling
        assert(~isempty(cells));
        pts = model.rock.krscale.drainage;
        [s_min, s_max] = getBounds(model, pts, cells);
    else
        if isfield(model.fluid, 'krPts')
            pts = model.fluid.krPts;
            [s_min, s_max] = getBounds(model, pts, satnum);
        else
            model = model.validateModel();
            sat_reg = model.FlowPropertyFunctions.RelativePermeability.regions(cells(1));
            [s_min, s_max] = getMinMaxPhaseSaturationsFromRelPerm(model, 1e-6, sat_reg);
        end
        s_min = repmat(s_min, nc, 1);
        s_max = repmat(s_max, nc, 1);
    end
    f = model.fluid;
    wix = model.getPhaseIndex('W');
    oix = model.getPhaseIndex('O');
    gix = model.getPhaseIndex('G');
    if model.gas
        % Minimum gas should always be zero
        s_min(:, gix) = 0;
    end
    
    if model.water
        if isfield(f, 'sWcon') && ~hasEndpoints
            % Minimum water = connate water
            if numel(f.sWcon) == 1
                s_min(:, wix) = f.sWcon;
            else
                s_min(:, wix) = f.sWcon(satnum);
            end
        end
        if isfield(model.rock, 'sw')
            assert(~isempty(cells));
            assert(hasEndpoints, 'Requires endpoint scaling');
            s_min(:, wix) = model.rock.sw(cells);
        end
        
        % Water can always be one
        s_max(:, wix) = 1;
        if model.gas
            % Account for minimum water saturation
            s_max(:, gix) = min(s_max(:, gix), 1 - s_min(:, wix));
        end
    end
end

function [s_min, s_max] = getBounds(model, pts, subs)
    phases = model.getPhaseNames();
    nph = numel(phases);
    nc = numel(subs);
    [s_min, s_max] = deal(ones(nc, nph));
    min_sub = 2;
    max_sub = 3;
    for i = 1:nph
        ph = lower(phases(i));
        if strcmp(ph, 'o') && ~isfield(pts, 'o')
            if model.water && model.gas
                s_min(:, i) = min(pts.ow(subs, min_sub), pts.og(subs, min_sub));
                s_max(:, i) = max(pts.ow(subs, max_sub), pts.og(subs, max_sub));
            elseif model.water
                s_min(:, i) = pts.ow(subs, min_sub);
                s_max(:, i) = pts.ow(subs, max_sub);
            elseif model.gas
                s_min(:, i) = pts.og(subs, min_sub);
                s_max(:, i) = pts.og(subs, max_sub);
            end
        else
            s_min(:, i) = pts.(ph)(subs, min_sub);
            s_max(:, i) = pts.(ph)(subs, max_sub);
        end
    end
end
