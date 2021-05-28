function [maps, u] = setupSimulationControlMappings(schedule, bounds, controllableWells)
% setup a map-structure to scaled control-vector from schedule (see OptimizationProblem)

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

[nw, ns] = deal(numel(schedule.control(1).W), numel(schedule.control));
if nargin < 3 || isempty(controllableWells)
    % use first control
    controllableWells = true(nw, 1);
end

if ~any(controllableWells)
    [maps.wellNo, maps.type, u] = deal([]);
    return;
end

isControl = getControls(schedule.control(1).W, controllableWells);
contrNms = fieldnames(isControl);
nc = ns*sum(structfun(@nnz, isControl));

nc_bnds = ns*sum(structfun(@(x)nnz(isfinite(x(:,1))), bounds));

if nc_bnds > nc
    warning('Bounds were given for %d controls/limits of non-finite value in schedule, these will be ignored', nc_bnds-nc);
elseif nc_bnds < nc
    warning('Bounds were missing for % d controls/limits of finite value in schedule, this will probably result in an error', nc-nc_bnds);
end

maps = struct('type', {repmat({''}, [nc,1])}, 'wellNo', nan(nc,1), ...
              'stepNo', rldecode((1:ns)', nc/ns), 'isTarget', false(nc,1), ...
              'bounds', nan(nc, 2));
u = nan(nc,1);

scale = @(v, bnds)(v-bnds(1))/(bnds(2)-bnds(1));
[ix, offset] = deal(0, nc/ns);
for k = 1:numel(contrNms)
    cnm = contrNms{k};
    for wno = 1:nw
        if isControl.(cnm)(wno)
            ix = ix +1;
            for step = 1:ns
                w = schedule.control(step).W(wno);
                ixs = ix + (step-1)*offset;
                maps.type{ixs}   = cnm;
                maps.wellNo(ixs) = wno;
                if ~isfield(bounds, cnm)
                    if ~strcmp(cnm, 'thp')
                        warning('Schedule contains finite limits on %s, but no bounds were given');
                    end
                else
                    bnds = bounds.(cnm)(wno,:);
                    assert(all(isfinite(bnds)), 'No bounds given for ...')
                    maps.bounds(ixs,:) = bnds;
                end
                if strcmp(w.type, cnm)
                    u(ixs) = scale(w.val, bnds);
                    maps.isTarget(ixs) = true;
                else
                    u(ixs) = scale(w.lims.(cnm), bnds);
                end
            end
        end
    end
end

assert(~any(isnan(u)), 'Unable to produce controls, probably missing some bounds');
if any(u < -sqrt(eps) | u > 1 + sqrt(eps))
    warning('Initial controls are not within bounds ...')
end
end

function contr = getControls(W, cix)
nw   = numel(W);
flds = {'bhp', 'wrat', 'orat', 'grat', 'lrat', 'rate'};
inp  = [flds; repmat({false(nw,1)}, [1, numel(flds)]) ];
contr = struct(inp{:});
for k = 1:nw
    if cix(k)
        contr.(W(k).type)(k) = true;
        if isfield(W(k), 'lims') && ~isempty(W(k).lims)
            nms = fieldnames(W(k).lims);
            for k1 = 1:numel(nms)
                if isfinite(W(k).lims.(nms{k1})) && ~strcmp(nms{k1}, 'thp')
                    contr.(nms{k1})(k) = true;
                end
            end
        end
    end
end
end
