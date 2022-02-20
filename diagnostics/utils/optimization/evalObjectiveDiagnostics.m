function [val, der, W, state, D] = evalObjectiveDiagnostics(u, obj, state, system, G, fluid, pv, T, W, scaling, varargin)
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

opt = struct('targets', (1:numel(W))', ...
             'verbose', mrstVerbose, ...
             'linsolve', @mldivide, ...
             'linsolveTOF', @mldivide, ...
             'msbasis', [], ...
             'minimize', true);
opt = merge_options(opt, varargin{:});
minu = min(u);
maxu = max(u);
%minu = cellfun(@min, u);
%maxu = cellfun(@max, u);
if or(minu < -sqrt(eps) , maxu > 1+sqrt(eps))
    warning('Controls are expected to lie in [0 1]')
end

boxLims = scaling.boxLims;
if isfield(scaling, 'obj')
    objScaling = scaling.obj;
else
    objScaling = 1;
end

% update wells:
W = control2well(u, W, 'scaling', scaling, 'targets', opt.targets);

[state, D, grd] = solveStationaryPressure(G, state, system, W, fluid, pv, T, 'objective', obj,...
                    'linsolve', opt.linsolve, 'linsolveTOF', opt.linsolveTOF, 'msbasis', opt.msbasis);
% scaled objective
sgn = -1;
if ~opt.minimize, sgn = 1; end
val = sgn*grd.objective.val/objScaling;
der = scaleGradient(sgn*grd.well(opt.targets), boxLims, objScaling);
end

% function W = control2well(u, W, boxLims, targets)
% for k = 1:numel(targets)
%     wnr = targets(k);
%     [umin, umax] = deal(boxLims(k,1), boxLims(k,2));
%     W(wnr).val = u(k)*(umax-umin)+umin;
% end
% end

function der = scaleGradient(grd, boxLims, objScaling)
dBox = boxLims(:,2)-boxLims(:,1);
der  = (dBox/objScaling).*grd;
end

