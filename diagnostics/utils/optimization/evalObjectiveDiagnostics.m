function [val, der, W, state, D] = evalObjectiveDiagnostics(u, obj, state, system, G, fluid, pv, T, W, scaling, varargin)
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

