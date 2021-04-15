function [eqs, cq_s, sol, Rw] = getWellContributions(W, sol, pBH, q_s, p, rho_s, b, r, rMax, m, varargin)
% [eqs, cq_s, sol, closedConns] = getWellContributions(W, sol, pBH, q_s, p, rho_s, b, r, rMax, m, varargin)
% INPUT:
% W     : well-structure
% sol   : well solution structure
% pBH   : well BHPs
% qt_s  : well total volumerates at standard conds
% p     : grid-cell pressures at well connections
% rho_s : vector of phase densities at std conds
% b     : cell array of phase b-factors at connections
% r     : cell array of phase solution-ratios at connections (rs,rv,...)
% rMax  : cell array of maximum phase solution-ratios at connections
% m     : cell array of phase mobilities at connections
% OPTIONAL INPUT
% model      : model beeing run (if empty, model is chosen based on
%              input), currently treats
%                OW: oil-water no mixing (r = {})
%                3P: oil-water-gas no mixing (r = {})
%                BO: black-oil, dry gas, live oil (r = {rs})
%                VO: black-oil, wet gas, live oil (r = {rs,rv})
% iteration  : current Newton iteration. Well pressure-drop will be updated
%              only if iteration == 1 (default -1)
% allowControlSwitching: if set to true, well controls switch whenever
%              supplied limits are violated (default true)
% allowWellSignChange: if set to true producers are allowed to become
%              injectors and vice versa (default false)
% allowCrossFlow: if set to false, connections with cross-flow will be shut
%              (default false)
%
% OUTPUT:
% eqs   : numPh+1 well equations
% cq_s  : cell array of connection phase volumerates at standard conds
% sol   : updated well solution structure (pressure-drop will be updated if
%         opt.iteration == 1, type will be updated if well limits are
%         violated.
% closedConns: logical vector with true for every closed connection (empty
%         if allowCrossFlow is set to true)

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

opt = struct('model',                   []      ,...
             'iteration',               -1      ,...
             'allowControlSwitching',   true    ,...
             'allowWellSignChange',     false   ,...
             'allowCrossFlow',          true    ,...
             'Verbose',                 mrstVerbose);
opt = merge_options(opt, varargin{:});

% If model is not given, decide based on input (this is a temporary
% solution, model should always be given in some way)
model = opt.model;
if isempty(opt.model)
    model = getModel(numel(b), numel(r));
end
nPh = numel(b); % # phases

%--------------------------------------------------------------------------
% Wellbore pressure-drop should be updated only first iteration in each
% time-step
if opt.iteration == 1
    sol = updateConnDP(W, sol, b, rMax, rho_s, model);
end
% For now, give warning if iteration number is not supplied (i.e.,
% iteration = -1)
if opt.iteration < 0
    warning(['Iteration number is not passed on to getWellContributions,', ...
             'this may indicate welbore pressure-drop will never be updated']);
end
% check for non-initialized bhp due to non-default reference depth:
if any(~isfinite(value(pBH)))
    [sol, pBH] = initializeBHP(sol, pBH, p);
end

%--------------------------------------------------------------------------
% Check and switch control limits if allowed
if opt.allowControlSwitching 
    [sol, withinLims] = updateControls(W, sol, pBH, q_s, model);
    [pBH, q_s] = updateVars(pBH, q_s, sol, W, withinLims);
end

%--------------------------------------------------------------------------
% if solveWellEqs:
%    bla bla
% end

% Compute well-flow and get first nPh well equations
[eqs, cq_s, mix_s, status, cstatus, Rw] = computeWellContributions(...
                    W, sol, pBH, q_s, p, b, r, m, model, ...
                    opt.allowWellSignChange, opt.allowCrossFlow);

% Finally impose supplied controls
eqs{nPh+1} = getControlEquations(sol, pBH, q_s, status, mix_s, model);
%--------------------------------------------------------------------------
% Update well properties which are not primary variables
nConn       = cellfun(@numel, {W.cells})'; % # connections of each well
perf2well   = rldecode((1:numel(W))', nConn);
toDb = @(x)cellfun(@value, x, 'UniformOutput', false);
cq_sDb = cell2mat(toDb(cq_s));
for wnr = 1:numel(sol)
    ix = perf2well==wnr;
    sol(wnr).cqs     = cq_sDb(ix,:);
    sol(wnr).cstatus = cstatus(ix);
    sol(wnr).status = status(wnr);
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function model = getModel(np, nr)
models = {''  , ''  , ''
          'OW', ''  , ''
          '3P', 'BO', 'VO'};
model = models{np, nr+1};
end

function  [pBH, q_s] = updateVars(pBH, q_s, sol, W, withinLims)
if all(withinLims)
    return
end
for k = 1:numel(sol)
    if ~withinLims(k)
        tp = sol(k).type;
        v  = sol(k).val;
        switch tp
            case 'bhp'
                pBH = assignVal(pBH,v,k);
            case 'rate'
                q_s{1} = assignVal(q_s{1}, v*W(k).compi(1), k);
                q_s{2} = assignVal(q_s{2}, v*W(k).compi(1), k);
                if numel(q_s)>2
                    q_s{3} = assignVal(q_s{3}, v*W(k).compi(1), k);
                end
            case 'orat'
                q_s{2} = assignVal(q_s{2}, v , k);
            case 'wrat'
                q_s{1} = assignVal(q_s{1}, v , k);
            case 'grat'
                q_s{3} = assignVal(q_s{3}, v , k);
        end % No good guess for qOs, etc...
    end
end
end

function [sol, pBH] = initializeBHP(sol, pBH, p)
% We should have sol(k).bhp = pBH(k) at this point
pd   = value(p);
for k = 1:numel(sol)
    if ~isfinite(sol(k).bhp)
        v = pd(k) + 5*sol(k).sign*barsa - sol(k).cdp(1);
        pBH = assignVal(pBH, v, k);
        sol(k).bhp = v;
    end
end
end

function x = assignVal(x, v, inx)
if isa(x, 'ADI')
    x.val(inx) = v;
else
    x(inx)=v;
end

end
