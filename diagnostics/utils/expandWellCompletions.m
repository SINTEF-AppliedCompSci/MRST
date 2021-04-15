function [xc,Wc]=expandWellCompletions(state, W, expansion, split)
%Pseudo-wells for computation of flow diagnostics for completions
%
% SYNOPSIS:
%   [xn, Wn] = expandWellCompletions(x, W, c)
%
% DESCRIPTION:
%   The routines in the flow diagnostic module compute time-of-flight,
%   tracer partition, and well-pair information associated with the whole
%   completion set of each individual well. To compute flow diagnostics for
%   subsets of the well completion, the corresponding well must be replaced
%   by a set of pseudo wells, one for each desired completion interval.
%
% REQUIRED PARAMETERS:
%   x - Reservoir and well solution structure either properly
%       initialized from functions 'initResSol' and 'initWellSol',
%       respectively, or the results from a call to function
%       'solveIncompFlow'.  Must contain valid cell interface fluxes,
%       'state.flux'.
%
%   W - Well structure as defined by function 'addWell'.
%
%   expansion - Either of two alternatives:
%
%     a)
%       nx2 vector that specifies how to expand individual wells into a
%       set of pseudo wells. That is, the completions of well number
%       c(i,1) are assigned to c(i,2) bins using a load-balanced
%       linear distribution and then the well is replaced with c(i,2)
%       pseudo wells, one for each bin of completions.
%
%     b)
%       Cell array of length n. Entry number i in this cell array should be
%       the same length as W(i).cells and contain the bins that will be
%       used to expand the well.
%
% RETURNS:
%   Wn  - Well structure containing the pseudo wells
%
%   xn  - Reservoir and solution structure corresponding
%
% SEE ALSO:
%   `computeTOFandTracer`, `computeWellPairs`

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

nw = numel(W);
if iscell(expansion)
    % Bins were provided to us, we simply count and use them.
    % split = cellfun(@(x) numel(unique(x)), expansion);
    bins = expansion;
else
    % split the completions into B bins
    splitNum = ones(1,nw);
    splitNum(expansion(:,1)) = expansion(:,2);
    
    bins = cell(nw, 1);
    for i = 1:nw
        M = numel(W(i).cells);
        b = 0:M-1;
        B = splitNum(i);
        L = floor(M ./ B);
        R = mod(M, B);
        bins{i} = max(floor(b ./(L+1)), floor((b - R)./L))+1;
    end
end
if nargin < 4
    split = cellfun(@max, bins);
end
state = validateStateForDiagnostics(state);

ind = rldecode(1:nw, split, 2);
xc = state;
xc.wellSol = state.wellSol(ind);
Wc = W(ind);
n = 0;

wfields = {'cells', 'dir', 'WI', 'dZ'};
for i=1:nw
    n = n+1;
    if split(i)==1, continue, end
    n = n-1;
    b = bins{i};
    ws = state.wellSol(i);
    for j = 1:split(i)
        n = n+1;
        ws_new = xc.wellSol(n);
        loc = b == j;
        % Setup new well
        for k = 1:numel(wfields);
            % handle well fields
            wf = wfields{k};
            Wc(n).(wf) = W(i).(wf)(loc, :);
        end
        % Make new well sol
        ws_new.flux = ws.flux(loc);
        if isfield(ws, 'cqs')
            % AD-like well state, we should sum up rates as well
            ws_new.cqs = ws.cqs(loc, :);
            ws_new.cdp = ws.cdp(loc, :);
            ws_new.cstatus = ws.cstatus(loc, :);
            % Assign new rates based on cqs field
            ws_new = setRates(ws_new);
        end
        % Well name
        name = [W(i).name ':' num2str(j)];
        Wc(n).name  = name;
        if ~isfield(ws_new, 'name'),
           ws_new.name = name;
        end
        % Write new, split well sol
        if ~isfield(xc.wellSol(n),'name'), 
           xc.wellSol(n).name = [];
        end
        xc.wellSol(n) = ws_new;
    end
end
end

function isActive = getActivePhases(ws)
    isActive = cellfun(@(x) ~isempty(x), {ws.qWs, ws.qOs, ws.qGs});
end

function ws = setRates(ws)
    isactive = getActivePhases(ws);
    names = {'qWs', 'qOs', 'qGs'};
    active = find(isactive);
    for i = 1:numel(active)
        ws.(names{active(i)}) = sum(ws.cqs(:, i));
    end
end
