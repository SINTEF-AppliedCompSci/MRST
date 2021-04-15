function [G, data, Gs, valid_ix] = readAndPrepareForPostProcessor(casenm, steps, info, precomp)
% Utility function to create input data for flow diagnostics postprocessor 
% GUI from ECLIPSE simulation output.

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

%if nargin < 4
%    precomp = [];
%end
init = readEclipseOutputFileUnFmt([casenm, '.INIT']);
grid = readEclipseOutputFileUnFmt([casenm, '.EGRID']); 


% more fields from init might be interesting
data = setStatic([], init, {'PORO', 'PERMX', 'PERMY', 'PERMZ', 'NTG', 'DEPTH'});

%[data.fluid, data.pvtDeck] = initFluidFromOutput(init);
% assume 1-1 between eclipse and mrst grid for now
G = initGridFromEclipseOutput(init, grid, 'outputSimGrid', true);
[G, Gs] = deal(G{1}, G{2});

% add porv to static props
pv = G.cells.PORV;
data.static(end+1) = struct('name', 'PORV', 'values', pv, 'limits', [min(pv), max(pv)]);
Gs.cells.PORV = G.cells.PORV;
% time in days (start and end of restart step)
startday = datenum(info.date(1, [3 2 1]));
data.time.prev = startday + info.time( max(steps-1,1) ) - info.time(1);
data.time.cur  = startday + info.time( steps ) - info.time(1);

if isempty(precomp)
    [data.states, rstrt] = convertRestartToStates(casenm, Gs, 'steps', steps, 'restartInfo', info, ...
                                                  'splitWellsOnSignChange', true);
%    if ~isfield(data.states{1}, 'flux') || isempty(data.states{1}.flux)
%        [~,~,~, info] = eclOutToGrids(init, grid);
%        data.states = addPhaseFluxes(data.states, Gs, info.T, data.fluid);
%    end
    %data.states = addFormationVolumeFactors(data.states, data.fluid, data.pvtDeck.RUNSPEC);
    %data.states = addConnectionPhaseFluxes(data.states, data.fluid, data.pvtDeck.RUNSPEC);
else
    rstrt = readEclipseRestartUnFmt(casenm, info, steps);
    data.states = cellfun(@(x)x.states{1}, precomp, 'UniformOutput', false);
end

for i = 1:numel(data.states)
    data.wells{i} = data.states{i}.wellSol;
end

valid_ix = isValidState(data.states);
if ~all(valid_ix)
    ns = nnz(~valid_ix);
    warning('Current version requires at least one open well.\n Skipping %2.0d of the %2.0d selected restart steps.', ns, numel(valid_ix));
    data.time.prev = data.time.prev(valid_ix);
    data.time.cur  = data.time.cur(valid_ix);
    data.states    = data.states(valid_ix);
    data.wells     = data.wells(valid_ix);
end

% include some more fields from restart later on
data = setDynamic(G, data, rstrt, valid_ix);

% set empty computed-field
data.computed = repmat(struct('name', '', 'values', [], 'limits', []), [0 1]);

end 

function data = setStatic(~, init, propnames)
for k = 1:numel(propnames)
    if isfield(init, propnames{k})
        v = init.(propnames{k}).values;
        data.static(k) = struct('name',    propnames{k}, ...
                                'values' , v, ...
                                'limits' , [min(v), max(v)]);
    end
end
end

function data = setDynamic(G, data, rstrt, valid_ix)
if ~any(valid_ix)
    data.dynamic = [];
else
    flds = {'PRESSURE', 'SWAT', 'SOIL', 'SGAS'};
    p = cellfun(@(x)x.pressure, data.states, 'UniformOutput', false);
    p = horzcat(p{:})/barsa;
    data.dynamic(1) = struct('name', 'PRESSURE', ...
        'values' , p, ...
        'limits' , [min(min(p)), max(max(p))]);
    
    for k = 1: size(data.states{1}.s, 2)
        nm = flds{k+1};
        vals = cellfun(@(x)x.s(:,k), data.states, 'UniformOutput', false);
        vals = horzcat(vals{:});
        % take min/max over all steps
        [minv, maxv] = deal(min(min(vals)), max(max(vals)));
        data.dynamic(k+1) = struct('name', nm, ...
            'values' , vals, ...
            'limits' , [minv, maxv]);
    end
    % additional fields of size corressponding to G.cells.num
    ix = structfun(@(fld)all(cellfun(@numel, fld)==G.cells.num), rstrt);
    fn = fieldnames(rstrt);
    % dont include sats, pressure, flux
    for pat = [flds, {'FLR', 'FLO'}]
        ix = and(ix, ~strncmp(fn, pat, numel(pat{1})));
    end
    ix = find(ix);
    k0 = numel(data.dynamic);
    for k = 1:numel(ix)
        nm = fn{ix(k)};
        vals = horzcat(rstrt.(nm){:});
        vals = vals(:, valid_ix);
        [minv, maxv] = deal(min(min(vals)), max(max(vals)));
        data.dynamic(k0+k) = struct('name', nm, ...
            'values' , vals, ...
            'limits' , [minv, maxv]);
    end
end
end

function flag = isValidState(states)
flag = true(numel(states), 1);
for k = 1:numel(states)
    ws = states{k}.wellSol;
    stat = and([ws.status], abs([ws.resv])>0);
    openPrd = any(and(stat, [ws.sign]<0));
    openInj = any(and(stat, [ws.sign]>0));
    flag(k) = openPrd || openInj;
end
end

function states = addConnectionPhaseFluxes(states, fluid, runspec)
% from output we typically only have total volume connection fluxes. 
% Approx water flux by using bW at bhp 

assert(runspec.OIL && runspec.WATER, 'Current code assumes both oil and water present')
for sk = 1:numel(states)
    ws = states{sk}.wellSol;
    for wk = 1:numel(ws)
        if size(ws(wk).flux, 2) == 1
            c = ws(wk).cells;
            resflux = ws(wk).flux;
            qwr = ws(wk).cqs(:,1)./states{sk}.b(c,1);
            qor = ws(wk).cqs(:,2)./states{sk}.b(c,2);
            if runspec.VAPOIL
                qor = qor - ws(wk).cqs(:,3).*states{sk}.rv(c)./states{sk}.b(c,2);
            end
            if ~runspec.GAS
                ws(wk).flux = [qwr, qor];
            else
                qgr = ws(wk).cqs(:,3)./states{sk}.b(c,3);
                if runspec.DISGAS
                    qgr = qgr - ws(wk).cqs(:,2).*states{sk}.rs(c)./states{sk}.b(c,3);
                    if runspec.DISGAS && runspec.VAPOIL
                        r = 1-states{sk}.rs(c).*states{sk}.rv(c);
                        [qor, qgr] = deal(qor./r, qgr./r);
                    end
                end
                ws(wk).flux = [qwr, qor, qgr];
            end
            % if remove crossflow ...
            %ws(wk).flux(resflux==0,:) = 0;
        end
    end
    states{sk}.wellSol = ws;
end
end

% function m = getModelDetails(init)
% [ih, lg] = deal(init.INTEHEAD.values, init.LOGIHEAD.values);
% unms = {'metric', 'field', 'lab'};
% m.units = unms{ih(3)};
% phnms = {'oil', 'water', 'gas'};
% phM   = [1 0 1 0 0 1 1
%          0 1 1 0 1 0 1
%          0 0 0 1 1 1 1];
% actPh = phM(:, ih(15));
% for k =1:3
%     m.(phnms{k}) = logical(actPh(k));
% end
% m.disgas = 
% m = 1;
% end
        
    

