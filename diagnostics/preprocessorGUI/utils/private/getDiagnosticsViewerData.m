function Data = getDiagnosticsViewerData(models,wells,initStates,varargin)
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

opt = struct('wellSolFields', [],...
    'state0', [],...
     'D',     [], ...
     'state', [], ...
     'WP',    [], ...  
     'wellCommunication',     [], ...
     'G',     []) ;
opt = merge_options(opt, varargin{:});

[state, D, WP, wellCommunication] = deal([]);
Data     = cell(1,numel(models));
fldnames = {'PORO', 'PERMX', 'PERMY', 'PERMZ', 'DEPTH', 'PORV'};
for k = 1:numel(models)
    model = models{k};
    state0 = initStates{k};
    
    if ~isempty(opt.state0),  state0 = opt.state0{k}; end
    if ~isempty(opt.D),       D      = opt.D{k};      end
    if ~isempty(opt.state),   state  = opt.state{k};  end  
    if ~isempty(opt.WP),      WP     = opt.WP{k};     end 
    if ~isempty(opt.wellCommunication)
        wellCommunication = opt.wellCommunication{k};
    end       

    tic
    [states, diagnostics] = ...
        computePressureAndDiagnostics(model, 'wells', wells{k}, ...
            'state0', state0, 'D', D, 'state', state, ...
            'WP', WP, 'wellCommunication', wellCommunication);
    toc
    init.PORO.values  = model.rock.poro;
    init.PERMX.values = model.rock.perm(:,1);
    init.PERMY.values = model.rock.perm(:,2);
    init.PERMZ.values = model.rock.perm(:,3);
    init.DEPTH.values = model.G.cells.centroids(:,3);
    init.PORV.values  = poreVolume(model.G,model.rock);
    
    [injColors,prodColors] = getWellPlotColours(diagnostics.wellCommunication);
    
    Data{k}.diagnostics = diagnostics;
    Data{k}.states     = states;
    Data{k}.static     = setStatic(init, fldnames);
    Data{k}.dynamic    = setDynamic(states);
    Data{k}.injColors  = injColors;
    Data{k}.prodColors = prodColors;
    Data{k}.wells      = wells{k};
    Data{k}.wellSol    = states.wellSol;
    Data{k}.wsdata     = getDataFromIncompTPFAWellSol(wells{k},states.wellSol);
    Data{k}.Gs         = model.G;
    if isempty(opt.G)
        Data{k}.G = model.G;
    else
        if iscell(opt.G)
            Data{k}.G = opt.G{k};
        else
            Data{k}.G = opt.G;
        end
    end
end

% update limits
for k=1:numel(Data{1}.dynamic)
    minv = inf;
    maxv = -inf;
    for j=1:numel(models)
        minv = min(minv, min(Data{j}.dynamic(k).values(:)));
        maxv = max(maxv, max(Data{j}.dynamic(k).values(:)));
    end
    for j = 1:numel(models)       
         Data{j}.dynamic(k).limits = [minv, maxv];
    end
end


function wsdata = getDataFromIncompTPFAWellSol(W,ws,varargin)
opt = struct('wellSolFields', []);
opt = merge_options(opt, varargin{:});

namelist = arrayfun(@(x) x.name, W, 'UniformOutput', false);
if isempty(opt.wellSolFields)
    keywordlist = fieldnames(ws(1));
else
    keywordlist = opt.wellSolFields;
end

% Generate data matrix
propIdx = 1;
nm      = cell(1,numel(keywordlist));
props   = cell(1,numel(keywordlist));

for i = 1:numel(keywordlist)
    kw = keywordlist{i};
    [nm{propIdx},unit] = getUnit(kw);
    if ~strcmp(unit, '-')
        % Sum all fluxes and cdp for each well.
        data = arrayfun(@(x) sum(sum(x.(kw))),ws);
        data = convertTo(data,unit);
        if strcmp(nm{propIdx},'m^3/day')
            data = abs(data);
        end
        wsdata.(kw)    = data;
        props{propIdx} = kw;
        propIdx        = propIdx + 1;
    end
end

wsdata.props     = props;
wsdata.wellNames = namelist;
wsdata.units     = nm;
wsdata.times     = 1;
wsdata.timesteps = 1;


function [nm,u] = getUnit(kw)
switch kw 
    case {'pressure', 'bhp'}
        nm  = 'barsa';  u = barsa();
    case 'cdp'
        nm  = 'barsa';  u = barsa();
    case 'flux'
        nm = 'm^3/day'; u = 1./day();
    otherwise 
        nm = '';        u = '-';
end


function static = setStatic(init, propnames)
%for k = 1:numel(propnames)
%    if isfield(init, propnames{k})
%        v = init.(propnames{k}).values;
%        static(k) = struct('name',    propnames{k}, ...
%            'values' , v, ...
%            'limits' , [min(v), max(v)]);
%    end
%end
names  = propnames(isfield(init, propnames));
vals   = cellfun(@(nm) init.(nm).values, names, 'UniformOutput', false);
lims   = cellfun(@(x) [min(x), max(x)],  vals,  'UniformOutput', false);
static = struct('name', names, 'values', vals, 'limits', lims);


function dynamic = setDynamic(state)

flds = {'PRESSURE', 'SWAT', 'SOIL', 'SGAS'};
p = state.pressure/barsa;
dynamic(1) = struct('name', 'PRESSURE', ...
    'values' , p, ...
    'limits' , [min(min(p)), max(max(p))]);

for k = size(state.s, 2):-1:1
    nm = flds{k+1};
    vals = state.s(:,k);
    % take min/max over all steps
    minv = min(vals(:));
    maxv = max(vals(:));
    dynamic(k+1) = struct('name', nm, ...
        'values' , vals, ...
        'limits' , [minv, maxv]);
end


% ------ Set unique colors for wells --------------------------
function [injColors,prodColors] = getWellPlotColours(wellCommunication)   
% Avoid first color (black) and add an extra color representing
% potential contributions from the reservoir
[nInj,nProd] = size(wellCommunication);
cmap = tatarizeMap(nInj+nProd+2);
injColors  = cmap([2:nInj+1 end],:);
prodColors = cmap(nInj+2:end,:);
    

