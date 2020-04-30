function Data = getDiagnosticsViewerData(models,wells,varargin)

opt = struct('wellSolFields', [],...
    'state0', [],...
    'fluid',  [], ...
     'D',     [], ...
     'state', [], ...
     'WP',     [], ...  
     'wellCommunication',        []);
opt = merge_options(opt, varargin{:});


for k = 1:numel(models)
    
    model = models{k};
    if ~isempty(opt.state0)
        state0 = opt.state0{k};
    else
        state0 = [];
    end
    if ~isempty(opt.fluid)
        fluid = opt.fluid{k};
    else
        fluid  = [];
    end    
    if ~isempty(opt.D)
        D = opt.D{k};
    else
        D  = [];
    end   
    if ~isempty(opt.state)
        state = opt.state{k};
    else
        state  = [];
    end   
    if ~isempty(opt.WP)
        WP = opt.WP{k};
    else
        WP  = [];
    end       
    if ~isempty(opt.wellCommunication)
        wellCommunication = opt.wellCommunication{k};
    else
        wellCommunication  = [];
    end       

    tic
    [states, diagnostics] = computePressureAndDiagnostics(model, 'wells', wells{k}, ...
        'state0', state0, 'fluid', fluid, 'D', D, 'state', state, ...
        'WP', WP, 'wellCommunication', wellCommunication);
    toc
    init.PORO.values = model.rock.poro;
    init.PERMX.values = model.rock.perm(:,1);
    init.PERMY.values = model.rock.perm(:,2);
    init.PERMZ.values = model.rock.perm(:,3);
    init.DEPTH.values = model.G.cells.centroids(:,3);
    init.PORV.values = poreVolume(model.G,model.rock);
    
    static = setStatic(init,  {'PORO', 'PERMX', 'PERMY', 'PERMZ', 'DEPTH', 'PORV'});
    dynamic = setDynamic(states);
    %                 wellSols{k} = setWellSols(states{k},wells{k},opt.wellSolFields,1);
    
    [injColors,prodColors] = getWellPlotColours(diagnostics.wellCommunication);
    
    Data{k}.diagnostics = diagnostics;
    Data{k}.states = states;
    Data{k}.static = static;
    Data{k}.dynamic = dynamic;
    Data{k}.injColors = injColors;
    Data{k}.prodColors = prodColors;
    Data{k}.wells = wells{k};
    Data{k}.wellSol = states.wellSol;
    
    Data{k}.wsdata = getDataFromIncompTPFAWellSol(wells{k},states.wellSol);
    
    Data{k}.G = model.G;
    
    
    
end


% update limits


for k = 1:numel(Data{1}.dynamic)
    minv = zeros(numel(models),1);
    maxv = zeros(numel(models),1);    
    for j = 1:numel(models)          
        minv(j) = min(min(Data{j}.dynamic(k).values));
        maxv(j) = max(max(Data{j}.dynamic(k).values));
    end
    minmax = [minv maxv];
    minall = min(min(minmax));
    maxall = max(max(minmax));    
    
    for j = 1:numel(models)       
         Data{j}.dynamic(k).limits = [minall, maxall];
    end
    
end


%             Data.wellInfo = wells;
%             Data.wellSols = wellSols;
%
%             stop = 'here'

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
for i = 1:numel(keywordlist)
    
    kw = keywordlist{i};
    
    % Sum all fluxes and cdp for each well.
    data = arrayfun(@(x) sum(x.(kw)),ws);

    [nm{propIdx},unit] = getUnit(kw);
    

     data =  convertTo(data,unit);
     if strcmp(nm{propIdx},'m^3/day')
         data = abs(data);
     end
     wsdata.(kw) = data;

    
    props{propIdx} = kw;
    propIdx = propIdx + 1;
    
end

wsdata.props = props;
wsdata.wellNames = namelist;
wsdata.units = nm;
wsdata.times = 1;
wsdata.timesteps = 1;



function [nm,u] = getUnit(kw)
switch kw 
    case 'pressure'
        nm  = 'barsa';
        u = barsa();
    case 'cdp'
        nm  = 'barsa';
        u = barsa();
    case 'flux'
        nm = 'm^3/day';
        u = 1./day();
    otherwise 
        nm = ''
        u = '-';
end



function static = setStatic(init, propnames)
for k = 1:numel(propnames)
    if isfield(init, propnames{k})
        v = init.(propnames{k}).values;
        static(k) = struct('name',    propnames{k}, ...
            'values' , v, ...
            'limits' , [min(v), max(v)]);
    end
end



function dynamic = setDynamic(state)

flds = {'PRESSURE', 'SWAT', 'SOIL', 'SGAS'};
p = state.pressure/barsa;
dynamic(1) = struct('name', 'PRESSURE', ...
    'values' , p, ...
    'limits' , [min(min(p)), max(max(p))]);

for k = 1: size(state.s, 2)
    nm = flds{k+1};
    vals = state.s(:,k);
    % take min/max over all steps
    [minv, maxv] = deal(min(min(vals)), max(max(vals)));
    dynamic(k+1) = struct('name', nm, ...
        'values' , vals, ...
        'limits' , [minv, maxv]);
end




function [injColors,prodColors] = getWellPlotColours(wellCommunication)
    
    % ------ Set unique colors for wells --------------------------
    % Avoid first color (black) and add an extra color representing
    % potential contributions from the reservoir
    [nInj,nProd] = size(wellCommunication);
    cmap = tatarizeMap(nInj+nProd+2);
    injColors  = cmap([2:nInj+1 end],:);
    prodColors = cmap(nInj+2:end,:);
    

function wellSols = setWellSols(state, W, wellSolFields, time)

ws = state.wellSol;
% Read names
% Assumes that all wellSols contain the same fields at the wellSol at
% timestep 1.

namelist = arrayfun(@(x) x.name, W, 'UniformOutput', false);


% Read keywords
if isempty(wellSolFields)
    keywordlist = fieldnames(ws(1));
else
    keywordlist = wellSolFields;
end

% Generate data matrix

propIdx = 1;
for i = 1:numel(keywordlist)
    
    kw = keywordlist{i};
    
    try
        dat = [ws(1).(kw)];
    catch
        continue
    end
    
    if ~isa([ws(1).(kw)],'double') || isempty([ws(1).(kw)]) || numel([ws(1).(kw)])~= numel(namelist)
        continue
    end
    
    wsdata.(kw) = zeros(numel(ws),numel(ws{1}));
    [nm{propIdx},unit] = getUnit(kw);
    
    for j = 1:numel(ws)
        data =  convertTo([ws(j).(kw)],unit);
        if strcmp(nm{propIdx},'m^3/day')
            data = abs(data);
        end
        wsdata.(kw)(j,:) = data;
    end
    
    props{propIdx} = kw;
    propIdx = propIdx + 1;
    
end

timesteps = 1:1:numel(time);


wsdata.props = props;
wsdata.wellNames = namelist;
wsdata.times = time;
wsdata.timesteps = timesteps;
wsdata.units = nm;




