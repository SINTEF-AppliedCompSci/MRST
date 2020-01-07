function deckmrst = model2Deck(model, schedule, varargin)
% Create deck-structure from MRST model and schedule. 
%
% SYNOPSIS:
%   deckmrst = model2Deck(model, schedule)
%   deckmrst = model2Deck(model, schedule, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function creates a deck-structure amenable for writing input files
%   to other simulators by e.g., writeDeck.
%
% REQUIRED PARAMETERS:
%   model     -  MRST model  
%
%   schedule  - MRST schedule
%
% OPTIONAL PARAMETERS:
%   'unit'   - Output deck-structure choice of unit ('metric', 'field' or 
%              'si'). Default is 'si'.
%   'state0' - Initial state. If not empty, written to SOLUTION section. 
%   'deck'   - Deck structure. All fields of deck that are not explicitly 
%              created in this function will be copied to deckmrst.
%
%   

opt = struct('state0',      [], ...
             'deck',        [], ...
             'unit',        'si');
opt  = merge_options(opt, varargin{:});
deck = opt.deck;
if ~isempty(deck)
    deck = convertDeckUnits(deck, 'outputUnit', 'si');
end

%% RUNSPEC
RUNSPEC.DIMENS      = [model.G.cells.num,1,1];
RUNSPEC.OIL         = model.oil;
RUNSPEC.GAS         = model.gas;
RUNSPEC.WATER       = model.water;
RUNSPEC.DISGAS      = model.disgas;
RUNSPEC.VAPOIL      = model.vapoil;
RUNSPEC.UNIFOUT     = 1;
RUNSPEC.SI          = 1;
RUNSPEC.WELLDIMS    = getWellDims(schedule);
%RUNSPEC.START       = 734813;
RUNSPEC.TABDIMS     = getTabDims(deck); 
%SATOPTS=[]
%SATOPTS.DIRECT=0;
%SATOPTS.HYSTER=0;
%SATOPTS.IRREVERS=0;
%RUNSPEC.SATOPTS=SATOPTS;
%RUNSPEC.GRIDOPTS={'YES'  [0]  [0]};
%RUNSPEC.REGDIMS= [1 1 1 0 0 1 0 0 0];
%RUNSPEC.EQLDIMS= [1 100 50 1 50];
if isfield(deck, 'RUNSPEC')
    RUNSPEC = addAdditional(RUNSPEC, deck.RUNSPEC);
end

%% GRID
nc = model.G.cells.num;
GRID.INIT = 1;
% make 1D grid where all connections are NNC
GRID.DX = model.G.cells.volumes;
GRID.DY = ones(model.G.cells.num, 1);
% aprox TOPS
depth = model.G.cells.centroids(:,3);
if isfield(model.G.faces, 'nodepos')
    [~,~,dz] = cellDims(model.G, (1:nc)');
else
    dz = ones(nc,1);
end
GRID.DZ = dz;
GRID.TOPS = depth - dz/2;

% permeabilities (not used)
GRID.PERMX = ones(nc,1);
GRID.PERMY = ones(nc,1);
GRID.PERMZ = ones(nc,1);
GRID.PORO  = ones(nc,1);
GRID.ACTNUM   =int32(ones(nc,1));
GRID.cartDims = RUNSPEC.DIMENS;
% make NNC - list: 
T = model.operators.T;
nf = numel(T);
GRID.NNC = [model.operators.N(:,1), ones(nf,2), model.operators.N(:,2), ones(nf,2), T];

%% EDIT
% set connections in 1D-grid to zero trans
EDIT.TRANX = zeros(nc,1);
EDIT.PORV  = model.operators.pv;
EDIT.DEPTH = depth;

%% PROPS (from input deck)
PROPS = [];
if isfield(deck, 'PROPS')
    PROPS = deck.PROPS;
end

%% REGIONS
REGIONS=struct();
if isfield(deck, 'REGIONS')
    REGIONS = deck.REGIONS;
end

%% SUMMARY
SUMMARY=struct();
if isfield(deck, 'SUMMARY')
    SUMMARY = addAdditional(SUMMARY, deck.SUMMARY);
end

%% SOLUTION
SOLUTION=struct();
if ~isempty(opt.state0)
    state0 = opt.state0;
    SOLUTION.PRESSURE = state0.pressure;
    wix = model.getPhaseIndex('W');
    if ~isempty(wix)
        SOLUTION.SWAT = state0.s(:,wix);
    end
    oix = model.getPhaseIndex('O');
    if ~isempty(oix)
        SOLUTION.SOIL = state0.s(:,oix);
    end
    gix = model.getPhaseIndex('G');
    if ~isempty(gix)
        SOLUTION.SGAS = state0.s(:,gix);
    end
    if isfield(state0, 'rs')
        SOLUTION.RS = state0.rs;
    end
    if isfield(state0, 'rv')
        SOLUTION.RV = state0.rv;
    end
elseif isfield(deck, 'SOLUTION')
    SOLUTION = addAdditional(SOLUTION, deck.SOLUTION);
end

%% SCHEDULE
SCHEDULE = convertScheduleToDeck(model, schedule, 'linearIndex', true, 'reduceOutput', true);


%% unhandeled ??
UnhandledKeywords=[];
UnhandledKeywords.SUMMARY = {'ALL'};
UnhandledKeywords.SCHEDULE= {'OPTIONS'  'RPTRST'};

%% assemble
deckmrst.GRID     = GRID;
deckmrst.RUNSPEC  = RUNSPEC;
deckmrst.EDIT     = EDIT;
deckmrst.PROPS    = PROPS;
deckmrst.REGIONS  = REGIONS;
deckmrst.SUMMARY  = SUMMARY;
deckmrst.SOLUTION = SOLUTION;
deckmrst.SCHEDULE = SCHEDULE;
deckmrst.UnhandledKeywords = UnhandledKeywords;

%% finally convert to requested units 
if ~isempty(opt.unit)
    deckmrst = convertDeckUnits(deckmrst, 'outputUnit', opt.unit);
end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function v = getWellDims(schedule)
%    1        2       3       4       5        6      7      8     9   10  11  12  13 14   
% nwells    nconn   ngrp  nwellgrp   ntfip    nrpvt  nrvpvt ntendp                    
v =  [10      20     10      20       5       10      5      4     3   0   1    1  10 201];
if isfield(schedule.control(end), 'W')
    W = schedule.control(end).W;
    v(1) = numel(W);
    v(2) = max(arrayfun(@(w)numel(w.cells), W));
    v(4) = v(1);
end
end

% -------------------------------------------------------------------------
function v = getTabDims(deck)
%    1        2       3       4       5        6      7      8     9   10  11  12  13 14 15  
% ntsfun    ntpvt   nssfun  nppvt   ntfip    nrpvt  nrvpvt ntendp                 ntrocc     
v = [1       1       50      50      16       30     30      1     1   1   10   1  -1  0 1];
if isempty(deck) || ~isfield(deck, 'PROPS')
    dispif(mrstVerbose, 'PROPS not provided, defaulting TABDIMS\n');
else
    nms = fieldnames(deck.PROPS);
    % get some saturation function
    ix = strncmp({'SWOF'}, nms, 4) | strncmp({'SWFN'}, nms, 4);
    if any(ix)
        satfld = nms{find(ix, 1, 'first')}; 
        v(1) = numel(deck.PROPS.(satfld));
    end
    % get some pvt
    ix = strncmp({'PVTO'}, nms, 4) | strncmp({'PVDO'}, nms, 4);
    if any(ix)
        pvtfld = nms{find(ix, 1, 'first')};
        v(2) = numel(deck.PROPS.(pvtfld));
    end
end
end

% -------------------------------------------------------------------------
function s = addAdditional(s, sa)
flds = setdiff(fieldnames(sa), fieldnames(s));
for k = 1:numel(flds)
    s.(flds{k}) = sa.(flds{k});
end
end

        
