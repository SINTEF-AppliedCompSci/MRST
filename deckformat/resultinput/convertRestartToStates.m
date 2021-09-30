function [states, rstrt] = convertRestartToStates(fn, G, varargin)
% states = convertRestartToStates(fn, G, varargin)
% Produce MRST-compatible states from eclipse restart (and maybe summary) data
% given by prefix fn.
% Optional input:
%         'additionalFields'    : not yet supported
%         'includeWellSols'     : (true)
%         'includeFluxes'       : (true) include wellSols if corresponding
%                                 fields are present in restart
%         'neighbors',          : ([]) use if other than G.faces.neighbors
%         'wellSolsFromRestart' : (true) read from restart rather than
%                                  summary
%         'consistentWellSols'  : (false) loop through to make number of
%                                 wells and perforations uniform
%
%         'steps',              : ([]) restart step numbers to read (1:n), if
%                                 empty, read all. should be used with
%                                 non-unified ouput. NOTE: Step number may
%                                 not correspond to file extension.

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

opt     = struct('additionalFields',    {{}}, ...
                 'includeWellSols',     true, ...
                 'includeFluxes',       true, ...
                 'includeMobilities',   true, ...
                 'includeAquifers',     false, ...
                 'neighbors',           [],   ...
                 'wellSolsFromRestart', true, ...
                 'consistentWellSols',  true,...
                 'splitWellsOnSignChange', false, ...
                 'steps',               [], ...
                 'restartInfo',         [], ...
                 'addTrajectory',       true, ...
                 'removeClosedWells',   true, ...
                 'removeCrossflow',     true, ...
                 'setToClosedTol',      0*meter^3/day);
opt     = merge_options(opt, varargin{:});

N    = opt.neighbors;
if isempty(N)
    N = G(1).faces.neighbors;
end

[dd, nm, ext] = fileparts(fn);
if strcmp(ext, '.UNRST')
    fn = fullfile(dd, nm);
end

% check if we have rsspec-file available
hasInfo = true;
if isempty(opt.restartInfo)
    if exist([fn, '.RSSPEC'], 'file')
        opt.restartInfo = processEclipseRestartSpec(fn);
    else
        hasInfo = false;
    end
end

if hasInfo
    rstrt = readEclipseRestartUnFmt(fn, opt.restartInfo, opt.steps);
else
    dispif(mrstVerbose, ...
        'Restart specification missing, this may potentially increase processing time\n');
    isUnified = exist([fn,'.UNRST'], 'file');
    if isUnified
        tmp = readEclipseOutputFileUnFmt([fn,'.UNRST']);
        ns  = numel(tmp.SEQNUM.values);
    else
       % assume the file is for a single restart step  
       tmp = readEclipseOutputFileUnFmt(fn);
       ns  = 1;
    end
    fnames = fieldnames(tmp);
    rstrt  = [];
    for i = 1 : numel(fnames)
        if ~strcmp(tmp.(fnames{i}).type, 'MESS')
            V = tmp.(fnames{i}).values;
            if  mod(numel(V), ns) ~= 0
                if mrstVerbose()
                    fprintf('Skipping entry %s\n', fnames{i});
                end
                continue
            end
            V = reshape(V,[],ns);
            if isempty(opt.steps)
                opt.steps=1:ns;
            end
            nsout = numel(opt.steps);
            V = V(:,opt.steps);
            V = mat2cell(V, size(V,1), ones([1, nsout]));
            rstrt.(fnames{i}) = reshape(V,[],1);
        end
    end
    % for old opm (???)
    if isfield(rstrt,'XCON') && size(rstrt.XCON{1}, 1) == 1
        rstrt.XCON = applyFunction(@transpose, rstrt.XCON);
    end
end

[opt, phn, unit, tr, na, isECLIPSE] = checkAndProcessInput(fn, rstrt, opt);

if opt.includeFluxes
    % check if magic index-map is present in grid, if not compute it
    if isfield(G.faces, 'ix')
        ix = G.faces.ix;
    else
        if ~opt.nncPresent
            NNC = [];
        else
            % Get NNC from grid-file (use EGRID if present, else GRID)
            grdFile = [fn,'.EGRID'];
            if isempty(dir(grdFile))
                grdFile = [fn,'.GRID'];
            end
            grd = readEclipseOutputFileUnFmt(grdFile);
            if isfield(grd, 'NNC1')
                NNC = [grd.NNC1.values, grd.NNC2.values];
            else
                % check if NNC is given in INIT-file
                initFile = [fn,'.INIT'];
                init  = readEclipseOutputFileUnFmt(initFile);
                if isfield(init, 'NNC1')
                    NNC = [init.NNC1.values, init.NNC2.values];
                else
                    error('NNC fluxes given in restart, but indices for NNC cells neither found in GRID or INIT files')
                end
            end
        end
        ix = eclFaceToFace(G, NNC, 'neighbors', opt.neighbors);
    end
 
    % Get all reservoir flux fieldnames in rstrt
    fluxNms = getResFluxNms(phn, isECLIPSE);
end

if opt.includeWellSols
    if opt.wellSolsFromRestart
        % Create ijk to active maping
        glob2act = zeros(prod(G.cartDims), 1);
        glob2act(G.cells.indexMap) = (1:G.cells.num)';
        ijk2act = @(ijk)glob2act( ...
            sub2ind(G.cartDims, ijk(:,1), ijk(:,2), ijk(:,3)) );
    else % get wellSols from summary
        [wellSols, tm] = convertSummaryToWellSols(fn, unit);
        wellSols = wellSols( ms2rs(tm, tr));
    end
end

% Set up default states
states = repmat({struct('pressure', zeros(na, 1),   ...
                        's', zeros(na, numel(phn)), ...
                        'flux', [], 'wellSol', [],  ...
                        'time', '-1', 'date', [])}, [numel(tr), 1]);

% Setup saturation names for current config
satnms = cellfun(@(x)['S',x], phn, 'UniformOutput', false);
satix  = isfield(rstrt, satnms); % occuring fields

% If not all saturation fields are present, assert only one is missing
assert(nnz(~satix)<=1, 'Saturation output found for less than nPh-1 phases');
dispif(mrstVerbose, 'Converting restart to mrst-states:     ')
for k = 1:numel(tr)
    dispif(mrstVerbose, '\b\b\b\b%3d%%', round(100*k/numel(tr)));
    % pressure
    p = rstrt.PRESSURE{k};
    if ~isempty(p), states{k}.pressure = convertFrom(p, unit.p); end
    % saturation
    for ii = find(satix)
        states{k}.s(:,ii) = rstrt.(satnms{ii}){k};
    end
    if any(~satix) % at most one
        states{k}.s(:,~satix) = 1 - sum(states{k}.s,2);
    end
    % time
    states{k}.time = convertFrom(tr(k), unit.t);

    % date
    states{k}.date = rstrt.INTEHEAD{k}(65:67)';
    if isfield(rstrt, 'SEQNUM')
        states{k}.index = rstrt.SEQNUM{k};
    end

    % rs/rv
    if isfield(rstrt, 'RS')
        states{k}.rs = convertFrom(rstrt.RS{k}, unit.qg/unit.ql);
    end

    if isfield(rstrt, 'RV')
        states{k}.rv = convertFrom(rstrt.RV{k}, unit.ql/unit.qg);
    end

    % polymer stuff
    if isfield(rstrt, 'POLYMER')
        states{k}.c = rstrt.POLYMER{k};
    end

    if isfield(rstrt, 'PADS')
        states{k}.ads = rstrt.PADS{k};
    end

    % compositional features
    states{k} = convertCompositions(states{k}, 'XMF', rstrt, k, 'x');
    states{k} = convertCompositions(states{k}, 'YMF', rstrt, k, 'y');
    states{k} = convertCompositions(states{k}, 'ZMF', rstrt, k, 'components');

    % fluxes
    if opt.includeFluxes
        f = zeros(size(N,1), numel(phn));
        if ~isempty(rstrt.(fluxNms(1).I){k})
            for phi = 1:numel(phn)
                f(ix.iI, phi) = ix.sI.*rstrt.(fluxNms(phi).I){k}(ix.iIe);
                f(ix.iJ, phi) = ix.sJ.*rstrt.(fluxNms(phi).J){k}(ix.iJe);
                if ~isempty(ix.iK)
                    f(ix.iK, phi) = ix.sK.*rstrt.(fluxNms(phi).K){k}(ix.iKe);
                end
                if isfield(ix, 'iN') && ~isempty(ix.iN)
                    f(ix.iN, phi) = ix.sN.*rstrt.(fluxNms(phi).N){k}(ix.iNe);
                end
            end
        end
        nbad = nnz(~isfinite(f));
        if nbad > 0
            dispif(mrstVerbose, ...
                'Reading restart resulted in %d non-finite flux-values.\n', nbad);
            f(~isfinite(f)) = 0;
        end
        states{k}.flux = convertFrom(f, unit.qr);
    end
    
    % wellSols
    if opt.includeWellSols
        if opt.wellSolsFromRestart
            states{k}.wellSol = createWellSol(rstrt, k, ijk2act);   
        else
            states{k}.wellSol = wellSols{k};
        end
    end
    
    % aquifers
    if opt.includeAquifers
        aq = createAquiSol(rstrt, k, ijk2act);
        if ~isempty(aq)
            states{k}.aquiSol= aq;
        end
    end
end
dispif(mrstVerbose, ',  done\n')
% check for eMap-field, and if present, reduce data

if isfield(G.cells, 'eMap')
    fns = {'pressure', 's', 'rs', 'rv'};
    for k1 = 1:numel(states)
        for k2 = 1:numel(fns)
            fn = fns{k2};
            if isfield(states{k1}, fn)
                states{k1}.(fn) = states{k1}.(fn)(G.cells.eMap, :);
            end
        end
    end
end

% make consistent wellSols
if opt.consistentWellSols && opt.includeWellSols
    states = makeWellSolsConsistent(states);
    % handle closing of wells and crossflow
    states = processWellStates(states, opt);
end

% if opt.addTrajectory
%     for k = 1:numel(states)
%         states{k}.wellSol = addTrajectories(states{k}.wellSol, G, 10);
%     end
% end

% Compute cell frac-flows as avg over cell influx + outflux and incorrectly
% "mob".
if opt.includeMobilities
    ie = ~any(N==0,2);
    n  = double(N(ie,:));
    nn = size(n,1);
    C  = sparse(n, repmat((1:nn)' , [1 2]), 1, G.cells.num, nn);
    for k = 1:numel(states)
        qa = C*abs(states{k}.flux(ie,:));
        if isfield(states{k}, 'wellSol') && ~isempty(states{k}.wellSol)
            wc = vertcat(states{k}.wellSol.cells);
            % need some tricks since only total well-flux is given
            % q_ph = q_ph^r + q_ph^w
            qwr    = qa(wc,:);
            qw_tot = abs(vertcat(states{k}.wellSol.flux));
            sqw    = max(sum(qwr, 2) + qw_tot, 100*eps*max(qw_tot));
            qa(wc,:) = bsxfun(@rdivide, qwr, 1-qw_tot./sqw);
        end
        sq = sum(qa,2);
        sq = max(sq, 100*eps*max(sq));
        states{k}.mob = bsxfun(@rdivide, qa, sq);
    end
end
end

%--------------------------------------------------------------------------

function u = getUnits(unit)
switch unit
    case 'metric'
        u.len = meter;
        u.p   = barsa;
        u.ql  = meter^3/day;
        u.qg  = meter^3/day;
        u.qr  = meter^3/day;
        u.t   = day;
    case 'field'
        u.p  = psia;
        u.ql  = stb/day;
        u.qg  = 1000*ft^3/day;
        u.qr  = stb/day;
        u.t  = day;
    otherwise
        error(['Unit ', unit, ' not supported']);
end
end

%--------------------------------------------------------------------------

function fluxNms = getResFluxNms(phns, isECLIPSE)
nph  = numel(phns);
fluxNms = repmat( struct('I', '', 'J', '', 'K', '', 'N', ''), [nph, 1] );
if isECLIPSE
    prefix  = 'FLR';
else % assume tNav
    prefix = 'FLO';
end

for phi = 1:nph
    fluxNms(phi).I = fixVarName( [prefix, phns{phi}, 'I+'] );
    fluxNms(phi).J = fixVarName( [prefix, phns{phi}, 'J+'] );
    fluxNms(phi).K = fixVarName( [prefix, phns{phi}, 'K+'] );
    fluxNms(phi).N = fixVarName( [prefix, phns{phi}, 'N+'] );
end

end

%--------------------------------------------------------------------------

function [opt, phn, unit, tr, na, isECLIPSE] = checkAndProcessInput(fn, rstrt, opt)
isECLIPSE = ismember(rstrt.INTEHEAD{1}(95), [100, 300, 500, 700]);
% Get first non-empty intehead-data for general info
ix = find(cellfun(@isempty, rstrt.INTEHEAD)==false, 1);
ih = rstrt.INTEHEAD{ix};

% number active gridcells
na  = ih(12);

% Units:
un  = {'metric', 'field', 'lab'};
ii  = ih(3);
if ~any(ii == 1:3)
    warning('Unknown unit indicator %d, assuming metric', ii);
    ii = 1;
end
un = un{ii};
unit = getUnits(un);

% Valid phase-configurations (w,o,g):
Mp = logical([0 1 0; 1 0 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1]);
ii = ih(15);
if ~any(ii == 1:7)
    warning('Unknown phase indicator %d, assuming three-phase', ii);
    ii = 7;
end
phns = {'WAT', 'OIL', 'GAS'};    % phase names
phn  = phns(logical(Mp(ii,:)));  % current phase names

% Get times
tr = cellfun(@(x)x(1), rstrt.DOUBHEAD).';
%try
%    tr  = convertFrom(rsspec.TIME.values, unit.t);
%catch
%    warning('Unsuccessful reading of RSSPEC-file, setting restart times to 0,1,...,(n-1).')
%    tr = (0 : (numel(rstrt.INTEHEAD)-1))';
%end


% Check that reservoir fluxes are present in output if desired
if opt.includeFluxes
    if isECLIPSE
        prefix = 'FLR';
    else
        prefix = 'FLO';
    end
    fluxesProvided = isfield(rstrt, fixVarName([prefix ,phn{1},'I+']));
    if ~fluxesProvided
        dispif(mrstVerbose, 'Reservoir fluxes not found in RSTRT-file\n');
        opt.includeFluxes = false;
    else
        % Check that GRID/EGRID and INIT files are present in case of NNC
        opt.nncPresent = isfield(rstrt, fixVarName(['FLR',phns{1},'N+']));
        if opt.nncPresent
            req = {'EGRID', 'GRID', 'INIT'};
            chk = cellfun(@(x)~isempty(dir([fn, '.', x])), req);
            if ~and(or(chk(1),chk(2)), chk(3))
                warning('Unable to process NNCs, could not find EGRID/GRID and/or INIT-files. No flux fields will be added to states.');
                opt.includeFluxes = false;
            end
        end
    end
end

% only compute mobilities if fluxes included
opt.includeMobilities = opt.includeMobilities && opt.includeFluxes;

% Check that wellSol-options are compatible with output-files
if opt.includeWellSols
    if opt.wellSolsFromRestart
        req = {'IWEL', 'SWEL', 'XWEL', 'ZWEL', 'ICON', 'SCON', 'XCON'};
        chk = isfield(rstrt, req);
        if ~all( chk )
            warning(['Unable to produce wellSols from RSTRT-file. Could not find fields: ', ...
                repmat('\n%s ', [1, nnz(~chk)]) ], req{~chk});
            opt.includeWellSols = false;
        end
    else % wellSols from summary
        req = {'SMSPEC', 'UNSMRY'};
        chk = cellfun(@(x)~isempty(dir([fn, '.', x])), req);
        if ~all(chk)
            warning(['Unable to produce wellSols from summary. Could not find files: ', ...
                repmat(['\n', fn, '.%s '], [1, nnz(~chk)]) ], req{~chk});
            opt.includeWellSols = false;
        end
    end
end
if ~isempty(opt.additionalFields)
    warning('Including addtional fields is not yet implemented')
end
end

%--------------------------------------------------------------------------

function ii = ms2rs(tm, tr)
% assumes numel(tm) >= numel(tr) and that tr is 'almost' subset of tm
% set threshold to 1/10 of smallest time-step
thr = min(diff(tm))/10;
ii = false(size(tm));
ir_next = 1;
for k = 1:numel(ii)
    if tr(ir_next) <= tm(k) + thr
        ii(k) = true;
        ir_next = ir_next + 1;
    end
end
assert(numel(tr)==nnz(ii), 'Ministeps not compatible with report times')
if max(abs(tr-tm(ii))) > thr
    warning('Mismatch observed between ministeps and reportsteps')
end
end

%--------------------------------------------------------------------------

function states = processWellStates(states, opt)
ns = numel(states);
nw = numel(states{1}.wellSol);
if opt.splitWellsOnSignChange
    ws = cellfun(@(x)x.wellSol, states, 'UniformOutput', false);
    sgn = cellfun(@(x)[x.sign], ws, 'UniformOutput', false);
    sgn = vertcat(sgn{:});
    for k = 1:nw
        if any(sgn(:,k)~= sgn(1,k))
            nw = nw + 1;
            firstSgn = sgn(1,k);
            [firstNm, secondNm] = deal(strcat(ws{1}(k).name, ' (inj)'), strcat(ws{1}(k).name, ' (prod)'));
            if firstSgn < 0
                [firstNm, secondNm] = deal(secondNm, firstNm);
            end
                
            for k1 = 1:ns
                ncon = numel(ws{k1}(k).cells);
                % copy into new last entry
                ws{k1}(nw) = ws{k1}(k);
                % rename and update status
                ws{k1}(k).name   = firstNm;
                ws{k1}(k).status = ws{k1}(k).status && (ws{k1}(k).sign == firstSgn);
                ws{k1}(k).sign   = firstSgn;
                if ~ws{k1}(k).status
                    ws{k1}(k).cstatus = false(ncon,1);
                    ws{k1}(k).resv    = 0;
                    ws{k1}(k).flux    = zeros(ncon,1);
                    ws{k1}(k).cqs     = zeros(ncon,3);
                end
                
                ws{k1}(nw).name = secondNm;
                ws{k1}(nw).status = ws{k1}(nw).status && (ws{k1}(nw).sign ~= firstSgn);
                ws{k1}(nw).sign   = -firstSgn;
                if ~ws{k1}(nw).status
                    ws{k1}(nw).cstatus = false(ncon,1);
                    ws{k1}(nw).resv    = 0;
                    ws{k1}(nw).flux    = zeros(ncon,1);
                    ws{k1}(nw).cqs     = zeros(ncon,3);
                end
            end
        end
    end
    if numel(states{1}.wellSol) ~= nw % some have been added
        for k = 1:ns
            states{k}.wellSol = ws{k};
        end
    end
end

alwaysClosed = true(nw,1);
for k = 1:numel(states)
    ws = states{k}.wellSol;
    for k1 = 1:numel(ws)
        if ws(k1).status
            if opt.removeCrossflow
                cf = ws(k1).flux*ws(k1).sign < 0;
                ws(k1).flux(cf) = 0;
            end
            if ws(k1).resv*ws(k1).sign <= opt.setToClosedTol
                ncon = numel(ws(k1).cells);
                ws(k1).status  = false;
                ws(k1).cstatus = false(ncon, 1);
                ws(k1).resv    = 0;
                ws(k1).flux    = zeros(ncon, 1);
                ws(k1).cqs     = zeros(ncon, 3);
            end
            alwaysClosed = and(alwaysClosed, ~vertcat(ws.status));
        end
    end
    states{k}.wellSol = ws;
end
if opt.removeClosedWells
    for k = 1:numel(states)
        states{k}.wellSol = states{k}.wellSol(~alwaysClosed);
    end
end
end

%--------------------------------------------------------------------------

function name = fixVarName(name)
    if ~isvarname(name)
        name = regexprep(name, {'+', '-'}, {'p', 'n'});
        name = genvarname(regexprep(name, '\W', '_'));
    end
end

function state = convertCompositions(state, prefix, rstrt, step, newfield)
    if isfield(rstrt, [prefix, '1'])
        ncomp = 1;
        while ncomp < 1000
            if ~isfield(rstrt, [prefix, num2str(ncomp+1)])
                break
            end
            ncomp = ncomp + 1;
        end
        tmp = repmat(rstrt.([prefix, '1']){step}, 1, ncomp);
        for i = 2:ncomp
            tmp(:, i) = rstrt.([prefix, num2str(i)]){step};
        end
        
        state.(newfield) = tmp;
    end
end
