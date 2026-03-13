function [T, W, auxiliaries] = upscaleTransGlobal(G, W, Tf, varargin)
% Calculate upscaled transmissibilities for a coarse model
%
% SYNOPSIS:
%   [Wc, Tc, aux] = upscaleTransGlobal(Gc, Wc, Tf,'pn1', pv1, ...)
%
% PARAMETERS:
%   Gc       - coarse grid
%   Wc       - coarse well structure
%   Tf       - transmissibility for fine grid
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%          - verbose -- Default value: Verbose = mrstVerbose
%          - globalFieldCases -- sets type of global fields to be used
%               'single'    : use single field as determined by Wc (default)
%               'revolving' : All combinations with one well rate = 1 and the
%                           rest bhp =0
%               'linear'    : use linear pressure bc in each coardinate
%                           direction
%               'custom'    : use custom well-configurations and/or bcs as
%                             given by 'wells', 'bc'.
%          - handleNegative -- handle negative transmissibilities:
%               'setToZero': (default)
%               'setToNan' :
%               'ignore'   : Leave as is
%          - handleUndetermined -- handle undetermined transmissibilities
%                          (due to fluxes below threshold)
%               'setToZero': default
%               'setToNan' :
%          - fluxThreshold -- Lower relative (wrt max of fineflux) value of
%                          coarse face flux in order for transmissibility to
%                          be computed. Default = 1e-8.
%          - mobility -- Cell of either one or two vectors representing coarse
%                        and/or fine mobility field, i.e., {m} or {m1, m2}.
%                        If only fine-grid mobilities are given, then
%                        coarse are computed from pv-average. If only
%                        coarse are given, then fine is computed from
%                        prolongation. Default = {} (all mobilities are set
%                        to one);
%          - clustering -- n_wells x 1 vector of indices corresponing to
%                         well clusters. Sets 'globalFieldCases' to 'revolving'
%                         and computes global fields from setting rate=1 for one
%                         well in each cluster and the rest bhp=0. Default []
%          - wells -- custom list of well structures for global field
%                         computations. Ignored unless 'globalFieldCases'
%                         is set to custom.
%          - bc    -- custom list of boundary conds. See 'wells' above.
%          - pv    -- fine grid pore volumes.
%                     Default: G.parent.cells.voumes
%
% RETURNS:
%   Tc     - Upscaled transmissibilities
%   Wc     - Well structure with updated coarse grid well indices
%   aux    - optional info/variables
%
% SEE ALSO:
%   coarsegrid module, grid_structure

opt = struct('verbose',              mrstVerbose,     ...
             'globalFieldCases',        'single',     ...
             'handleNegative',       'setToZero',     ...
             'handleUndetermined',   'setToZero',     ...
             'fluxThreshold',               1e-8,     ...
             'mobility',                    {{}},     ...
             'match_method',          'max_flux',     ...
             'bc',                            [],     ...
             'wells',                         [],     ...
             'pv',        G.parent.cells.volumes,     ...
             'linsolve',               @mldivide);


opt = merge_options(opt, varargin{:});

[G, W, opt] = checkInput(G, W, opt);

% extract parent well structure
Wf = horzcat(W.parent);

% handle input mobilities:
opt = handleMobility(G, opt);

% compute global fields as instructed by 'globalFieldCases':
solsf = computeGlobalFields(G.parent, Wf, Tf, opt);

% project solutions to coarse grid:
sols = cell(1, numel(solsf));
for k = 1:numel(sols)
    sols{k} = fine2coarse(solsf{k}, G, W, opt);
end

% compute upscaled trans for all cases:
[Tl, TWl] = computeTransGlobal(sols, G, W, opt);

% Pick single trans/WI (match_method):
[T, TW] = pickUniqueTrans(Tl, TWl, sols, opt.match_method);

% Treat negative/undetermined trans
[T, TW] = handleNonPositiveTrans(T, TW, opt);

% update W
for k = 1:numel(W),
    W(k).WI = TW{k};
end

% Give axiliary variables as output if required
if nargout > 2
    auxiliaries = struct('trans', {Tl}, 'wi', {TWl}, 'sols', {sols});
end;
end

% ------------------------------------------------------------------------
function sols = computeGlobalFields(G, W, T, opt)
switch opt.globalFieldCases
    case 'single'
        W_list  = {W};
        bc_list = {[]};
    case 'revolving'
        nw = numel(W);
        [W_list, bc_list] = deal(cell(1,nw), cell(1,nw));
        for ci = 1:nw
            w = W;
            for wi = 1:nw
                if wi==ci
                    w(wi).type = 'rate';
                    w(wi).val  = 1;
                else
                    w(wi).type = 'bhp';
                    w(wi).val  = 0;
                end
            end
            W_list{ci} = w;
        end
    case 'linear'
        dim = G.griddim;
        [W_list, bc_list] = deal(cell(1,dim), cell(1,dim));
        [pmin, pmax] = deal(0, 1);

        dd  = [1 1 2 2 3 3];
        fac = {'West', 'East', 'North', 'South', 'Top', 'Bottom'};
        val = [pmin,   pmax,   pmin,    pmax,    pmin,  pmax];
        for k = 1:min(6, 2*dim)
            bc_list{dd(k)} = pside(bc_list{dd(k)}, G, fac{k}, val(k));
        end
    case 'custom'
        W_list  = opt.wells;
        bc_list = opt.bc;
end

state0 = initState(G, [], 0);
fluid0 = initSingleFluid('mu', 1, 'rho', 1);
sols   = cell(1, numel(W_list));
for k = 1:numel(W_list)
    sols{k} = incompTPFA(state0, G, T, fluid0, ...
                                      'wells',      W_list{k}, ...
                                      'bc',         bc_list{k}, ...
                                      'LinSolve',   opt.linsolve, ...
                                      'gravity',         [0 0 0], ...
                                      'use_trans',       true);
end
end
% ------------------------------------------------------------------------
function [Tl, TWl] = computeTransGlobal(sols, G, W, opt)
% Compute coarse transmissibilities/WI for given coarse single-phase
% solution solsc

% nf x nsols matrix of candidate transmissibilities
Tl = nan(G.faces.num, numel(sols));

% cell of candidate well indices
TWl = arrayfun(@(x)nan(numel(x.cells), numel(sols)), W, ...
                                    'UniformOutput', false);
N   = G.faces.neighbors;
mob = opt.mobility{2}; % coarse mobility
ix_int = all(N, 2);    % interior faces


for k = 1:numel(sols);
    sc = sols{k};
    minFlux = opt.fluxThreshold*max(abs(sc.flux));
    ix_f = abs(sc.flux) > minFlux;
    % interior faces:
    ix = and(ix_int, ix_f);
    if nnz(ix) > 0
        dp = sc.pressure(N(ix,:))*[-1; 1];
        m  = sum(mob(N(ix,:)), 2)/2; % use average mob
        Tl(ix, k) = - sc.flux(ix)./(m.*dp);
    end
    % exterior faces:
    ix = and(~ix_int, ix_f);
    if nnz(ix) > 0
        sgn = 2*(N(ix, 1) > 0) -1;
        c   = max(N(ix,:), [], 2);
        dp  = sc.facePressure(ix) - sc.pressure(c);
        m   = mob(c);
        Tl(ix, k) = sgn.*sc.flux(ix)./(m.*dp);
    end
    % wells
    for kw = 1: numel(W)
        c    = W(kw).cells;
        flux = sc.wellSol(kw).flux;
        dp   = sc.wellSol(kw).pressure - sc.pressure(c);
        mw   = mob(c);
        ix   = abs(flux) > minFlux;
        TWl{kw}(ix, k) = flux(ix)./(mw(ix).*dp(ix));
    end
end
end
% ------------------------------------------------------------------------
function sol = fine2coarse(solf, G, W, opt)
% project fine scale solution to coarse grid by volume averageing pressure
% and summing fluxes

% coarse grid pressure
pvf = opt.pv;
pc = accumarray(G.partition, solf.pressure.*pvf)./...
     accumarray(G.partition, pvf);

% coarse grid fluxes
fluxc = coarsenFlux(G, solf.flux);

% coarse grid face pressure:
% Since we have no half-trans, there is no unambigous choice for interior
% face pressure, pressure for exterior faces can be choosen 'freely', e.g.,
% use areal averaging
if isfield(solf, 'facePressure'),
    sfp  = solf.facePressure(G.faces.fconn); % sub face pressures
    sfa  = G.parent.faces.areas(G.faces.fconn);   % sub face areas
    % standard incantation:
    fcno = rldecode(1:G.faces.num, diff(G.faces.connPos), 2)';
    % finally coarse face pressures:
    pfc  = accumarray(fcno, sfp.*sfa)./accumarray(fcno, sfa);
    % set interior face pressures to nan to avoid using these for anything
    pfc( all(G.faces.neighbors, 2) ) = nan;
end

% coarse well fluxes (bhp remain the same)
if ~isempty(W),
    wsc = solf.wellSol;
    for i = 1 : numel(W)
        % standard ...
        pno = rldecode(1 : numel(W(i).cells), diff(W(i).fcellspos), 2).';
        wsc(i).flux = accumarray(pno, solf.wellSol(i).flux(W(i).fperf));
    end
else
    wsc = [];
end

% update solc
sol = struct('pressure',      pc,     ...
              'flux',          fluxc,         ...
              'facePressure',  pfc, ...
              'wellSol',       wsc);
end
% ------------------------------------------------------------------------
function [T, TW] = pickUniqueTrans(Tl, TWl, sols, match_method)
% check for nans/inf and replace by 0
ix  = ~isfinite(Tl);
ixT = all(ix,2); % set these to nan
Tl(ix) = 0;
% same for wells:
ix  = cellfun(@(x)~isfinite(x), TWl, 'UniformOutput', false);
ixW = cellfun(@(x)all(x,2), ix, 'UniformOutput', false);
for k = numel(TWl), TWl{k}(ix{k})=0; end

% concatenate absolute fluxes and well-fluxes
flux = cellfun(@(x)abs(x.flux), sols, 'UniformOutput', false);
flux = horzcat(flux{:});
wellflux = cell(numel(TWl),1);
for k = 1:numel(wellflux)
    wellflux{k} = cellfun(@(x)abs(x.wellSol(k).flux), sols, ...
                          'UniformOutput', false );
    wellflux{k} = horzcat(wellflux{k}{:});
end

switch match_method
    case 'max_flux'
        [~, choice] = max(flux, [], 2);
        nf = numel(choice);
        T  = Tl( (1:nf)' + (choice-1)*nf );

        TW = cell(1, numel(TWl));
        for k = 1:numel(TW)
            [~, choice] = max(wellflux{k}, [], 2);
            np    = numel(choice);
            TW{k} = TWl{k}( (1:np)' + (choice-1)*np );
        end
    case 'flux_mean'
        error('Not implemented')
end

% finally set all nondetermined to nan
T(ixT) = nan;
for k = 1:numel(TW)
    TW{k}(ixW{k}) = nan;
end
end
% ------------------------------------------------------------------------
function [T, TW] = handleNonPositiveTrans(T, TW, opt)
% handle undetermined trans:
if strcmp(opt.handleUndetermined, 'setToZero')
    ix = isnan(T);
    T(ix) = 0;
    dispif(opt.verbose & nnz(ix)>0, ...
           '%d undetermined transmissiblities were set to zero\n', nnz(ix));
    % wells
    nz = 0;
    for k = 1:numel(TW)
        ixw = isnan(TW{k});
        TW{k}(ixw) = 0;
        nz = nz + nnz(ixw);
    end
    if nz > 0
        warning('%d undetermined well indices were set to zero\n', nz);
    end
end
% handle negative trans:
hn = opt.handleNegative;
if ~strcmp(hn, 'ignore')
    setVal = 0;
    if strcmp(hn, 'setToNan'), setVal = nan; end
    ix = (T<0);
    T(ix) = setVal;
    dispif(opt.verbose & nnz(ix)>0, ...
        ['%d negative transmissiblities: ', hn, '\n'], nnz(ix));
    % wells
    nn = 0;
    for k = 1:numel(TW)
        ixw = TW{k} < 0;
        TW{k}(ixw) = setVal;
        nn = nn + nnz(ixw);
    end
    if nz > 0
        warning('%d negative well indices: ', hn,  '\n', nn);
    end
end
end
% ------------------------------------------------------------------------
function [G, W, opt] = checkInput(G, W, opt)
% check custom case
if ~isempty(or(opt.wells,opt.bc))
    if ~strcmp(opt.globalFieldCases, 'custom')
        warning('Input wells and/or bc ignored since globalFieldCases is not set to custom')
    end
    if or(numel(opt.wells)==numel(opt.bc), numel(opt.wells)*numel(opt.bc)==0)
        error('Number of elements in lists for wells and bc must match when both are non-empty')
    end
end
end

% ------------------------------------------------------------------------
function opt = handleMobility(G, opt)
mobs = opt.mobility;
if isempty(mobs)
    mobf = ones(G.parent.cells.num,1);
    mobc = ones(G.cells.num, 1);
elseif numel(mobs) == 1
    m = mobs{1};
    if numel(m) == G.parent.cells.num  % fine mob given
        mobf = m;
        % produce coarse mob by pv-averaging
        pvf  = opt.pv;
        mobc = accumarray(G.partition, m.*pvf)./...
               accumarray(G.partition, pvf);
    elseif numel(m) == G.cells.num
        mobc = m;
        % produce fine mob by prolonging
        mobf = m(G.partition);
    else
        error('Input mobility-vector should match either number of fine or coarse grid cells')
    end
elseif numel(mobs) == 2
    ix = [1,2];
    if numel(mobs{1}) <= numel(mobs{2}), ix = [2,1]; end
    [mobf, mobc] = deal(mobs{ix});
    if or(numel(mobf)~=G.parent.cells.num, numel(mobc)~=G.cells.num)
        error('Mobility vectors should match number of cells in fine and coarse grid')
    end
else
    error('Number of mobility vectors > 2')
end
opt.mobility = {mobf, mobc};
end
