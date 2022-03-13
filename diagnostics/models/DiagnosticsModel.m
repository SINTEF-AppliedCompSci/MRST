classdef DiagnosticsModel < WrapperModel
    properties
        pressureNonLinearSolver
        pressureTol = 1e-4;
        %ctRatio         = 0;      % Compressibility over time ratio (weakly compressible system)
        state0 = [];
        computeForward  = true;
        computeBackward = true;
        tracerWells     = 'none';  % 'all', 'none' or list of well-names
        computeWellTOFs = false;
        maxTOF          = 1000*year;
        processCycles   = false;
        useMaxTOF       = true;
        dSolver         = [];
        pSolver         = [];
        phaseWeights    = [];
        pvScale         = 1;
    end
    
    methods
        function model = DiagnosticsModel(parent, varargin)
            opt = struct('convertFluid', false, 'state0', []);
            assert(isa(parent, 'PressureModel'));
            model = model@WrapperModel(parent);
            %model.useIncTol = false;
            model.pressureTol = 1e-4;
            model = merge_options(model, varargin{:});
            if isempty(model.pressureNonLinearSolver)
                model.pressureNonLinearSolver = NonLinearSolver();
            end
            if isempty(model.phaseWeights)
                model.phaseWeights = model.state0.s;
            end
        end
        % -----------------------------------------------------------------
        
        function [state, names, origin] = getStateAD(model, state, varargin)
            state = model.reduceState(state, true);
            opt = struct('initPressure', false, 'initTOFForward', false, 'initTOFBackward', false);
            opt = merge_options(opt, varargin{:});
            state = model.validateState(state);
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
            % place 'generic' TOF last
            vars   = [vars, {state.tof(:,1), state.tof(:,2)}];
            names  = [names, {'tof_forward', 'tof_backward'}];
            origin = [origin, {'DiagnosticsModel', 'DiagnosticsModel'}];
            % don't keep saturations
            isP = strcmp(names, 'pressure');
            origP = origin{isP};
            keep = isP | cellfun(@(x) ~strcmp(x, origP), origin);
            % select which variables to initialize
            init = keep;
            if ~opt.initPressure
                init(1:end-2) = false;
            end
            if ~opt.initTOFForward
                init(end-1) = false;
            end
            if ~opt.initTOFBackward
                init(end) = false;
            end
            if any(init)
                [vars{init}] = model.parentModel.AutoDiffBackend.initVariablesAD(vars{init});
             end
            state = model.initStateAD(state, vars, names, origin);
            names = names(init);
            origin = origin(init);
        end
        % -----------------------------------------------------------------
        
        function state = validateState(model, state)
            state = model.parentModel.validateState(state);
            nc = model.G.cells.num;
            if ~isfield(state, 'tof')
                state.tof = zeros(nc,2);
            else
                sz = size(state.tof);
                assert(all(sz==[nc,2]), 'Unexpeted size of state.tof');
            end
        end
        % -----------------------------------------------------------------
        
        function state = initStateAD(model, state, vars, names, origin)
            ix = strcmp('DiagnosticsModel', origin);
            for i = 1:numel(names)
                if ix(i)
                    state = model.setProp(state, names{i}, vars{i});
                end
            end
            state = model.parentModel.initStateAD(state, vars(~ix), names(~ix), origin(~ix));
        end
        % -----------------------------------------------------------------
        
        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
                case 'tof_forward'
                    fn = 'tof';
                    index = 1;
                case 'tof_backward'
                    fn = 'tof';
                    index = 2;
                otherwise
                    [fn, index] = model.parentModel.getVariableField(model, name, varargin{:});
            end
        end
        % -----------------------------------------------------------------
        
        function system = getCoupledSystem(model, state, forces, varargin)
            opt = struct('resOnly',   false, ...
                         'validate',  false);
            opt = merge_options(opt, varargin{:});
            resOnly = opt.resOnly;
            
            if opt.validate
                model  = model.validateModel(forces);
            else
                model.parentModel.parentModel.FacilityModel = model.parentModel.parentModel.FacilityModel.validateModel(forces);
            end
            
            system = struct('pressure', [], 'forward', [], 'backward', []);
            % pressure system
            assert(~isempty(model.state0))
            state = model.reduceState(state, true);
            state = model.validateState(state);
            [fstate, var_names, origin] = getStateAD(model, state, 'initPressure', ~resOnly);
            [eqs, eq_names, eq_types] = model.parentModel.getModelEquations(model.state0, fstate, 1*day, forces);
            system.pressure = struct('eqs', {eqs}, 'eqNames', {eq_names}, 'eqTypes', {eq_types}, ...
                'varNames', {var_names}, 'varOrig', {origin});
            
            % forward TOF system
            if model.computeForward
                [fstate, var_names, origin] = getStateAD(model, state, 'initPressure', ~resOnly, 'initTOFForward', ~resOnly);
                [eqs, eq_names, eq_types, src] = model.getTOFEquations(fstate, 'forward');
                system.forward = struct('eqs', {eqs}, 'src', {src}, 'eqNames', {eq_names}, 'eqTypes', {eq_types}, ...
                    'varNames', {var_names}, 'varOrig', {origin});
            end
            
            % backward TOF system
            if model.computeBackward
                [fstate, var_names, origin] = getStateAD(model, state, 'initPressure', ~resOnly, 'initTOFBackward', ~resOnly);
                [eqs, eq_names, eq_types, src] = model.getTOFEquations(fstate, 'backward');
                system.backward = struct('eqs', {eqs}, 'src', {src}, 'eqNames', {eq_names}, 'eqTypes', {eq_types}, ...
                    'varNames', {var_names}, 'varOrig', {origin});
            end
        end
        % -----------------------------------------------------------------
        
        function psys = getPressureSystem(model, state, forces, resOnly)
            if nargin < 4
                resOnly = false;
            end
            [fstate, var_names, origin] = getStateAD(model, state, 'initPressure', ~resOnly);
            [eqs, eq_names, eq_types] = model.parentModel.getModelEquations(model.state0, fstate, 1*day, forces);
            psys = struct('eqs', {eqs}, 'eqNames', {eq_names}, 'eqTypes', {eq_types}, ...
                'varNames', {var_names}, 'varOrig', {origin});
        end
        % -----------------------------------------------------------------
        
        function r = getPressureResiduals(model, state, W)
            psys = getPressureSystem(model, state, W, true);
            r = psys.eqs;
        end
        % -----------------------------------------------------------------
        
        function r = getTOFResiduals(model, state, W, direction)
            forces = model.getValidDrivingForces();
            forces.W = W;
            model  = model.validateModel(forces);
            state = model.validateState(state);
            [fstate, var_names, origin] = model.getStateAD(state);
            r = model.getTOFEquations(fstate, forces, direction); 
        end
        % -----------------------------------------------------------------    
            
        function [tof_eqs, tof_names, tof_types, src] = getTOFEquations(model, state, direction)
            if nargin < 3
                direction = 'forward';
            end
            % reset tof to zero in case of thresholding
            
            % setup single forward or backward tof-equation to create
            % system matrix/derivatives
            gbmodel = model.parentModel.parentModel;
            tmp   = gbmodel.getProps(state, 'PhaseFlux');
            flux = tmp{1};
            for k = 2:numel(tmp)
                flux = flux + tmp{k};
            end
            tmp = gbmodel.FacilityModel.getProps(state, 'PhaseFlux');
            wellflux = tmp{1};
            for k = 2:numel(tmp)
                wellflux = wellflux + tmp{k};
            end
            
            [wc, p2w] = gbmodel.FacilityModel.getActiveWellCells(state.wellSol);
            W    = gbmodel.FacilityModel.getWellStruct();
            wsgn = vertcat(W.sign);
            psgn = wsgn(p2w);
            
            perfInj = (wellflux > 0);
            [Div, Up] = deal(gbmodel.operators.Div, gbmodel.operators.faceUpstr);
            pv = (sum(model.phaseWeights, 2).*gbmodel.operators.pv)./model.pvScale;
            % for now, only include non-crossflow connections
            if strcmp(direction, 'forward')
                tof = model.getProp(state, 'tof_forward');
                injCells  = perfInj & (psgn>0);
                src       = struct('val', wellflux(injCells), 'cells', wc(injCells));
                dd = -Div(flux);
                dd(src.cells) = dd(src.cells) + src.val;
                tof_eqs{1} = Div(flux.*Up(flux>0, tof)) + dd.*tof - pv;
                tof_names = {'tof_forward'};
                tof_types = {'cell'};
            else
                tof = model.getProp(state, 'tof_backward');
                prodCells = ~perfInj & (psgn<0);
                src       = struct('val', wellflux(prodCells), 'cells', wc(prodCells));
                dd = Div(flux);
                dd(src.cells) = dd(src.cells) - src.val;
                tof_eqs{1} = -Div(flux.*Up(flux<0, tof)) + dd.*tof - pv;
                if model.useMaxTOF && ~isa(tof_eqs{1}, 'double')
                    prob = getProblemCells(model, -flux, src, pv, model.maxTOF);
                    tof_eqs{1}(prob) = (pv(prob)/model.maxTOF).*(tof(prob)-model.maxTOF*ones(nnz(prob),1));
                    badCells = thresholdConnectedComponents(tof_eqs{1}.jac{end}, pv, max(-value(src.val)), model.maxTOF);
                    tof_eqs{1}(badCells) = (pv(badCells)/model.maxTOF).*(tof(badCells)-model.maxTOF*ones(nnz(badCells),1));
                end
                tof_names = {'tof_backward'};
                tof_types = {'cell'};
            end
        end
        % -----------------------------------------------------------------
        
        function [state, report] = solvePressure(model, W, varargin)
            opt = struct('state0', []);
            %dt = model.ctRatio*day;
            [opt, extra] = merge_options(opt, varargin{:});
            if isempty(opt.state0)
                opt.state0 = model.state0;
            end
             ls = BackslashSolverAD('keepNumber', model.G.cells.num);
             nls = NonLinearSolver('LinearSolver', ls, 'useLineSearch', true, 'maxIterations', 10);
            [state, report] = standaloneSolveAD(opt.state0, model.parentModel, 1*day, 'W', W, extra{:}, 'NonLinearSolver', nls);
        end
        % -----------------------------------------------------------------
        
        function [state, D] = solveDiagnostics(model, state, varargin)
            opt = struct('wells', []);
            opt = merge_options(opt, varargin{:});
            if ~isempty(opt.wells)
                forces   = model.getValidDrivingForces();
                forces.W = opt.wells;
                model    = model.validateModel(forces);
            elseif isempty(model.parentModel.parentModel.FacilityModel)
                error('Need validated model or well-input');
            end
            
            %gbmodel = model.parentModel.parentModel;
            D = getDiagnosticsStruct(model);
            state = model.validateState(state);
            if model.computeForward
                fstate = model.getStateAD(state, 'initTOFForward', true);
                [eqs, ~, ~, src] = getTOFEquations(model, fstate, 'forward');
                D = computeDiagnosticsQuantities(model, eqs, src, D, 'forward');
                state.tof(:,1) = D.tof(:,1);
            end
            if model.computeBackward
                fstate = model.getStateAD(state, 'initTOFBackward', true);
                [eqs, ~, ~, src] = getTOFEquations(model, fstate, 'backward');
                D = computeDiagnosticsQuantities(model, eqs, src, D, 'backward');
                state.tof(:,2) = D.tof(:,2);
            end
        end
        % -----------------------------------------------------------------
        
        function L = solveAdjoint(model, system, objPartials)
            % system:      from getCoupledSystem
            % objPartials: struct wrt pressure, forward and/or backward variables
            L = struct('forward', [], 'backward', []);
            nv = system.pressure.eqs{end}.numelValue;
            p_rhs = [];
            
            % adjoint of forward tof/tracer system
            if model.computeForward && ~isempty(objPartials.forward)
                [L.forward, p_rhs] = computeAdjointQuantities(model, system.forward, objPartials.forward, p_rhs);
            end
            % adjoint of backward tof/tracer system
            if model.computeBackward && ~isempty(objPartials.backward)
                [L.backward, p_rhs] = computeAdjointQuantities(model, system.backward, objPartials.backward, p_rhs);
            end
            % adjoint of pressure system
            L.pressure = computeAdjointPressure(model, system.pressure, objPartials.pressure, p_rhs);
        end
    end
end
% -----------------------------------------------------------------
% -----------------------------------------------------------------

function Lp = computeAdjointPressure(model, system, partials, p_rhs)
% remove last eq and last var
eqs = system.eqs;

% transpose and reset rhs (remember minus)
rhs = p_rhs' - horzcat(partials{:});
% fake ADI
tmp = initVariablesADI(0);
tmp.jac{1} = rhs; 
ls = BackslashSolverAD('keepNumber', model.G.cells.num);
problem = LinearizedProblem(eqs, system.eqTypes, system.eqNames, system.varNames, []);
%[Lp, r, rep] = model.pSolver.solveAdjointProblem([], problem, [], tmp, model);
[Lp, r, rep] = ls.solveAdjointProblem([], problem, [], tmp, model);
%Lp = model.pSolver.solveLinearProblem(problem, model);
end
% -----------------------------------------------------------------

function [Ld, p_rhs] = computeAdjointQuantities(model, system, partials, p_rhs)
% for now, last variable is tof
[eqs, src] = deal(system.eqs, system.src);
nc = model.G.cells.num;
A = eqs{end}.jac{end}.';
% adjoint tof-eq: A*lt = -DJ/Dtau
% adjoint tracer eq: A'*lc = -DJ/Dci
doTOF    = ~isempty(partials.tof);
doTracer = ~isempty(partials.tracer);

rhs = -full([partials.tof, partials.tracer]).';
solver = model.dSolver;
if isempty(solver)
    X = A\rhs;
else
    X = zeros(size(rhs));
    for k = size(X,2)
        X(:,k) = solver(A, rhs(:,k));
    end
end
if doTOF
    Ld.tof = X(:,1);
    ix = 2;
else
    Ld.tof = [];
    ix = 1;
end
if doTracer
    Ld.tracer = X(:, ix:end);
else
    Ld.tracer = [];
end
% get rhs-contribution to pressure system
% tof-contribution:  [Dtof_eq/D(p,q,..)]' * lt
% tracer-contribution:
%    [Dtof_eq/D(p,q,..)]' * sum(lci) - [Dsrc/D(p,q,..)]' * sum(Di*lci)
% where Di picks the source cells for tracer i
J = eqs{1}.jac;
B = horzcat(J{1:end-1})';
if isempty(p_rhs)
    p_rhs = zeros(size(B,1),1);
end
if doTOF
    p_rhs = p_rhs - B*Ld.tof;
end
if doTracer
    Js = src.val.jac;
    C = horzcat(Js{1:2})';
    tmp = zeros(nc,1);
    for trNo = 1:size(Ld.tracer, 2)
        wc = src.wellCells{trNo};
        tmp(wc) = tmp(wc) + Ld.tracer(wc, trNo);
    end
    p_rhs = p_rhs - ( B*sum(Ld.tracer, 2) - C*tmp );
end
end
% -----------------------------------------------------------------

function  D = computeDiagnosticsQuantities(model, eqs, src, D, direction)
if strcmp(direction, 'forward')
    wells = D.inj;
    ntr   = size(D.itracer, 2);
    nwtof = 0;
    if isfield(D, 'itof')
        nwtof = size(D.itof, 2);
    end
    [tofix, trnm, tofnm] = deal(1, 'itracer', 'itof');
else
    wells = D.prod;
    ntr   = size(D.ptracer, 2);
    nwtof = 0;
    if isfield(D, 'ptof')
        nwtof = size(D.ptof, 2);
    end
    [tofix, trnm, tofnm] = deal(2, 'ptracer', 'ptof');
end
A = eqs{1}.jac{1};
gboModel = model.parentModel.parentModel;
op = gboModel.operators;
pv = (sum(model.phaseWeights, 2).*op.pv)./model.pvScale;
nc = size(D.tof,1);

q = zeros(nc,1);
q(src.cells) = src.val;
TrRHS = zeros(nc, ntr);
wc = gboModel.FacilityModel.getWellCells();
for k = 1:ntr
    TrRHS(wc{k},k) = q(wc{k});
end
solver = [];

if isempty(solver)
    T  = A \ [pv TrRHS];
    D.tof(:,tofix) = T(:, 1);
    D.(trnm)  = T(:, 2:end);
else % if other solver, iterate over RHSs
    T = zeros(size(TrRHS)+[0, 1]);
    D.tof(:,tofix) = solver(A, pv);
    for k = 1:ntr
        D.(trnm)(:,k) = solver(A, TrRHS(:, k));
    end
end

if nwtof > 0
    C   = T(:, 2:end);
    pvi = bsxfun(@times, C, pv);
    pvi(pvi<0) = 0;
    if isempty(solver)
        X  = A \ pvi;
    else % if other solver, iterate over RHSs
        X = zeros(size(pvi));
        for k = 1:size(X, 2)
            X(:, k) = solver(A, pvi(:, k));
        end
    end
    X(X<0) = 0;
    % disregard tracer values < sqrt(eps)
    ix     = and(pvi*opt.maxTOF > X, C > sqrt(eps));
    X(~ix) = opt.maxTOF;
    X(ix)  = X(ix)./C(ix);
    D.(tofnm) = X;
end
end
% -----------------------------------------------------------------

function D = getDiagnosticsStruct(model)
gboModel = model.parentModel.parentModel;
W = gboModel.FacilityModel.getWellStruct();
iwells = vertcat(W.sign) > 0;
if ~strcmp(model.tracerWells, 'none')
    error('Handle well subsets here');
end
sub = false;
[D.inj, D.prod] = deal(find(iwells & sub), find(~iwells & sub));
nc = model.G.cells.num;
D.tof     = nan(nc,2);
D.itracer = nan(nc, numel(D.inj));
D.ptracer = nan(nc, numel(D.prod));
if model.computeWellTOFs
    D.itof = nan(nc, numel(D.inj));
    D.ptof = nan(nc, numel(D.prod));
end
end
% -----------------------------------------------------------------

function L = getLambdaStruct(model, dep)
nc = model.G.cells.num;
L = struct('tof_forward', dep.tof(1), 'tof_backward', dep.tof(2),'inj', dep.inj, 'prod', dep.prod);
L.itracer = nan(nc, numel(L.inj));
L.ptracer = nan(nc, numel(nProd));
end
% -----------------------------------------------------------------

function c = getProblemCells(model, flux, src, pv, maxTOF)
[flux, src.val] = deal(value(flux), value(src.val));
N  = model.parentModel.parentModel.operators.N;
nc = model.G.cells.num;
pos = flux >= 0;
influx = accumarray(N(:,1), -flux.*(~pos), [nc, 1]) + ...
         accumarray(N(:,2),  flux.*pos,    [nc, 1]);
influx(src.cells) = influx(src.cells) + src.val;
tm     = pv./influx;
tm(~isfinite(tm)) = inf;
c = tm > maxTOF;
end
% -----------------------------------------------------------------

function [badCells] = thresholdConnectedComponents(A, pv, maxIn, maxTOF)
    % Find strongly connected components in flux-matrix:
    [p,r,r]= dmperm(A); %#ok
    % Pick components containing more than a single cell
    ix = find(diff(r)>1);
    if ~isempty(ix)
        nc = numel(p);
        % Retrieve cell-indices to components, and construct sparse index-mapping
        c  = arrayfun(@(b,e)p(b:e)', r(ix), r(ix+1)-1, 'UniformOutput', false);
        rc = rldecode( (1:numel(c))', cellfun(@numel, c)');
        C  = sparse(vertcat(c{:}), rc, 1, nc, numel(c));
        % Compute influx to each component
        q_in = full(diag(C'*A*C));
        % Threshold
        compAboveMax = full((C'*pv)./ (max(q_in, eps*maxIn)) ) > maxTOF;

        if any(compAboveMax)
            badCells = (vertcat(c{compAboveMax}));
            dispif(mrstVerbose, 'Found %d strongly connected components, ', nnz(compAboveMax));
            dispif(mrstVerbose, 'total of %d cells, with influx below threshold.\n', numel(badCells));
        else
            badCells = [];
        end
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
