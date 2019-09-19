classdef TransportModelDG < TransportModel
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = TransportModelDG(parent, varargin)
           
            model = model@TransportModel(parent);
            model.disc = [];
            [model, discArgs] = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, discArgs{:});
            end
            model.parentModel.disc = model.disc;
            
            model.parentModel.operators = setupOperatorsDG(model.disc, model.parentModel.G, model.parentModel.rock);
            model.parentModel.outputFluxes = false;
            
            
        end
        
        %-----------------------------------------------------------------%
        function id = isDof(model, name)
            id = numel(name) > 3 && strcmp(name(end-2:end), 'dof');
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
                case {'swdof'}
                    index = model.satVarIndex('sw');
                    fn = 'sdof';
                case {'sodof'}
                    index = model.satVarIndex('so');
                    fn = 'sdof';
                case {'sgdof'}
                    index = model.satVarIndex('sg');
                    fn = 'sdof';
                case {'sdof'}
                    index = ':';
                    fn = 'sdof';
                case {'stdof'}
                    index = ':';
                    fn = 'sTdof';
                case {'rsdof'}
                    index = ':';
                    fn = 'rsdof';
                case {'rvdof'}
                    index = ':';
                    fn = 'rvdof';
                case {'pressuredof'}
                    index = ':';
                    fn = 'pressuredof';
                otherwise
                    [fn, index] = getVariableField@TransportModel(model, name, varargin{:});
            end
        end
        
        %-----------------------------------------------------------------%
        function index = satVarIndex(model, name)
            index = model.parentModel.satVarIndex(name);
        end
        
        %-----------------------------------------------------------------%
        function groupings = getStateFunctionGroupings(model)
            groupings = model.parentModel.getStateFunctionGroupings();
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            
            state.degree = repmat(model.disc.degree, model.G.cells.num, 1);
            wm = model.parentModel.FacilityModel.WellModels;
            for i = 1:numel(wm)
                state.degree(wm{i}.W.cells) = 0;
            end
            state = validateState@TransportModel(model, state);
            state = assignDofFromState(model.disc, state);
        end
        
        %-----------------------------------------------------------------%
        function [dofvars, dofnames, names, origin] = getPrimaryVariables(model, state)
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
            isParent = strcmp(origin, class(model.parentModel));
            vars = vars(isParent);
            names = names(isParent);
            dofnames = cellfun(@(bn) [bn, 'dof'], names, 'UniformOutput', false);
            origin = origin(isParent);
            dofvars = cell(1, numel(vars));
            isBO = strcmpi(origin, 'GenericBlackOilModel');
            for i = 1:numel(dofnames)
                [fn, ~] = model.getVariableField(dofnames{i}, false);
                if ~isempty(fn)
                    dofvars{i} = model.getProp(state, dofnames{i});
                elseif strcmpi(names{i}, 'x') && isBO(i)
                    if model.parentModel.water
                        sW = model.getProp(state, 'sW');
                    else
                        sW = deal(0);
                    end
                    [sG, sGdof] = model.getProps(state, 'sG', 'sGdof');
                    st = model.parentModel.getCellStatusVO(state,  1-sW-sG, sW, sG);
                    for j = 1:numel(st)
                        if numel(st{j}) == model.G.cells.num
                            st{j} = rldecode(st{j}, state.nDof, 1);
                        end
                    end
                    [rsdof, rvdof] = model.getProps(state, 'rsdof', 'rvdof');
                    xdof = st{1}.*rsdof + st{2}.*rvdof + st{3}.*sGdof;
                    dofvars{i} = xdof;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            parent = model.parentModel;
            % Get the AD state for this model
            [basevars, dofbasenames, basenames, baseorigin] = model.getPrimaryVariables(state);
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                [f, ix] = model.getVariableField(dofbasenames{i}, false);
                if strcmp(f, 'sdof') || strcmpi(dofbasenames{i}, 'xdof')
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            % Figure out saturation logic
            isP    = strcmp(dofbasenames, 'pressuredof');
            vars   = basevars;
            names  = dofbasenames;
            origin = baseorigin;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                    && sum(isS) == nph - 1;
%             assert(useTotalSaturation, 'DG currently only supports total saturation formulation!');
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sTdof';
                sTdof = model.getProp(state, replacement);
                % Replacing
                vars{isP} = sTdof;
                names{isP} = replacement;
                origin{isP} = class(model);
            else
                % Remove pressure and skip saturation closure
                vars = vars(~isP);
                names = names(~isP);
                origin = origin(~isP);
            end
            if init
                [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            end
            if useTotalSaturation
                basevars(~isP) = vars(~isP);
            else
                basevars(~isP) = vars;
            end
            % Let parent model handle state initialization
            state = model.initStateAD(state, basevars, basenames, baseorigin);
            fixedSat = false;
            for i = 1:numel(dofbasenames)
                if any(strcmpi(basenames{i}, {'sw', 'so', 'sg', 'x'}))
                    if fixedSat
                        continue
                    end
                    basenames{i} = 's';
                    dofbasenames{i}  = 'sdof';
                    fixedSat = true;
                end
                v     = model.getProp(state, dofbasenames{i});
                vm    = model.disc.getCellMean(state, value(v));
                state = model.setProp(state, basenames{i}, vm); 
            end
            if useTotalSaturation
                % Set total saturation as well
                sTdof       = vars{isP};
                state.sTdof = sTdof;
                % Evaluate at cell cubature points
                cellValue         = model.disc.evaluateProp(state, sTdof, 'cell');
                state.cellStateDG = model.setProp(state.cellStateDG, 'sT', cellValue);
                % Evaluate mean
                cellMean          = model.disc.getCellMean(state, sTdof);
                state.wellStateDG = model.setProp(state.wellStateDG, 'sT', cellMean);
                % Evaluate at face cubature points
                faceValue         = model.disc.evaluateProp(state, sTdof, 'face');
                state.faceStateDG = model.setProp(state.faceStateDG, 'sT', faceValue);
                % Set mean in state
                state = model.setProp(state, 'sT', cellMean);
            end
        end
        
        function state = initStateAD(model, state, vars, names, origin)
              
            pmodel = model.parentModel;
            isBO = strcmpi(class(pmodel), 'GenericBlackOilModel');
            if isBO && (pmodel.disgas || pmodel.vapoil)
                removed = false(size(vars));
                if pmodel.water
                    isw     = strcmpi(names, 'sw');
                    sWdof   = vars{isw};
                    sW      = model.getProp(state, 'sW');
                    removed = removed | isw;
                else
                    sW = 0;
                end

                isx  = strcmpi(names, 'x');
                xdof = vars{isx};
                sG   = model.getProps(state, 'sG');
                st   = pmodel.getCellStatusVO(state, 1-sW-sG, sW, sG);
                cells = rldecode((1:model.G.cells.num)', state.nDof, 1);
                for j = 1:numel(st)
                    if numel(st{j}) == model.G.cells.num
                        st{j} = st{j}(cells);
                    end
                end
                sGdof = st{2}.*(1-sWdof) + st{3}.*xdof;
                fill = model.disc.getFillSat(state);
                sOdof = st{1}.*(fill-sWdof) + ~st{1}.*(fill - sWdof - sGdof);
                if pmodel.water
                    sat = {sWdof, sOdof, sGdof};
                else
                    sat = {sOdof, sGdof};
                end
                removed(isx) = true;
                state = model.setProp(state, 'sdof', sat);
                names = names(~removed);
                vars  = vars(~removed);
                for i = 1:numel(names)
                    state = model.setProp(state, [names{i}, 'dof'], vars{i});
                end
                state = model.initStateFunctionContainers(state);
                
                if not(isempty(pmodel.FacilityModel))
                    % Select facility model variables and pass them off to attached
                    % class.
                    fm = class(pmodel.FacilityModel);
                    isF = strcmp(origin, fm);
                    state = pmodel.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                    removed = removed | isF;
                end
                if pmodel.disgas
                    rsSat = pmodel.getProp(state, 'RsMax');
                    rsSat = rsSat(cells);
                    rsdof = ~st{1}.*rsSat + st{1}.*xdof;
                    % rs = rs.*(value(sO) > 0);
                    state = model.setProp(state, 'rsdof', rsdof);
                end

                if pmodel.vapoil
                    rvSat = pmodel.getProp(state, 'RvMax');
                    rvSat = rvSat(cells);
                    rvdof = ~st{2}.*rvSat + st{2}.*xdof;
                    % rv = rv.*(value(sG) > 0);
                    state = model.setProp(state, 'rvdof', rvdof);
                    % No rv, no so -> zero on diagonal in matrix
                    rv = model.disc.getCellMean(state, value(rvdof));
                    sO = model.disc.getCellMean(state, value(sOdof));
                    bad_oil = value(sO) == 0 & value(rv) == 0;
                    
                    if any(bad_oil)
                        sOdof(bad_oil) = 1 - sWdof(bad_oil) - value(sGdof(bad_oil));
                        state = model.setProp(state, 'sOdof', sOdof);
                    end
                end
            else
                state = initStateAD@TransportModel(model, state, vars, names, origin);
                v     = model.getProp(state, 's');
                state = model.setProp(state, 'sdof', v);
            end
            
            state = model.evaluateBaseVariables(state);
            
        end
        
        function state = evaluateBaseVariables(model, state)
             
            [cellStateDG, faceStateDG, wellStateDG] = deal(state);
%             names = {'pressure', 's'};
            names = fieldnames(state);
            for k = 1:numel(names)
                name = names{k};
                if numel(name) > 3 && strcmp(name(end-2:end), 'dof')
                    % Get dofs
                    dof = model.getProp(state, name);
                    n = name(1:end-3);
                    % Evaluate at cell cubature points
                    cellValue = model.disc.evaluateProp(state, dof, 'cell');
                    cellStateDG = model.setProp(cellStateDG, n, cellValue);
                    % Get cell mean
                    cellMean = model.disc.getCellMean(state, dof);
                    wellStateDG = model.setProp(wellStateDG, n, cellMean);
                    % Evaluate at face cubature points
                    faceValue = model.disc.evaluateProp(state, dof, 'face');
                    faceStateDG = model.setProp(faceStateDG, n, faceValue);
                end
            end
            
            cellStateDG.sT = getTotalSaturation(cellStateDG.s);
            wellStateDG.sT = getTotalSaturation(wellStateDG.s);
            faceStateDG.sT = getTotalSaturation(faceStateDG.s);
            
            [~, ~, cells] = model.disc.getCubature((1:model.G.cells.num)', 'volume');
            [~, ~, ~, faces] = model.disc.getCubature(find(model.parentModel.operators.internalConn), 'face');
            fcells = [model.G.faces.neighbors(faces,1); model.G.faces.neighbors(faces,2)];
            
            cellStateDG.type  = 'cell';
            cellStateDG.cells = cells;
            cellStateDG.fcells = fcells;
            cellStateDG.faces = faces;
            wellStateDG.type  = 'cell';
            cellStateDG.cells = cells;
            faceStateDG.type  = 'face';
            faceStateDG.cells = fcells;
            faceStateDG.faces = faces;
            
            state.cellStateDG = cellStateDG;
            state.wellStateDG = wellStateDG;
            state.faceStateDG = faceStateDG;
            
        end
        
        function state = assignBaseVariables(model, state)
            
            names = {'s'};
            for name = names
                if isfield(state, name{1}) && isfield(state, [name{1}, 'dof'])
                    dof = model.getProp(state, [name{1}, 'dof']);
                    v   = model.disc.getCellMean(state, dof);
                    state.(name{1}) = v;
                end
            end
             
            if strcmpi(model.formulation, 'totalSaturation')
                if isfield(state, 'sT') && isfield(state, 'sTdof')
                    dof = model.getProp(state, 'stdof');
                    v   = model.disc.getCellMean(state, dof);
                    state.sT = v;
                end
            end
             
        end
        
        function model = validateModel(model, varargin)
            model = validateModel@TransportModel(model, varargin{:});
                        
            model.parentModel.FluxDiscretization = FluxDiscretizationDG(model.parentModel);
            fp = model.parentModel.FlowPropertyFunctions;
            pvt = fp.getRegionPVT(model.parentModel);
            fp = fp.setStateFunction('PoreVolume', MultipliedPoreVolumeDG(model.parentModel, pvt));
            fp = fp.setStateFunction('GravityPermeabilityGradient', GravityPermeabilityGradientDG(model.parentModel));
            model.parentModel.FlowPropertyFunctions = fp;
            
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(tmodel, state0, state, dt, drivingForces)
            state0 = tmodel.evaluateBaseVariables(state0);
            model = tmodel.parentModel;
            [acc, flux, cellflux, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            state.wellStateDG = rmfield(state.wellStateDG, 'FlowProps');
            state.wellStateDG = rmfield(state.wellStateDG, 'FluxProps');
            src = model.FacilityModel.getComponentSources(state.wellStateDG);
            % Assemble equations and add in sources
            if strcmpi(tmodel.formulation, 'missingPhase')
                % Skip the last phase! Only mass-conservative for
                % incompressible problems
                acc = acc(1:end-1);
                flux = flux(1:end-1);
                names = names(1:end-1);
                types = types(1:end-1);
            end
            d        = tmodel.disc;
            d.nDof   = state.nDof;
            d.dofPos = state.dofPos;
            psi = d.basis.psi;
            grad_psi = d.basis.grad_psi;
            ix    = d.getDofIx(state, 1, src.cells);
            cells = rldecode((1:model.G.cells.num)', d.nDof, 1);
            d.sample = acc{1}(d.getDofIx(state, Inf));
            eqs = cell(1, numel(acc));
            for i = 1:numel(acc)
                eqs{i} = d.inner(acc{i}     , psi     , 'dV') ...
                       - d.inner(cellflux{i}, grad_psi, 'dV') ...
                       + d.inner(flux{i}    , psi     , 'dS');
                if ~isempty(src.cells)
                    eqs{i}(ix) = eqs{i}(ix) - src.value{i};
                end
                if ~model.useCNVConvergence
                    pv     = model.operators.pv(cells);
                    v      = model.G.cells.volumes(cells);
                    eqs{i} = eqs{i}.*(dt./(pv.*v));
                end    
            end
        end
        
        %-----------------------------------------------------------------%
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            [model, state] = prepareTimestep@TransportModel(model, state, state0, dt, drivingForces);
            state = assignDofFromState(model.disc, state, {'pressure'});
        end
        
        %-----------------------------------------------------------------%
        function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
            vars = cellfun(@(n) n(1:end-3), vars, 'UniformOutput', false);
            [restVars, satVars, wellVars] = model.parentModel.splitPrimaryVariables(vars);
            restVars = cellfun(@(n) [n, 'dof'], restVars, 'UniformOutput', false);
            satVars = cellfun(@(n) [n, 'dof'], satVars, 'UniformOutput', false);
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Remove DG states
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
            
            if strcmpi(class(model.parentModel), 'GenericBlackOilModel') ...
                    && model.parentModel.disgas || model.parentModel.vapoil
                [state, report] = model.updateStateBO(state, problem, dx, drivingForces);
            else
                s = state;
                [restVars, satVars] = model.splitPrimaryVariables(problem.primaryVariables);
                % Update saturation dofs
                state = model.updateSaturations(state, dx, problem, satVars);
                % Update non-saturation dofs
                state = model.updateDofs(state, dx, problem, restVars);
                % Update cell averages from dofs
                state0 = state;
                state  = model.assignBaseVariables(state);
                report = [];

                if 1
                % Compute dx for cell averages
                dx0 = model.getMeanIncrement(state, state0, problem);
                % Let parent model do its thing
                problem0 = problem;
                problem0.primaryVariables = cellfun(@(n) n(1:end-3), problem0.primaryVariables, 'UniformOutput', false);
                [state0_corr, report] = updateState@TransportModel(model, state0, problem0, dx0, drivingForces);
                % Correct updates in dofs according to parent model
                dx0_corr = model.getMeanIncrement(state0_corr, state0, problem);
                
                cells    = rldecode((1:model.G.cells.num)', state.nDof, 1);
                frac     = cellfun(@(x,y) x(cells)./y(cells), dx0_corr, dx0, 'UniformOutput', false);
                for i = 1:numel(frac)
                    frac{i}(~isfinite((frac{i}))) = 0;
                end
                dx_corr  = cellfun(@(dx, f) dx.*f, dx, frac, 'UniformOutput', false);
                % Update saturation dofs
                state = model.updateSaturations(s, dx_corr, problem, satVars);
                % Update non-saturation dofs
                state = model.updateDofs(state, dx_corr, problem, restVars);
                % Update cell averages from dofs
                state = model.assignBaseVariables(state);
                end
            
            end
            
        end
        
        % --------------------------------------------------------------------%
        function [state, report] = updateStateBO(model, state, problem, dx, drivingForces)
            vars = problem.primaryVariables;
            removed = false(size(vars));
            pmodel = model.parentModel;
            cells = rldecode((1:model.G.cells.num)', state.nDof, 1);
            if pmodel.disgas || pmodel.vapoil
                % The VO model is a bit complicated, handle this part
                % explicitly.
                state0 = state;
                state = model.initStateFunctionContainers(state);

                state = pmodel.updateStateFromIncrement(state, dx, problem, 'pressure', pmodel.dpMaxRel, pmodel.dpMaxAbs);
                state = pmodel.capProperty(state, 'pressure', pmodel.minimumPressure, pmodel.maximumPressure);

                [vars, ix] = model.stripVars(vars, 'pressure');
                removed(~removed) = removed(~removed) | ix;

                % Black oil with dissolution
                [so, sg, sodof, sgdof] = model.getProps(state, 'so', 'sg', 'sodof', 'sgdof');
                if pmodel.water
                    sw = model.getProp(state, 'sw');
                    dsw = model.getIncrement(dx, problem, 'swdof');
                else
                    sw = 0;
                    dsw = 0;
                end
                % Magic status flag, see inside for doc
                st0 = pmodel.getCellStatusVO(state0, so, sw, sg);
                st = st0;
                for j = 1:numel(st)
                    if numel(st{j}) == model.G.cells.num
                        st{j} = st{j}(cells);
                    end
                end
                dr = model.getIncrement(dx, problem, 'xdof');
                % Interpretation of "gas" phase varies from cell to cell, remove
                % everything that isn't sG updates
                dsg = st{3}.*dr - st{2}.*dsw;

                if pmodel.disgas
                    rsMax = pmodel.getProp(state, 'rsMax');
                    rsMax = rsMax(cells);
                    drs_rel = rsMax.*pmodel.drsMaxRel;
                    drs = min(pmodel.drsMaxAbs, drs_rel);
                    state = model.updateStateFromIncrement(state, st{1}.*dr, problem, ...
                                                           'rsdof', inf, drs);
                end

                if pmodel.vapoil
                    rvMax = pmodel.getProp(state, 'rvMax');
                    drv_rel = rvMax.*pmodel.drsMaxRel;
                    drs = min(pmodel.drsMaxAbs, drv_rel);
                    state = model.updateStateFromIncrement(state, st{2}.*dr, problem, ...
                                                           'rvdof', inf, drs);
                end

                dso = -(dsg + dsw);
                nPh = nnz(pmodel.getActivePhases());

                ds = zeros(numel(dso), nPh);
                phIndices = pmodel.getPhaseIndices();
                if pmodel.water
                    ds(:, phIndices(1)) = dsw;
                end
                if pmodel.oil
                    ds(:, phIndices(2)) = dso;
                end
                if pmodel.gas
                    ds(:, phIndices(3)) = dsg;
                end

                state = model.updateStateFromIncrement(state, ds, problem, 'sdof', inf, pmodel.dsMaxAbs);
                state.s = model.disc.getCellMean(state, state.sdof);
                if pmodel.vapoil
                    state.rs = model.disc.getCellMean(state, state.rsdof);
                end
                if pmodel.disgas
                    state.rv = model.disc.getCellMean(state, state.rvdof);
                end
                kr = pmodel.FlowPropertyFunctions.RelativePermeability;
                state = kr.applyImmobileChop(model, state, state0);

                % We should *NOT* be solving for oil saturation for this to make sense
                assert(~any(strcmpi(vars, 'sodof')));
                
                problem0 = problem;
                problem0.primaryVariables = {'swdof', 'sodof', 'sgdof'};
                if pmodel.disgas
                    problem0.primaryVariables{end+1} = 'rsdof';
                end
                if pmodel.vapoil
                    problem0.primaryVariables{end+1} = 'rvdof';
                end
                dx0 = model.getMeanIncrement(state, state0, problem0);
                
                state_corr = computeFlashBlackOil(state, state0, pmodel, st0);
                state_corr.s  = bsxfun(@rdivide, state_corr.s, sum(state_corr.s, 2));

                dx0_corr = model.getMeanIncrement(state_corr, state0, problem0);
                
                frac = cell(1, numel(dx0));
                dx_corr = cell(1, 3);
                dx{1} = dsw; dx{2} = dso; dx{3} = dsg;
                if pmodel.disgas
                    dx{end+1} = st{1}.*dr;
                end
                if pmodel.vapoil
                    dx{end+1} = st{2}.*dr;
                end
                for i = 1:numel(frac)
                    if  ~isempty(dx0{i})
                        f = dx0_corr{i}(cells)./dx0{i}(cells);
                        f(~isfinite(f)) = 0;
                        dx_corr{i} = dx{i}.*f;
                    end
                end
                % Update saturation dofs
                state = model.updateSaturations(state0, dx_corr, problem0, {'swdof', 'sodof', 'sgdof'});
                % Update non-saturation dofs
                restVars = {};
                if pmodel.disgas
                    restVars{end+1} = 'rsdof';
                end
                if pmodel.vapoil
                    restVars{end+1} = 'rvDof';
                end
                state = model.updateDofs(state, dx_corr, problem0, restVars);
                % Update cell averages from dofs
                state = model.assignBaseVariables(state);
                
                %  We have explicitly dealt with rs/rv properties, remove from list
                %  meant for autoupdate.
                [vars, ix] = model.stripVars(vars, {'sw', 'so', 'sg', 'rs', 'rv', 'x'});
                removed(~removed) = removed(~removed) | ix;
            end

            % We may have solved for a bunch of variables already if we had
            % disgas / vapoil enabled, so we remove these from the
            % increment and the linearized problem before passing them onto
            % the generic reservoir update function.
            problem.primaryVariables = vars;
            dx(removed) = [];

            report = [];
        end
        
        %-----------------------------------------------------------------%
        function dx = getMeanIncrement(model, state, state0, problem)
            
            vars = problem.primaryVariables;
            dx   = cell(numel(vars),1);
            for i = 1:numel(vars)
                vn = vars{i}(1:end-3);
                v  = model.getProp(state, vn);
                v0 = model.getProp(state0, vn);
                dx{i} = v - v0;
            end

        end
        
        % ----------------------------------------------------------------%
        function state = updateDofs(model, state, dx, problem, dofVars)
            
            for i = 1:numel(dofVars)
                state = updateStateFromIncrement(model, state, dx{i}, problem, dofVars{i}, inf, inf);
            end
            
        end
        
        % ----------------------------------------------------------------%
        function state = updateSaturations(model, state, dx, problem, satVars)

            if nargin < 5
                % Get the saturation names directly from the problem
                [~, satVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            if isempty(satVars)
                % No saturations passed, nothing to do here.
                return
            end
            % Solution variables should be saturations directly, find the
            % missing link
            saturations0 = lower(model.parentModel.getSaturationVarNames);
            saturations  = cellfun(@(n) [n, 'dof'], saturations0, 'uniformOutput', false);
            fillsat = setdiff(saturations, lower(satVars));
            nFill = numel(fillsat);
            assert(nFill == 0 || nFill == 1)
            if nFill == 1
                % Fill component is whichever saturation is assumed to fill
                % up the rest of the pores. This is done by setting that
                % increment equal to the negation of all others so that
                % sum(s) == 0 at end of update
                fillsat = fillsat{1};
                solvedFor = ~strcmpi(saturations, fillsat);
            else
                % All saturations are primary variables. Sum of saturations is
                % assumed to be enforced from the equation setup
                solvedFor = true(numel(saturations), 1);
            end
            ds = zeros(sum(state.nDof), numel(saturations));
            
            tmp = 0;
            ix = model.disc.getDofIx(state, Inf);
            for phNo = 1:numel(saturations)
                if solvedFor(phNo)
                    v = model.getIncrement(dx, problem, saturations{phNo});
                    ds(ix, phNo) = v;
                    if nFill > 0
                        % Saturations added for active variables must be subtracted
                        % from the last phase
                        tmp = tmp - v;
                    end
                end
            end
            ds(ix, ~solvedFor) = tmp;
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
            state   = model.updateStateFromIncrement(state, ds, problem, 'sdof', Inf, Inf);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TransportModel(model, state0, state, dt, drivingForces);
            state = rmfield(state, 'cellStateDG');
            state = rmfield(state, 'faceStateDG');
            state = rmfield(state, 'wellStateDG');
            
            propfn = model.parentModel.getStateFunctionGroupings();
            d = model.disc;
            d.nDof = state.nDof;
            d.dofPos = state.dofPos;
            ix = d.getDofIx(state, 1, Inf);
            psi    = model.disc.basis.psi(1);
            d.sample = state.sdof(:,1);
            for i = 1:numel(propfn)
                p = propfn{i};
                struct_name = p.getStateFunctionContainerName();
                names = p.getNamesOfStateFunctions();
                if isfield(state, struct_name)
                    for j = 1:numel(names)
                        name = names{j};
                        if ~isempty(state.(struct_name).(name))
                            v = state.(struct_name).(name);
                            nph = numel(v);
                            for ph = 1:nph
                                v{ph} = d.inner(v{ph}, psi, 'dV');
                                v{ph} = v{ph}(ix);
                            end
                            state.(struct_name).(name) = v;
                        end
                    end
                end
            end
            
            if 1
                 state = model.disc.limiter(state, 's');
            else
                d = model.disc;
                d.jumpTolerance = 1e-3;
                d.jumpLimiter = 'tvb';
                d.plotLimiterProgress = true;
                w = WENOUpwindDiscretization(model.parentModel, model.G.griddim);
                [C, pts, cells, basis, supports, linear_weights, scaling] = w.getTriangulation(model.parentModel);

                interp_setup.tri_cells = cells;
                interp_setup.tri_basis = basis;
                interp_setup.tri_points = pts;
                interp_setup.linear_weights = linear_weights;
                interp_setup.cell_support = supports;
                interp_setup.scaling = scaling;
                interp_setup.C = C;
                interp_setup.cell_support_count = cellfun(@numel, interp_setup.cell_support);

                d.interp_setup = interp_setup;
                state = d.limiters(model, state, state, false);
            end

        end
        
    end
    
end

function sT = getTotalSaturation(s)
    if iscell(s)
        sT  = 0;
        nph = numel(s);
        for i = 1:nph
            sT = sT + s{i};
        end
    else
        sT = sum(s,2);
    end
end