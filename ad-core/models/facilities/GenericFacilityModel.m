classdef GenericFacilityModel < FacilityModel
    properties
        T = 288.15; % Metric standard conditions
        pressure = 101.325*kilo*Pascal; % Metric standard pressure
        SeparatorGroup
        outputFluxes = true;
        doPostUpdate
    end
    
    methods
        function model = GenericFacilityModel(varargin)
            % A generic facility model to go with generic reservoir models
            model@FacilityModel(varargin{:});
            model.setWellTargets = false;
        end
        
        function n = getNumberOfComponents(fm)
            n = fm.ReservoirModel.getNumberOfComponents();
        end
        
        function n = getNumberOfPhases(fm)
            n = fm.ReservoirModel.getNumberOfPhases();
        end
        
        function src = getComponentSources(facility, state)
            fp = facility.FacilityFlowDiscretization;
            map = fp.get(facility, state, 'FacilityWellMapping');
            if isempty(map.W)
                val = [];
            else
                val = fp.get(facility, state, 'ComponentTotalFlux');
            end
            src = struct('value', {val}, 'cells', map.cells);
        end
                
        function [surfaceRates, surfaceDensity] = getSurfaceRates(facility, state)
            fp = facility.FacilityFlowDiscretization;
            cflux = fp.get(facility, state, 'ComponentTotalFlux');
            map = fp.get(facility, state, 'FacilityWellMapping');
            model = facility.ReservoirModel;
            for c = 1:numel(cflux)
                % Sum over each well
                cflux{c} = map.perforationSum*cflux{c};
            end
            % We use a simple, but fast approach based on the
            % individual components' preference at different conditions
            [p, temp] = facility.getSurfaceConditions();
            surfaceDensity = fp.get(facility, state, 'InjectionSurfaceDensity');
            nph = model.getNumberOfPhases();
            surfaceRates = cell(1, nph);
            [surfaceRates{:}] = deal(zeros(numelValue(cflux{1}), 1));
            for c = 1:numel(cflux)
                component = model.Components{c};
                if ~isa(component, 'ConcentrationComponent')
                    composition = component.getPhaseCompositionSurface(model, state, p, temp);
                    for ph = 1:nph
                        if any(composition{ph})
                            ci = composition{ph}./surfaceDensity{ph};
                            surfaceRates{ph} = surfaceRates{ph} + ci.*cflux{c};
                        end
                    end
                end
            end
            isProd = ~map.isInjector;
            if ~isempty(facility.SeparatorGroup) && any(isProd)
                % We have separators for the producers. Perform a
                % potentially costly separation to figure out phase rates.
                cfluxProd = cellfun(@(x) x(isProd), cflux, 'UniformOutput', false);
                [surfaceRatesProd, surfaceDensityProd] = facility.SeparatorGroup.getSurfaceRates(model, cfluxProd);
                for ph = 1:nph
                    if ~isa(surfaceDensity{ph}, 'ADI') && isa(surfaceDensityProd{ph}, 'ADI')
                        surfaceDensity{ph} = model.AutoDiffBackend.convertToAD(surfaceDensity{ph}, surfaceDensityProd{ph});
                    end
                    surfaceRates{ph}(isProd) = surfaceRatesProd{ph};
                    surfaceDensity{ph}(isProd) = surfaceDensityProd{ph};
                end
            end
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end

        function [eqs, names, types, state] = getModelEquations(facility, state0, state, dt, drivingForces)
            model = facility.ReservoirModel;
            fp = facility.FacilityFlowDiscretization;
            map = fp.get(facility, state, 'FacilityWellMapping');
            primary_choice = facility.primaryVariableSet;
            varNone = strcmpi(primary_choice, 'none');
            if isempty(map.W) || varNone
                [eqs, names, types] = deal({});
                return
            end
            if strcmpi(primary_choice, 'explicit')
                nph = model.getNumberOfPhases();
                ph = model.getPhaseNames();
                [q_s, bh] = facility.getProps(state, 'q_s', 'bhp');
                
                [bhp] = facility.getProps(state, 'BottomHolePressure');
                [surfaceRates] = facility.getSurfaceRates(state);
                [eqs, names, types] = deal(cell(1, nph+1));
                eqs{end} = bhp - bh;
                names{end} = 'bh closure';
                types{end} = 'well';
                for i = 1:nph
                    eqs{i} = surfaceRates{i} - q_s{i};
                    names{i} = sprintf('q%s', ph(i));
                    types{i} = 'well';
                end
                return;
            end
            nph = model.getNumberOfPhases();
            wsum = map.perforationSum;
            % Get surface rate equations
            [surfaceRates, rhoS] = facility.getSurfaceRates(state);
            % Approximate scaling to get the tolerance in terms of
            % reservoir rates. Otherwise, e.g. gas phases may be very
            % differently scaled.
            rhoR = value(model.getProps(state0, 'Density'));
            rhoScale = bsxfun(@rdivide, value(rhoS), mean(rhoR, 1));
            % One equation for each phase corresponding to the volumetric
            % rate at surface conditions
            [sn, phnames] = model.getPhaseNames();
            switch lower(primary_choice)
                case 'standard'
                    q_s = facility.getProp(state, 'q_s');
                    [eqs, names, types] = deal(cell(1, nph+1));
                    for ph = 1:nph
                        eqs{ph} = (q_s{ph} - surfaceRates{ph}).*rhoScale(:, ph);
                        names{ph} = [phnames{ph}, 'Wells'];
                        types{ph} = 'perf';
                    end
                    targetRates = q_s;
                    ctrl_index = nph+1;
                case {'bhp_massfractions', 'bhp'}
                    varBHP = strcmpi(facility.primaryVariableSet, 'bhp');
                    varMF = strcmpi(facility.primaryVariableSet, 'bhp_massfractions');
                    assert(varBHP || varMF || varNone);
                    % We need to actually store the surface rates in wellSol
                    % here, since there are no corresponding primary variables
                    if varMF
                        massFractions = state.FacilityState.massfractions;
                        cnames = model.getComponentNames();
                        componentSources = facility.getProps(state, 'ComponentTotalFlux');
                        ncomp = model.getNumberOfComponents();
                        [eqs, names, types] = deal(cell(1, ncomp));
                        total = 0;
                        for i = 1:ncomp
                            componentSources{i} = map.perforationSum*componentSources{i};
                            total = total + componentSources{i};
                        end
                        for i = 1:ncomp-1
                            eqs{i+1} = massFractions{i} - componentSources{i}./total;
                            names{i+1} = ['well_', cnames{i}];
                            types{i+1} = 'wellcomposition';
                        end
                    elseif varBHP
                        % Just BHP
                        [eqs, names, types] = deal(cell(1, 1));
                    end
                    ctrl_index = 1;
                    targetRates = surfaceRates;
                    qSurf = value(surfaceRates);
                    for ph = 1:numel(sn)
                        fld = ['q', sn(ph), 's'];
                        for i = 1:numel(map.active)
                            ix = map.active(i);
                            state.wellSol(ix).(fld) = qSurf(i, ph);
                        end
                    end
                    clear qSurf;
                otherwise
                    error('Unknown variable set for facility %s', facility.primaryVariableSet);
            end
            % Set up AD for control equations
            nact = numel(map.active);
            bhp = facility.getProp(state, 'bhp');
            backend = model.AutoDiffBackend;
            % Equation for well matching correct control
            ctrl_eq = backend.convertToAD(zeros(nact, 1), bhp);
            % Current control types (strings)
            well_controls = {state.wellSol(map.active).type}';
            % Current targets (numerical values)
            targets = vertcat(state.wellSol(map.active).val);

            % Handle bhp
            % Since control equation convergence tollerance is given by
            % toleranceWellRate, we scale bhp-control eqs accordingly
            is_bhp  = strcmp(well_controls, 'bhp');
            if isfinite(facility.toleranceWellRate) && isfinite(facility.toleranceWellBHP)
                eqScale = facility.toleranceWellRate/facility.toleranceWellBHP;
            else
                eqScale = 1/(day()*barsa());
            end
            ctrl_eq(is_bhp) = eqScale*(bhp(is_bhp) - targets(is_bhp));

            % Following: Different types of rate controls.
            % Surface total rates
            is_rate = strcmp(well_controls, 'rate') | strcmpi(well_controls, 'vrat');
            % Surface oil rates
            is_orat = strcmp(well_controls, 'orat');
            % Surface water rates
            is_wrat = strcmp(well_controls, 'wrat');
            % Surface gas rates
            is_grat = strcmp(well_controls, 'grat');
            % Surface liquid rates (water + oil)
            is_lrat = strcmp(well_controls, 'lrat');
            % Reservoir rates (at averaged conditions from previous step)
            is_resv = strcmp(well_controls, 'resv');
            % Reservoir rates (at current conditions for each perf.)
            is_volume = strcmp(well_controls, 'volume');

            phases = model.getPhaseNames();
            is_surface_control = false(nact, 1);
            wrates = backend.convertToAD(zeros(nact, 1), bhp);
            wrates_actual = zeros(nact, 1); % Actual rates into perforations - as doubles
            surfaceRatesValue = value(surfaceRates);
            for i = 1:nph
                switch phases(i)
                    case 'W'
                        act = is_rate | is_wrat | is_lrat;
                    case 'O'
                        act = is_rate | is_orat | is_lrat;
                    case 'G'
                        act = is_rate | is_grat;
                end
                is_surface_control(act) = true;
                wrates(act) = wrates(act) + targetRates{i}(act);
                wrates_actual(act) = wrates_actual(act) + surfaceRatesValue(act, i);
            end
            ctrl_eq(is_surface_control) = wrates(is_surface_control) - targets(is_surface_control);
            % RESV controls are special
            if any(is_resv)
                rho = cellfun(@(x) x.ControlDensity, facility.WellModels(map.active), 'UniformOutput', false);
                rho = vertcat(rho{is_resv});
                resv_rates = 0;
                for ph = 1:nph
                    resv_rates = resv_rates + q_s{ph}(is_resv).*rho(:, ph);
                end
                ctrl_eq(is_resv) = resv_rates - targets(is_resv);
            end

            % Zero surface rate conditions
            wsign = vertcat(map.W.sign);
            surface_value = value(wrates);
            
            zeroTarget = targets == 0 & (is_surface_control | is_resv);
            zeroBHP = (is_bhp & sign(surface_value) ~= wsign & wsign ~= 0 & surface_value ~= 0);
            zeroRates = zeroTarget | ... % Actual target zero
                        zeroBHP | ... % Bad flow direction
                        (wrates_actual == 0 & is_surface_control); % Rates are zero for controlling rate -> no solution
            if any(zeroRates)
                q_t = 0;
                for i = 1:nph
                    q_t = q_t + targetRates{i}(zeroRates);
                end
                ctrl_eq(zeroRates) = q_t;
            end

            if any(is_volume)
                phase_flux = fp.get(facility, state, 'PhaseFlux');
                total_flux = 0;
                for i = 1:numel(phase_flux)
                    total_flux = total_flux + phase_flux{i};
                end
                well_total_flux = wsum*total_flux;
                ctrl_eq(is_volume) = well_total_flux(is_volume) - targets(is_volume);
            end

            assert(all(is_surface_control | is_bhp | is_volume | is_resv));

            eqs{ctrl_index} = ctrl_eq;
            names{ctrl_index} = 'closureWells';
            types{ctrl_index} = 'well';
        end
        
        function state = applyWellLimits(fm, state)
            active = fm.getIndicesOfActiveWells(state.wellSol);
            for i = 1:numel(active)
                w = active(i);
                well = fm.WellModels{w};
                state.wellSol(w) = fm.applyWellLimitsWellSol(well, state.wellSol(w));
            end
        end

        function model = validateModel(model, varargin)
            if isempty(model.doPostUpdate)
                model.doPostUpdate = strcmpi(model.primaryVariableSet, 'none');
            end
            model = validateModel@FacilityModel(model, varargin{:});
        end

        function [model, state] = prepareReportstep(model, state, state0, dt, drivingForces)
            [model, state] = prepareReportstep@FacilityModel(model, state, state0, dt, drivingForces);
            [model, state] = model.updateRESVControls(state, state0);
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Update state.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.updateState`

            [state, report] = updateState@FacilityModel(model, state, problem, dx, drivingForces);
            if model.doPostUpdate
                [surfaceRates, surfaceDensity] = model.getSurfaceRates(state);
                names = model.ReservoirModel.getPhaseNames();
                map = model.getProp(state, 'FacilityWellMapping');
                bhp = vertcat(state.wellSol(map.active).bhp);
                for j = 1:numel(map.active)
                    ix = map.active(j);
                    state.wellSol(ix).bhp = bhp(j);
                    for i = 1:numel(surfaceRates)
                        state.wellSol(ix).(['q', names(i), 's']) = surfaceRates{i}(j);
                    end
                end
            end
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            % Update pressure drop
            wellSol = state.wellSol;
            actWellIx = model.getIndicesOfActiveWells(wellSol);
            nw = numel(actWellIx);
            if nw > 0
                rho = model.ReservoirModel.getProps(state, 'Density');
                rho = value(rho);
                [wc, perf2well] = model.getActiveWellCells(wellSol);
                rho = rho(wc, :);
                % Use mobility in well-cells if no connection fluxes are
                % available (typically first step for well)
                if ~isfield(wellSol(1), 'ComponentTotalFlux') || ...
                    any(arrayfun(@(x) sum(sum(abs((x.ComponentTotalFlux)))) < 1e-20, wellSol))
                    mob = model.ReservoirModel.getProps(state, 'Mobility');
                    mob = [mob{:}];
                    mob = mob(wc,:);
                else
                    [mob, mob_i] = deal([]); % not used
                end
                    
                for i = 1:nw
                    wellNo = actWellIx(i);
                    wm = model.WellModels{wellNo};
                    % Possible index error for i here - should it be
                    % wellno?
                    rho_i = rho(perf2well == wellNo, :);
                    if ~isempty(mob)
                        mob_i = mob(perf2well == wellNo, :);
                    end
                    wellSol(wellNo) = wm.updateConnectionPressureDropState(model.ReservoirModel, wellSol(wellNo), rho_i, rho_i, mob_i);
                    ctrl = wm.W.type;
                    if ~(strcmpi(ctrl, 'bhp') || strcmpi(ctrl, 'thp'))
                        wellSol(wellNo).status = wm.W.val ~= 0;
                    end
                end
            end
            state.wellSol = wellSol;
            state = model.applyWellLimits(state);
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % Generic update function for reservoir models containing wells.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.updateAfterConvergence`

            [state, report] = updateAfterConvergence@FacilityModel(model, state0, state, dt, drivingForces);
            nw = model.getNumberOfActiveWells(state.wellSol);
            if nw > 0
                map = model.FacilityFlowDiscretization.get(model, state, 'FacilityWellMapping');
                phaseq = value(model.getProp(state, 'PhaseFlux'));
                cf = model.getProp(state, 'ComponentTotalFlux');
                if iscell(cf)
                    cf = reshape(cf, 1, []);
                end
                cf = value(cf);
                nwt = numel(state.wellSol);
                active = false(nwt, 1);
                active(map.active) = true;
                actIndex = zeros(nwt, 1);
                actIndex(map.active) = (1:numel(map.active))';
                nph = model.getNumberOfPhases();
                ncomp = model.getNumberOfComponents();
                for wi = 1:nwt
                    if active(wi)
                        act = map.perf2well == actIndex(wi);
                        flux = phaseq(act, :);
                        ctf = cf(act, :);
                    else
                        nc = numel(drivingForces.W(wi).cells);
                        flux = zeros(nc, nph);
                        ctf = zeros(nc, ncomp);
                    end
                    state.wellSol(wi).flux = flux;
                    state.wellSol(wi).ComponentTotalFlux = ctf;
                end
            end
        end
        
        function [groups, names, models] = getStateFunctionGroupings(model)
            [groups, names, models] = getStateFunctionGroupings@PhysicalModel(model);
            name = 'FacilityFlowDiscretization';
            ffd = model.(name);
            if ~isempty(ffd)
                groups = [groups, {ffd}];
                if nargout > 1
                    names = [names, {name}];
                    models = [models, {model}];
                end
            end
        end
        
        function [p, T] = getSurfaceConditions(fm)
            p = fm.pressure;
            T = fm.T;
        end
        
        function wellSol = applyWellLimitsWellSol(fm, well, wellSol)
            % Update solution variables and wellSol based on the well
            % limits. If limits have been reached, this function will
            % attempt to re-initialize the values and change the controls
            % so that the next step keeps within the prescribed ranges.
            withinLimits = true;
            if ~well.allowControlSwitching
                % We cannot change controls, so we return
                return
            end
            if isfield(well.W, 'status') && ~well.W.status
                % Well is inactive
                return
            end
            lims = well.W.lims;
            model = fm.ReservoirModel;

            if ~isnumeric(lims)
                phases = model.getPhaseNames();
            
                qs_t = 0;
                for i = 1:numel(phases)
                    qs_t = qs_t + wellSol.(['q', phases(i), 's']);
                end
                bhp = wellSol.bhp;
            
                if well.isInjector()
                    % Injectors have three possible limits:
                    % bhp:  Upper limit on pressure.
                    % rate: Upper limit on total surface rate.
                    % vrat: Lower limit on total surface rate.
                    modes   = {'bhp', 'rate', 'vrat'};
                    lims = well.setMissingLimits(lims, modes, inf);
                    if ~isfinite(lims.vrat)
                        % VRAT is lower limit, switch default sign
                        lims.vrat = -inf;
                    end

                    flags = [value(bhp)  > lims.bhp, ...
                              qs_t       > lims.rate, ...
                              qs_t       < lims.vrat];
                else
                    % Producers have several possible limits:
                    % bhp:  Lower limit on pressure.
                    % orat: Lower limit on surface oil rate
                    % lrat: Lower limit on surface liquid (water + oil) rate
                    % grat: Lower limit on surface gas rate
                    % wrat: Lower limit on surface water rate
                    % vrat: Upper limit on total volumetric surface rate

                    modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};
                    lims = fm.setMissingLimits(lims, modes, -inf);
                    if ~isfinite(lims.vrat)
                        % VRAT is upper limit, switch default sign
                        lims.vrat = inf;
                    end
                    [q_w, q_o, q_g, q_sl] = deal(0);
                    if model.water
                        q_w = wellSol.qWs;
                    end
                    if model.oil
                        q_o = wellSol.qOs;
                    end
                    if model.gas
                        q_g = wellSol.qGs;
                    end
                    if isfield(wellSol, 'qSs')
                        if ~isempty(wellSol.qSs)
                            q_sl = wellSol.qSs;
                        end
                    end
                    flags = [value(bhp) < lims.bhp,  ...
                        q_o          < lims.orat, ...
                        q_w + q_o    < lims.lrat, ...
                        q_g + q_sl   < lims.grat, ...
                        q_w          < lims.wrat, ...
                        qs_t         > lims.vrat];
                end
            else
                modes = {};
                flags = false;
                assert(isempty(lims) || isinf(lims))
            end
            % limits we need to check (all others than w.type):
            chkInx = ~strcmp(wellSol.type, modes);
            vltInx = find(flags(chkInx), 1);
            if ~isempty(vltInx)
                withinLimits = false;
                modes  = modes(chkInx);
                switchMode = modes{vltInx};
                fprintf('Well %s: Control mode changed from %s to %s.\n', wellSol.name, wellSol.type, switchMode);
                wellSol.type = switchMode;
                wellSol.val  = lims.(switchMode);
            end

            if ~withinLimits
                v  = wellSol.val;
                switch wellSol.type
                    case 'bhp'
                        wellSol.bhp = bhp;
                    case 'rate'
                        for ix = 1:numel(phases)
                            wellSol.(['q', phases(ix), 's']) = v*well.W.compi(ix);
                        end
                    case 'orat'
                        wellSol.qOs = v;
                    case 'wrat'
                        wellSol.qWs = v;
                    case 'grat'
                        wellSol.qGs = v;
                end % No good guess for qOs, etc...
            end
        end
        
        function lims = setMissingLimits(fm, lims, modes, val)
            missing_fields = {modes{~cellfun(@(x) isfield(lims, x), modes)}};
            for f = missing_fields
               lims = setfield(lims, f{:}, val);
            end
        end
        
        function [model, state] = updateRESVControls(model, state, state0)
            % Treat RESV
            activeWellMask = model.getWellStatusMask(state.wellSol);
            isRESVHist = cellfun(@(x) strcmpi(x.W.type, 'resv_history'), model.WellModels(activeWellMask));
            isRESVNow = cellfun(@(x) strcmpi(x.W.type, 'resv'), model.WellModels(activeWellMask));
            isRESV = isRESVHist | isRESVNow;
            if any(isRESV)
                % Local index
                W = model.getWellStruct(activeWellMask);
                
                isHist = isRESVHist(isRESV);
                compi = vertcat(W.compi);
                compi = compi(isRESV, :);
                
                rates = vertcat(W(isRESV).val);
                qs = bsxfun(@times, rates, compi);
                
                rmodel = model.ReservoirModel;
                disgas = isprop(rmodel, 'disgas') && rmodel.disgas;
                vapoil = isprop(rmodel, 'vapoil') && rmodel.vapoil;
                oix = rmodel.getPhaseIndex('O');
                gix = rmodel.getPhaseIndex('G');
                
                pvt_reg = rmodel.PVTPropertyFunctions.Density.regions;
                if isempty(pvt_reg)
                    pvt_reg = ones(rmodel.G.cells.num, 1);
                end
                cells = arrayfun(@(x) x.cells(1), W(isRESV));
                nc = numel(cells);
                regNo = pvt_reg(cells);
                
                rs = zeros(nc, 1);
                rv = zeros(nc, 1);
                pm = zeros(nc, 1);
                pv = rmodel.getProp(state0, 'PoreVolume');
                if rmodel.water
                    sw = rmodel.getProp(state0, 'sw');
                    pv = pv.*(1-sw);
                end
                for reg = 1:max(regNo)
                    subs = regNo == reg;
                    local = pvt_reg == reg;
                    pvi = pv.*local;
                    pm(subs) = sum(state0.pressure.*pvi)/sum(pvi);
                    if disgas
                        pvi = pv.*local;
                        rs(subs) = sum(state0.rs.*pvi)/sum(pvi);
                    end
                    if vapoil
                        pvi = pv.*local;
                        rv(subs) = sum(state0.rv.*pvi)/sum(pvi);
                    end
                end
                substate = struct('pressure', pm, ...
                                  's', repmat([1, 0, 0], nc, 1), ...
                                  'rs', rs, ...
                                  'rv', rv);
                rs0 = rs;
                rv0 = rv;
                if any(isHist)
                    if disgas
                        rs(isHist) = min(qs(isHist, gix)./qs(isHist, oix), rs);
                    end
                    if vapoil
                        rv(isHist) = min(qs(isHist, oix)./qs(isHist, gix), rv);
                    end
                end
                minimodel = model.ReservoirModel;
                minimodel.PVTPropertyFunctions = minimodel.PVTPropertyFunctions.subset(cells);
                minimodel.FlowPropertyFunctions = minimodel.FlowPropertyFunctions.subset(cells);
                % Avoid using flag for interpolation
                if isprop(minimodel.PVTPropertyFunctions.ShrinkageFactors, 'useSaturatedFlag')
                    minimodel.PVTPropertyFunctions.ShrinkageFactors.useSaturatedFlag = true;
                end
                b = minimodel.getProp(substate, 'ShrinkageFactors');
                shrink = 1 - rs.*rv;
                shrink0 = 1 - rs0.*rv0;
                newRates = 0;
                nph = model.getNumberOfPhases();
                factors = zeros(nc, nph);
                if rmodel.water
                    wix = 1;
                    bW = b{wix};
                    factors(:, wix) = 1./bW;
                    newRates = newRates + qs(:, wix)./bW;
                end
                if rmodel.oil
                    bO = b{oix};
                    orat = qs(:, oix);
                    factors(:, oix) = 1./(bO.*shrink0);
                    if vapoil
                        orat = orat - rv.*qs(:, gix);
                        factors(:, gix) = factors(:, gix) - rv./(bO.*shrink0);
                    end
                    newRates = newRates + orat./(bO.*shrink);
                end
                if rmodel.gas
                    bG = b{gix};
                    grat = qs(:, gix);
                    factors(:, gix) = factors(:, gix) + 1./(bG.*shrink0);
                    if vapoil
                        grat = grat - rs.*qs(:, oix);
                        factors(:, oix) = factors(:, oix) - rs./(bG.*shrink0);
                    end
                    newRates = newRates + grat./(bG.*shrink);
                end
                resvIx = find(isRESV);
                actIx = find(activeWellMask);
                for i = 1:numel(resvIx)
                    I = resvIx(i);
                    global_well_ix = actIx(I);
                    if isRESVHist(I)
                        model.WellModels{global_well_ix}.W.val = newRates(i);
                        state.wellSol(global_well_ix).val = newRates(i);
                        model.WellModels{global_well_ix}.W.type = 'resv';
                        state.wellSol(global_well_ix).type = 'resv';
                    end
                    model.WellModels{global_well_ix}.ControlDensity = factors(i, :);
                end
            end
        end
    end
end

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
