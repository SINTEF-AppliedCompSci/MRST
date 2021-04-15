classdef SimpleWell < PhysicalModel
    % Base class implementing a single, instantaneous equilibrium well model
    %
    % SYNOPSIS:
    %   wm = SimpleWell(W)
    %
    % DESCRIPTION:
    %   Base class for wells in the AD-OO framework. The base class is also
    %   the default well implementation. For this will model, the
    %   assumptions are that the well-bore flow is rapid compared to the
    %   time-steps taken by the reservoir simulator, making instantaneous
    %   equilibrium and mixing in the well-bore a reasonable assumption.
    %
    % PARAMETERS:
    %   W - Well struct. See `addWell` and `processWells`.
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to the specified value.
    %
    % RETURNS:
    %   model - Class instance of `SimpleWell`.
    %
    % SEE ALSO:
    %   `FacilityModel`, `MultisegmentWell`,

    properties
        W % Struct defining the well (see `addWell`)
        allowCrossflow % Boolean indicating if the well perforations allow cross-flow
        allowSignChange % BHP-controlled well is allowed to switch between production and injection
        allowControlSwitching % Limit reached results in well switching to another control 
        dpMaxRel % Maximum allowable relative change in well pressure
        dpMaxAbs % Maximum allowable absolute change in well pressure
        dsMaxAbs % Maximum allowable change in well composition/saturation
        VFPTable % Vertical lift table. EXPERIMENTAL.
        doUpdatePressureDrop
        simplePressureDrop
        ControlDensity
    end

    methods
        function well = SimpleWell(W, varargin)
            % Class constructor
            well = well@PhysicalModel([]);
            well.W = W;
            well.allowCrossflow = true;
            well.allowSignChange = false;
            well.allowControlSwitching = true;
            well.doUpdatePressureDrop = true;
            well.simplePressureDrop = false;

            well.dpMaxRel = inf;
            well.dpMaxAbs = inf;
            well.dsMaxAbs = inf;
            if nargin > 1
                well = merge_options(well, varargin{:});
            end
        end

        function well = updateWell(well, W)
            % Update well with a new control struct.
            %
            % SYNOPSIS:
            %   well = well.updateWell(W);
            %
            % PARAMETERS:
            %   model - Well model.
            %   W     - Well struct representing the same wells, but with
            %           changed controls, active perforations and so on.
            %
            % RETURNS:
            %   model - Updated well model.
            %
            well.W = W;
        end

        function wsol = validateWellSol(well, resmodel, wsol, state) %#ok
            % Validate wellSol for simulation
            %
            % SYNOPSIS:
            %   wellSol = well.validateWellSol(model, wellSol, resSol)
            %
            % DESCRIPTION:
            %   Validate if wellSol is suitable for simulation. This function
            %   may modify the wellSol if the errors are fixable at runtime,
            %   otherwise it should throw an error. Note that this function
            %   is analogous to validateState in the base model class.
            %
            % PARAMETERS:
            %   well    - Well model class instance.
            %   model   - `ReservoirModel` class instance.
            %   wellSol - Well solution to be updated.
            %   resSol  - Reservoir state
            %
            % RETURNS:
            %   wellSol - Updated well solution struct.

            [names, fromResModel] = well.getExtraPrimaryVariableNames(resmodel);
            for i = 1:numel(names)
                if fromResModel(i)
                    m = resmodel;
                else
                    m = well;
                end
                
                fn = m.getVariableField(names{i});
                if ~isfield(wsol, 'fn')
                    % Just set it to zero since this usually works further
                    % on down in the model. It is somewhat dangerous,
                    % though.
                    wsol = m.setProp(wsol, fn, 0);
                end
            end
        end

        function counts = getVariableCounts(wm, fld)
            % Get number of primary variables of a specific type for well
            %
            % SYNOPSIS:
            %   counts = wellmodel.getVariableCounts('bhp')
            %
            % DESCRIPTION:
            %   Should return the number of primary variables added by this
            %   well for field "fld". The simple base class only uses a
            %   single variable to represent any kind of well field, but in
            %   e.g. `MultisegmentWell`, this function may return values
            %   larger than 1.
            %
            % NOTE:
            %   A value of zero should be returned for a unknown field.
            %
            % PARAMETERS:
            %   wellmodel - Well model class instance.
            %   fld       - Primary variable name.
            %
            % RETURNS:
            %   counts - Number of variables of this type required by the
            %            well model.

            try
                fn = wm.getVariableField(fld);
            catch
                fn = [];
            end
            if isempty(fn)
                counts = 0;
            else
                counts = 1;
            end
        end

        function [names, fromResModel] = getExtraPrimaryVariableNames(well, resmodel)
            % Get additional primary variables added by this well.
            %
            % SYNOPSIS:
            %   [names, fromRes] = well.getExtraPrimaryVariableNames(model)
            %
            % DESCRIPTION:
            %   Additional primary variables in this context are variables
            %   that are not the default MRST values (surface rates for each
            %   pseudocomponent/phase and bottomhole pressure).
            %
            %   In addition, this function returns a indicator per variable
            %   if it was added by the reservoir model, or the well model.
            %
            % PARAMETERS:
            %   well     - Well class instance
            %   resmodel - `ReservoirModel` class instance.
            %
            % RETURNS:
            %   names   - Names of additional primary variables.
            %   fromRes - Boolean array indicating if the added variables
            %             originate from the well, or the reservoir.

            names = resmodel.getExtraWellPrimaryVariableNames();
            fromResModel = true(size(names));
        end

        function [names, types] = getExtraEquationNames(well, resmodel)
            % Returns the names and types of the additional equation names
            % this well model introduces.
            %
            % SYNOPSIS:
            %   [names, types] = well.getExtraEquationNames(model)
            %
            % DESCRIPTION:
            %   We have two options: Either the
            %   well itself can add additional equations (modelling e.g. flow
            %   in the well-bore) or the reservoir can add additional
            %   equations (typically for additional components)
            %
            % PARAMETERS:
            %   well     - Well class instance
            %   resmodel - `ReservoirModel` class instance.
            %
            % RETURNS:
            %   names   - Names of additional equations.
            %   types   - Type hints for the additional equations.
            %

            [names, types] = resmodel.getExtraWellEquationNames();
        end

        function [vars, names] = getExtraPrimaryVariables(well, wellSol, resmodel)
            % Returns the values and names of extra primary variables added
            % by this well.
            %
            % SYNOPSIS:
            %   [names, types] = well.getExtraEquationNames(model)
            %
            % PARAMETERS:
            %   well     - Well class instance
            %   resmodel - `ReservoirModel` class instance.
            %
            % RETURNS:
            %   vars    - Cell array of extra primary variables
            %   names   - Cell array with names of extra primary variables
            % 
            % SEE ALSO:
            %   `getExtraPrimaryVariableNames`
            [names, fromResModel] = well.getExtraPrimaryVariableNames(resmodel);
            vars = cell(size(names));
            [vars{~fromResModel}] = well.getProps(wellSol, names{~fromResModel});
            [vars{fromResModel}] = resmodel.getProps(wellSol, names{fromResModel});
        end

        function [weqs, ctrlEq, extra, extraNames, qMass, qVol, wellSol] = computeWellEquations(well, wellSol0, wellSol, resmodel, q_s, bh, packed, dt, iteration)
            % Compute well equations and well phase/pseudocomponent source terms
            [weqs, qMass, mix_s, status, cstatus, qVol] = computeWellContributionsSingleWell(well, wellSol, resmodel, q_s, bh, packed);
            ctrlEq =  setupWellControlEquationsSingleWell(well, wellSol0, wellSol, bh, q_s, status, mix_s, resmodel);

            % Update well properties which are not primary variables
            cq_sDb = qMass;
            for i = 1:numel(qMass)
                cq_sDb{i} = value(qMass{i});
            end
            cq_sDb = cell2mat(cq_sDb);

            wellSol.cqs     = bsxfun(@rdivide, cq_sDb, resmodel.getSurfaceDensities);
            wellSol.cstatus = cstatus;
            wellSol.status  = status;
            extra = {};
            extraNames = {};
        end

        function [compEqs, compSrc, compNames, wellSol] = computeComponentContributions(well, wellSol0, wellSol, resmodel, q_s, bh, packed, qMass, qVol, dt, iteration)
            % Compute component equations and component source terms
            % 
            % SEE ALSO:
            %   `ad_core.models.ReservoirModel.getExtraWellContributions`
            [compEqs, compSrc, compNames, wellSol] = resmodel.getExtraWellContributions(well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            return
        end

        function isInjector = isInjector(well)
            % Check if well is an injector
            %
            % SYNOPSIS:
            %   isInj = well.isInjector();
            %
            % PARAMETERS:
            %   well - Well class instance
            %
            % RETURNS:
            %   isInjector - boolean indicating if well is specified as an
            %                injector. Wells with sign zero is in this
            %                context defined as producers.
            isInjector = well.W.sign > 0;
        end

        function [names, types] = getWellEquationNames(well, resmodel)
            % Get the names and types for the well equations of the model.
            [~, longnames] = resmodel.getPhaseNames();
            names = cellfun(@(x) [x, 'Wells'], longnames, 'UniformOutput', false);
            types = cell(size(names));
            [types{:}] = deal('perf');
        end

        function wellSol = updateConnectionPressureDrop(well, wellSol0, wellSol, model, q_s, bhp, packed, dt, iteration)
            % Update the pressure drop within the well bore, according to a
            % hydrostatic pressure distribution from the bottom-hole to the
            % individual perforations.
            %
            % To avoid dense linear systems, this update only happens at
            % the start of each nonlinear loop.
            if iteration ~= 1 || ~well.doUpdatePressureDrop
                return
            end
            [p, mob, rho, dissolved, comp, wellvars] = unpackPerforationProperties(packed);
            for i = 1:numel(rho)
                rho{i} = value(rho{i});
                if ~isempty(dissolved)
                    for j = 1:numel(dissolved{i})
                        dissolved{i}{j} = value(dissolved{i}{j});
                    end
                end
            end
            rho     = cell2mat(rho);
            active = model.getActivePhases();
            numPh = nnz(active);

            rhoS = model.getSurfaceDensities();
            b = phaseDensitiesTobfactor(rho, rhoS, dissolved);
            w = well.W;
            if ~isfield(w, 'topo')
                nperf = numel(w.cells);
                w.topo = [(0:(nperf-1))', (1:nperf)'];
            end
            qs = wellSol.cqs; %volumetric in-flux at standard conds
            qs(qs == 0) = w.sign*1e-12;
            C = well.wb2in(w);            % mapping wb-flux to in-flux
            wbqs  = abs(C\qs);       % solve to get well-bore fluxes at surface conds
            wbqst = sum(wbqs, 2);   % total wb-flux at std conds
            % if flux is zero - just use compi
            zi = wbqst == 0;
            if any( zi )
                wbqs(zi,:)  = ones(nnz(zi),1)*w.compi;
                wbqst(zi,:) = sum(wbqs(zi,:), 2);
            end
            % Compute mixture at std conds:
            mixs = wbqs ./ (wbqst*ones(1,numPh));
            % compute volume ratio Vr/Vs
            volRat = well.compVolRat(mixs, p, b, model);
            % Mixture density at connection conds (by using static b's)
            rhoMix = (mixs*rhoS')./volRat;
            % rhoMix is now density between neighboring segments given by
            % topo(:,1)-topo(:,2) computed by using conditions in well-cell
            % topo(:,2). This is probably sufficiently accurate.

            % get dz between segment nodes and bh-node1. This is a simple
            % hydrostatic distribution.
            dpt = [0; w.dZ];
            dz  = diff(dpt);
            g   = norm(gravity);
            ddp = g*rhoMix.*dz; % p-diff between connection neighbors
            % well topology assumes we can traverse from top down, but we
            % use a loop just in case of crazy ordering. If this loop does
            % not converge, the solver will throw an error.
            cdp    = nan(size(ddp));
            cdp(1) = ddp(1);
            its = 0; maxIts = 100;
            while and(any(isnan(cdp)), its<maxIts)
                its = its +1;
                % Traverse from top node and down with the pressure
                % differentials found earlier.
                for cnr = 2:numel(cdp)
                    cdp(w.topo(cnr,2)) = cdp(w.topo(cnr,1)) + ddp(cnr);
                end
            end
            if its == maxIts
                % If this loop did not converge, something is wrong with
                % either the densities or the well itself. Regardless of
                % reason, we throw an error.
                error(['Problem with topology for well: ', wellSol.name, '. Segments appear not to be connected'])
            end
            wellSol.cdp = cdp;
        end
        
        function wellSol = updateConnectionPressureDropState(well, model, wellSol, rho_res, rho_well, mob_res)
            % Simpler version
            if ~well.doUpdatePressureDrop
                return
            end
            w = well.W;
            if ~isfield(w, 'topo')
                nperf = numel(w.cells);
                w.topo = [(0:(nperf-1))', (1:nperf)'];
            end
            if w.sign < 1 && isfield(wellSol, 'ComponentTotalFlux') && any(wellSol.ComponentTotalFlux(:))
                qVol = sum(wellSol.flux, 2);
                qMass = sum(wellSol.ComponentTotalFlux, 2);
            else
                % typically initial step, assume uniform drawdown and use compi/mobility accordingly
                sgn = w.sign;
                if sgn == 0
                    sgn = 1;
                end
                if sgn == -1 && nargin == 6 && ~isempty(mob_res)
                    qVol  = bsxfun(@times, w.WI,  mob_res);
                else
                    qVol  = w.WI*w.compi;
                end
                qMass = sum(qVol.*rho_res, 2);
                qVol = sum(qVol, 2);
            end
            if well.allowCrossflow && wellSol.sign ~= 0
                % Ignore density contribution from crossflow connections,
                % since these flows may have unphysical densities.
                xFlow = sign(qMass) ~= wellSol.sign;
                xFlow = xFlow & ~all(xFlow);
                qMass(xFlow) = 0;
                qVol(xFlow) = 0;
            end
            isShut = abs(sum(qMass)) < 1e-20 | ~wellSol.status;
            g   = norm(model.gravity);
            if well.simplePressureDrop
                if isShut
                    T = max(well.W.WI, 1e-16);
                    q_approx = bsxfun(@times, T, mob_res);
                    rhoMix = sum(q_approx.*rho_res, 2)./sum(q_approx, 2);
                else
                    rhoMix = sum(qMass)/sum(qVol);
                end
                cdp = g.*w.dZ.*rhoMix;
            else
                C = well.wb2in(w);      % mapping wb-flux to in-flux
                wbMassFlux  = abs(C\qMass);  % solve to get well-bore mass flux
                wbVolumeFlux  = abs(C\qVol); % solve to get well-bore volume flux
                % Approximate mixture density by averaging mass and volume
                % fluxes
                if isShut
                    T = max(well.W.WI, 1e-16);
                    q_approx = bsxfun(@times, T, mob_res);
                    rhoMix = sum(q_approx.*rho_res, 2)./sum(q_approx, 2);
                else
                    rhoMix = wbMassFlux./wbVolumeFlux;
                    rhoMix(isnan(rhoMix)) = 0;

                    nc = numel(rhoMix);
                    for i = (nc-1):-1:1
                        % If the density ends up being close to zero, we might be
                        % dealing with either disabled perforations or a well with
                        % zero rate. Back-traverse the mixture density and insert
                        % the values. If the values are at the end, it doesn't
                        % really matter for the pressure drop model.
                        if rhoMix(i) < 1e-8
                            rhoMix(i) = rhoMix(i+1);
                        end
                    end
                end
                bad = rhoMix == 0 | ~isfinite(rhoMix);
                if all(bad)
                    warning(['Unable to compute wellbore mixture density for', ...
                             ' well %s. Pressure drop may be incorrect.'], wellSol.name);
                else
                    rhoMix(bad) = mean(rhoMix(~bad));
                end
                % get dz between segment nodes and bh-node1. This is a simple
                % hydrostatic distribution.
                dpt = [0; w.dZ];
                dz  = diff(dpt);%.*w.cstatus;
                ddp = g*rhoMix.*dz; % p-diff between connection neighbors
                % well topology assumes we can traverse from top down, but we
                % use a loop just in case of crazy ordering. If this loop does
                % not converge, the solver will throw an error.
                cdp    = nan(size(ddp));
                cdp(1) = ddp(1);
                its = 0; maxIts = 100;
                while and(any(isnan(cdp)), its<maxIts)
                    its = its +1;
                    % Traverse from top node and down with the pressure
                    % differentials found earlier.
                    for cnr = 2:numel(cdp)
                        cdp(w.topo(cnr,2)) = cdp(w.topo(cnr,1)) + ddp(cnr);
                    end
                end
                if its == maxIts
                    % If this loop did not converge, something is wrong with
                    % either the densities or the well itself. Regardless of
                    % reason, we throw an error.
                    error(['Problem with topology for well: ', wellSol.name, '. Segments appear not to be connected'])
                end
            end
            wellSol.cdp = cdp;
        end

        function [q_s, bhp, wellSol, withinLimits] = updateLimits(well, wellSol0, wellSol, model, q_s, bhp, wellvars, p, mob, rho, dissolved, comp, dt, iteration)
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
            qs_double = zeros(size(q_s));
            for i = 1:numel(qs_double)
                qs_double(i) = value(q_s{i});
            end
            qs_t = sum(qs_double);

            actPh = model.getActivePhases();
            if isprop(model, 'solvent') && model.solvent
                solIx = sum(actPh(1:4));
            end
            gasIx = sum(actPh(1:3));
            oilIx = sum(actPh(1:2));
            watIx = 1;

            if ~isnumeric(lims)
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
                    lims = well.setMissingLimits(lims, modes, -inf);
                    if ~isfinite(lims.vrat)
                        % VRAT is upper limit, switch default sign
                        lims.vrat = inf;
                    end
                    [q_w, q_o, q_g, q_sl] = deal(0);
                    if model.water
                        q_w = qs_double(watIx);
                    end
                    if model.oil
                        q_o = qs_double(oilIx);
                    end
                    if model.gas
                        q_g = qs_double(gasIx);
                    end
                    if isprop(model, 'solvent') && model.solvent
                        q_sl = qs_double(solIx);
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
                        bhp = assignValue(bhp, v, 1);
                    case 'rate'
                        for ix = 1:numel(q_s)
                            q_s{ix} = assignValue(q_s{ix}, v*well.W.compi(ix), 1);
                        end
                    case 'orat'
                        q_s{oilIx} = assignValue(q_s{oilIx}, v, 1);
                    case 'wrat'
                        q_s{watIx} = assignValue(q_s{watIx}, v, 1);
                    case 'grat'
                        q_s{gasIx} = assignValue(q_s{gasIx}, v, 1);
                end % No good guess for qOs, etc...
            end
        end

        function wellSol = updateWellSol(well, wellSol, variables, dx, resmodel)
            % Update function for the wellSol struct
            isBHP = strcmpi(variables, 'bhp');
            if any(isBHP)
                dv = dx{isBHP};
                dv = well.limitUpdateRelative(dv, wellSol.bhp, well.dpMaxRel);
                dv = well.limitUpdateAbsolute(dv, well.dpMaxAbs);
                hasLim = isfield(well.W.lims, 'bhp');
                if hasLim
                    lim = well.W.lims.bhp;
                    dist = wellSol.bhp - lim;
                end
                wellSol.bhp = wellSol.bhp + dv;
                tol = 0.1*barsa;
                if hasLim && sign(wellSol.bhp - lim) ~= sign(dist) && abs(dist) > 2*tol
                    wellSol.bhp = well.W.lims.bhp + tol*sign(dist);
                end
                variables = variables(~isBHP);
                dx = dx(~isBHP);
            end
            isMF = strcmpi(variables, 'well_massfractions');
            if any(isMF)
                f = 'well_massfractions';
                dv = dx{isMF};
                wellSol = well.updateStateFromIncrement(wellSol, dv, [], f);
                wellSol = well.capProperty(wellSol, f, 0, 1);
                wellSol.massfractions = wellSol.massfractions./sum(wellSol.massfractions);
                variables = variables(~isMF);
                dx = dx(~isMF);
            end
            [names, fromResModel] = well.getExtraPrimaryVariableNames(resmodel);
            for i = 1:numel(dx)
                vname = variables{i};
                esub = strcmpi(names, vname);
                if any(esub) && fromResModel(esub)
                    % This variable is known by the reservoir model (e.g
                    % added polymer components)
                    wellSol = resmodel.updateStateFromIncrement(wellSol, dx{i}, [], variables{i});
                else
                    % This variable is known by the well model (e.g.
                    % multisegment pressures or fluxes)
                    wellSol = well.updateStateFromIncrement(wellSol, dx{i}, [], variables{i});
                end
            end
        end

        function [wellSol, well_shut] = updateWellSolAfterStep(well, resmodel, wellSol, wellSol0)
            % Updates the wellSol after a step is complete (book-keeping)
            well_shut = false;
            w = well.W;
            if ~w.status
                return
            end
            % Check if producers are becoming injectors and vice versa. The indexes
            % of such wells are stored in inx.
            wsg = w.sign;
            ssg = sign(well.getTotalRate(wellSol));
            closed = wsg ~= ssg && wsg ~= 0;
            % A well can be set to zero rate without being shut down. We update inx
            % to take into account this fact.
            closed = closed & ~strcmpi(w.type, 'bhp') & w.val ~= 0;
            if closed && ~well.allowSignChange && well.allowControlSwitching
                dispif(well.verbose, 'Well %s shut down.\n', w.name);
                wellSol.status = false;
                well_shut = true;
            end

            switched = ~strcmpi(wellSol.type, wellSol0.type);
            if switched && well.verbose && ...
                    ~(strcmpi(wellSol0.type, 'resv_history') && strcmpi(wellSol.type, 'resv'))
                fprintf('Step complete: Well %s was switched from %s to %s controls.\n',...
                                                                 w.name, ...
                                                                 wellSol0.type, ...
                                                                 wellSol.type);
            end
        end

        function [fn, index] = getVariableField(model, name, varargin)
            index = 1;
            switch(lower(name))
                case 'bhp'
                    fn = 'bhp';
                case 'qos'
                    fn = 'qOs';
                case 'qgs'
                    fn = 'qGs';
                case 'qws'
                    fn = 'qWs';
                case 'qss'
                    fn = 'qSs';
                case {'well_temperature', 't'}
                    fn = 'T';
                case 'well_massfractions'
                    fn = 'massfractions';
                    index = ':';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end
        end

        function ws = ensureWellSolConsistency(well, ws) %#ok
            % Run after the update step to ensure consistency of wellSol
        end
        
        function lims = setMissingLimits(well, lims, modes, val)
            missing_fields = {modes{~cellfun(@(x) isfield(lims, x), modes)}};
            for f = missing_fields
               lims = setfield(lims, f{:}, val);
            end
        end
    end
    
    methods (Static, Access = protected)
        function C = wb2in(w)
            conn = w.topo(2:end, :);
            % Number of connections between perforations
            nconn = size(conn, 1);
            % Number of perforations
            nperf = numel(w.cells);

            if nconn + 1 ~= nperf
                warning(['Mismatch between connection count (', num2str(nconn+1),...
                        ') and perforation count (', num2str(nperf), '). Well model', ...
                        'Does not appear to be a tree.']);
            end

            id = (1:nperf)';
            % First identity, then honor topology.
            ii = [id; conn(:, 1)];
            jj = [id; conn(:, 2)];
            vv = [ones(nperf, 1); -ones(nconn, 1)];

            C = sparse(ii, jj, vv, nperf, nperf);
        end

        function volRat = compVolRat(mixs, p, b, model)
        % Compute volume ratio Vr/Vs. Only complicated for blackoil.
        x = mixs;
        dg = isprop(model, 'disgas') && model.disgas;
        vo = isprop(model, 'vapoil') && model.vapoil;

        if dg || vo
            [~, isgas] = model.getVariableField('sg');
            [~, isoil] = model.getVariableField('so');

            both = find(isgas | isoil);

            g = mixs(:, isgas);
            o = mixs(:, isoil);

            if dg
                if iscell(model.fluid.rsSat)
                    rsSatFn = model.fluid.rsSat{1};
                else
                    rsSatFn = model.fluid.rsSat;
                end
                rsMax = rsSatFn(value(p));
            else
                rsMax = 0;
            end
            if isa(model, 'ThreePhaseBlackOilModel')
                % Vapoil/disgas
                if vo
                    if iscell(model.fluid.rvSat)
                        rvSatFn = model.fluid.rvSat{1};
                    else
                        rvSatFn = model.fluid.rvSat;
                    end
                    rvMax = rvSatFn(value(p));
                else
                    rvMax = 0;
                end

                gor = abs(g./o);
                gor(isnan(gor)) = inf;
                rs = min(rsMax, gor);
                ogr = abs(o./g);
                ogr(isnan(ogr)) = inf;
                rv = min(rvMax, ogr);
                d = 1-rs.*rv;
                x(:,isgas) = (x(:,isgas) - rs.*o)./d;
                x(:,isoil) = (x(:,isoil) - rv.*g)./d;
                x(:,both) = x(:,both).*(x(:,both)>0);
            else
                % Only gas dissolution
                x(:,isgas) = x(:,isgas) - rsMax.*o;
                x(:,isgas) = x(:,isgas).*(x(:,isgas)>0);
            end
        end

        volRat = sum(x./b ,2);
        end

        function qt = getTotalRate(sol)
           ns = numel(sol);
           qt = zeros([ns, 1]);
           if ns == 0
               return
           end
           typelist = {'qWs', 'qOs', 'qGs', 'qSs'};
           types    = typelist(isfield(sol(1), typelist));
           for w = 1:ns
              for t = reshape(types, 1, []),
                 x = sol(w).(t{1});
                 if ~isempty(x),
                    qt(w) = qt(w) + x;
                 end
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
