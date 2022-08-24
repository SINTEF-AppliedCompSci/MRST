classdef WellboreModel < WrapperModel
    
    properties
        W
        water, oil, gas
        InletPropertyFunctions
        thermal = false
    end
   
    methods
        
        %-----------------------------------------------------------------%
        function model = WellboreModel(rmodel, trajectories, varargin)
            
            opt = struct('names', {{}}, ...
                         'WI'   , []);
            [opt, varargin] = merge_options(opt, varargin{:});
            % Get reservoir model subgrid consisting of all perforated cells
            if (iscell(trajectories) && isnumeric(trajectories{1})) ...
                    || isstruct(trajectories)
                % Well trajectories given as 
                [W, G] = processWellboreTrajectories( ...
                            rmodel.G, rmodel.rock, trajectories, ...
                            'names', opt.names, varargin{:});
            elseif iscell(trajectories) && isfield(trajectories{1}, 'cells')
                % Well trajectories given as grids
                error('Not supported yet')
            end
            
            % Get subset of rock
            rock = extractSubrock(rmodel.rock, G.cells.global);
            % Update reservoir model operators
            rmodel = rmodel.removeStateFunctionGroupings();
            rmodel.G = G; rmodel.rock = rock;
            rmodel = rmodel.setupOperators(G, rock);
            rmodel.useCNVConvergence = false;
            rmodel.nonlinearTolerance = 1e-3;
            phases = rmodel.getPhaseNames();
            for ph = phases
                rmodel.fluid.(['kr', ph']) = @(s) s;
            end
            
            % Construct WrapperModel
            model = model@WrapperModel(rmodel);
            model.W = W;

            % Set trajectories, default refDepths and names
            refDepth = zeros(model.numWells(), 1);
            names    = cell(model.numWells,1);
            for i = 1:model.numWells()
                cells       = G.cells.wellNo == i;
                refDepth(i) = min(G.cells.centroids(cells,3));
                names{i}    = sprintf('W%d', i);
            end
            model = merge_options(model, varargin{:});
            
            % Set phase and thermal flag
            model.water = model.parentModel.water;
            model.oil   = model.parentModel.oil;
            model.gas   = model.parentModel.gas;
            if isprop(model.parentModel, 'thermal')
                model.thermal = model.parentModel.thermal;
            end
            
            % Set up operators
            model = model.setupOperators(G, rock);
            if ~isempty(opt.WI)
                model.parentModel.operators.WI = opt.WI(model.G.cells.order);
            end

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function G = computeWellGridGeometry(model, G)
            
            % Get cell segment lengths and radii and set to G.cells for
            % easy access later
            G.cells.length = model.getSegmentLength();
            G.cells.radius = model.getSegmentRadius();
            % Comute cross-sectional area of each cell segment
            area = pi.*G.cells.radius.^2;
            
            % Compute cell volumes
            G.cells.volumes = area.*G.cells.length;
            
            % Compute face area
            area = [nan; area];
            N    = G.faces.neighbors + 1;
            % Remember to rescale area-weighted normals
            G.faces.normals = G.faces.normals./G.faces.areas;
            % Compute as mean of well cross-section area of adjacent cells
            G.faces.areas   = mean(area(N),2,'omitnan');
            G.faces.normals = G.faces.normals.*G.faces.areas;
            
            % Set face segment lengths and radii for easy access later
            G.faces.length = model.getSegmentLength(true);
            G.faces.radius = model.getSegmentRadius(true);
            
        end
        
        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            
            % This will take care of parent model validation
            model = validateModel@WrapperModel(model, varargin{:});
            % Remove facility model
            model.parentModel.FacilityModel = [];
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            
            % Parent model state validation
            state = model.parentModel.validateState(state);
            % Add mass flux if not already set
            if ~isfield(state, 'massFlux')
                nf = nnz(model.parentModel.operators.internalConn);
                state.massFlux = zeros(nf, 1);
            end
            % Add bhp if not already set
            if ~isfield(state, 'bhp')
                bhp       = model.getProp(state, 'pressure');
                state.bhp = bhp(model.G.cells.topCell);
            end
            % Add surface rates if not already set
            rnames = model.getSurfaceRateNames();
            for rname = rnames
                if ~isfield(state, rname{1})
                    state.(rname{1}) = zeros(model.numWells(), 1);
                end
            end
            % Add bht if not already set
            if isprop(model.parentModel, 'thermal') && model.parentModel.thermal
                if ~isfield(state, 'bht')
                    bht       = model.getProp(state, 'temperature');
                    state.bht = bht(model.G.cells.topCell);
                end
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function rnames = getSurfaceRateNames(model)

            rnames = num2cell(model.parentModel.getPhaseNames());
            rnames = cellfun(@(name) ['q', name, 's'], rnames, 'UniformOutput', false); 

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)

            switch(lower(name))
                case {'massflux'}
                    fn = 'massFlux';
                    index = ':';
                case {'bhp'}
                    fn = 'bhp';
                    index = ':';
                case {'qws'}
                    fn = 'qWs';
                    index = ':';
                case {'qos'}
                    fn = 'qOs';
                    index = ':';
                case {'qgs'}
                    fn = 'qGs';
                    index = ':';
                case {'bht'}
                    fn = 'bht';
                    index = ':';
                case {'effect'}
                    fn = 'effect';
                    index = ':';
                otherwise
                    % Pass on to WrapperModel
                    [fn, index] = getVariableField@WrapperModel(model, name, varargin{:});
            end

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)

            % Parent model setup
            model = setupStateFunctionGroupings@WrapperModel(model, varargin{:});
            % Replace component phase flux state function
            fd = model.parentModel.FlowDiscretization;
            fd = fd.setStateFunction('ComponentPhaseFlux', WellboreComponentPhaseFlux(model));
            % Replace phase puwind flag state function
            fd = fd.setStateFunction('PhaseUpwindFlag', WellborePhaseUpwindFlag(model));
            model.parentModel.FlowDiscretization = fd;
            % Add inlet property functions
            if isempty(model.InletPropertyFunctions)
                ip = StateFunctionGrouping('InletPropertyFunctions');
                ip = ip.setStateFunction('InletPressure', WellboreInletPressure(model));
                model.InletPropertyFunctions = ip;
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [groupings, names, models] = getStateFunctionGroupings(model)

            [groupings, names, models] = model.parentModel.getStateFunctionGroupings();
            % TODO: Set linear relperms inside wellbore
            groupings = [groupings, {model.InletPropertyFunctions}];
            names     = [names    , {'InletPropertyFunctions'}    ];
            models    = [models   , class(model)                  ];

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            
            % Make wellbore grid
            GW = model.computeWellGridGeometry(G);
            % Make rock structure to represent wellbore
            rock = model.parentModel.rock;
            % Euivalent permeability from laminar pipe flow assumption
            r = vertcat(model.W.r); r = r(model.G.cells.order);            
            rock.perm = min(r.^2/8, 1e-5);
            rock.poro = ones(GW.cells.num, 1);
            model.parentModel.rock = rock;
            model.parentModel.G = GW;
            % Compute operators in standard way
            model.parentModel = model.parentModel.setupOperators(GW, rock, varargin{:});
            % Set well indices
            WI = vertcat(model.W.WI);
            WI = WI(model.G.cells.order);
            model.parentModel.operators.WI = WI;
            
            if model.thermal
                WIth = vertcat(model.W.WIth);
                WIth = WIth(model.G.cells.order);
                model.parentModel.operators.WIth = WIth;
            end
            
            % Set well grid to model
            model.G = GW;
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
            
            % Get primary variable set of parent model
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
            % Add mass flux
            vm = model.getProps(state, 'massFlux');
            vars = [vars, {vm}];
            names  = [names, 'massFlux'];
            origin = [origin, class(model)];
            % Add bhp and inlet/outlet rates
            rnames = model.getSurfaceRateNames();
            qs     = cell(size(rnames));
            [bhp, qs{:}] = model.getProps(state, 'bhp', rnames{:});
            wvars = [{bhp}, qs];
            wnames = [{'bhp'}, rnames];
            worigin = repmat({class(model)}, 1, 1 + numel(rnames));
            % Add bht if applicable
            if isprop(model.parentModel, 'thermal') && model.parentModel.thermal
                bht = model.getProp(state, 'bht');
                wvars{end+1}   = bht;
                wnames{end+1}  = 'bht';
                worigin{end+1} = class(model);
            end
            vars   = [vars, wvars];
            names  = [names, wnames];
            origin = [origin, worigin];

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)

            parent = strcmpi(origin, class(model.parentModel));
            state = model.parentModel.initStateAD(state, vars(parent), names(parent), origin(parent));

            vars   = vars(~parent);
            origin = origin(~parent); %#ok
            names  = names(~parent);
            
            for i = 1:numel(names)
                state = model.setProp(state, names{i}, vars{i});
            end
            state = model.initStateFunctionContainers(state);
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            
            [eqs, names, types, state] = model.parentModel.getModelEquations(state0, state, dt, drivingForces);
            if isempty(drivingForces.W), return; end
            Wd = drivingForces.W;
            
            if ~isfield(state, 'pu')
                state.pu = state.FluxDisc.PhaseUpwindFlag{1};
            else
                state.pu(:,end+1) = state.FluxDisc.PhaseUpwindFlag{1};
            end
            
            % Add inlet fluxes
            [Q, qs] = model.computeInletFlux(state, Wd);
            cnames = model.parentModel.getComponentNames();
            topCells = model.G.cells.topCell;
            for name = cnames
                mix = strcmpi(name{1}, names);
                eqs{mix}(topCells) = eqs{mix}(topCells) - Q{mix};
            end
            
            eix = strcmpi('energy', names);
            if any(eix)
                Qh = model.computeInletHeatFlux(state, Wd);
                eqs{eix}(topCells) = eqs{eix}(topCells) - Qh;
            end
            
            % Set up well equations, one per phase + one control equation
            nph = model.parentModel.getNumberOfPhases();
            weqs = cell(1, nph + 1);
            % Set rate equations
            rnames = model.getSurfaceRateNames();
            qsw = cell(size(rnames));
            [bhp, qsw{:}] = model.getProps(state, 'bhp', rnames{:});
            scale = (meter^3/day)^(-1);
            weqs(1:nph) = cellfun(@(qs, qph) (qs - qph).*scale, qsw, qs, 'UniformOutput', false);

            ctrl = model.AutoDiffBackend.convertToAD(zeros(model.numWells,1), bhp);
            ctrlTarget = vertcat(Wd.val);
            ctrlType = {Wd.type}';
            % Bottom-hole pressure control
            is_bhp  = strcmpi(ctrlType, 'bhp');
            % Surface total rates
            is_rate = strcmp(ctrlType, 'rate') | strcmpi(ctrlType, 'vrat');
            % Surface oil rates
            is_orat = strcmp(ctrlType, 'orat');
            % Surface water rates
            is_wrat = strcmp(ctrlType, 'wrat');
            % Surface gas rates
            is_grat = strcmp(ctrlType, 'grat');
            % Surface liquid rates (water + oil)
            is_lrat = strcmp(ctrlType, 'lrat');
            % Reservoir rates (at averaged conditions from previous step)
            is_resv = strcmp(ctrlType, 'resv');
            % Reservoir rates (at current conditions for each perf.)
            is_volume = strcmp(ctrlType, 'volume');
            % Check for unsupported well controls
            assert(~any(is_resv | is_volume), 'Not implemented yet');

            % Set BHP control equations
            ctrl(is_bhp) = ctrlTarget(is_bhp) - bhp(is_bhp);

            % Set rate control equations
            phases = model.parentModel.getPhaseNames();
            wrates = bhp*0;
            is_surface_rate = false(model.numWells(), 1);
            for ph = 1:nph
                switch phases(ph)
                    case 'W'
                        act = is_rate | is_wrat | is_lrat;
                    case 'O'
                        act = is_rate | is_orat | is_lrat;
                    case 'G'
                        act = is_rate | is_grat;
                end
                wrates(act) = wrates(act) + qsw{ph}(act);
                is_surface_rate(act) = true;
            end
            % Set BHP control equations
            ctrl(is_surface_rate) = ctrlTarget(is_surface_rate) - wrates(is_surface_rate);
            
            % Append to model equations
            weqs{end} = ctrl;
            wnames = [rnames, 'control'];
            wtypes = [repmat({'rate'}, 1, nph), 'control'];
            if any(eix)
                bht = model.getProp(state, 'bht');
                Tw  = model.getProp(state, 'T');
                Tw = Tw(topCells);
                Twinj = vertcat(Wd.T);
                flag = Qh > 0;
                Tw(flag) = Twinj(flag);
                Teq = bht - Tw;
                weqs{end+1} = Teq;
                
                wnames{end+1} = 'temperature';
                wtypes{end+1} = 'control';
            end
            
            roughness = 1e-4;
            
            d = model.G.faces.radius.*2;
            l = model.G.faces.length;
            
            % o---o---o
            [v, rho, mu, pot, flag] = model.parentModel.getProps(state, ...
                                    'ComponentPhaseFlux', ...
                                    'Density', ...
                                    'Viscosity', ...
                                    'PhasePotentialDifference', ...
                                    'PhaseUpwindFlag');
            rho{1} = model.parentModel.operators.faceAvg(rho{1});
            
            mu{1} = model.parentModel.operators.faceUpstr(flag{1}, mu{1});
            dp = wellBoreFriction(v{1}, rho{1}, mu{1}, d, l, roughness, 'massRate');
            
            weqs{end+1} = (pot{1} - dp)./mean(value(bhp));
            wnames{end+1} = 'flux';
            wtypes{end+1} = 'flux';
            
            eqs   = [eqs, weqs];
            names = [names, wnames];
            types = [types, wtypes];

            state.inletFlux = Q;
            state.inletHeatFlux = Qh;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, forces)

            parent = strcmpi(problem.types, 'cell');

            p = problem; p.primaryVariables = p.primaryVariables(parent);
            [state, report] = updateState@WrapperModel(model, state, p, dx, forces);

            for i = (nnz(parent) + 1):numel(problem.primaryVariables)
                 p = problem.primaryVariables{i};
                 % Update the state
                 state = model.updateStateFromIncrement(state, dx, problem, p);
            end

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%        
        function [model, state] = prepareReportstep(model, state, state0, dt, drivingForces) %#ok
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [Q, qph, stateInlet] = computeInletFlux(model, state, W)

            pi = model.getProps(state, 'InletPressure');
            [p, mob, rho] = model.parentModel.getProps(state, 'pressure', 'Mobility', 'Density');
            p   = p(model.G.cells.topCell);
            mob = applyFunction(@(x) x(model.G.cells.topCell), mob);
            rho = applyFunction(@(x) x(model.G.cells.topCell), rho);

            % Get top face transmissibility
            T = model.parentModel.operators.T_all(model.G.faces.topFace);

            nph   = model.parentModel.getNumberOfPhases();
            ncomp = model.parentModel.getNumberOfComponents();

            stateInlet = state;
            stateInlet.pressure(model.G.cells.topCell) = pi;
            if isprop(model.parentModel, 'thermal') && model.parentModel.thermal
                bht = model.getProp(state, 'bht');
                stateInlet.T(model.G.cells.topCell) = bht;
            end
            if ~iscell(stateInlet.s)
                stateInlet.s(model.G.cells.topCell,:) = vertcat(W.compi);
            else
                compi = vertcat(W.compi);
                for ph = 1:nph
                    stateInlet.s{ph}(model.G.cells.topCell) = compi(:,ph);
                end
            end
            if isfield(stateInlet, 'components')
                if ~iscell(stateInlet.components)
                    stateInlet.components(model.G.cells.topCell,:) = vertcat(W.components);
                else
                    components = vertcat(W.components);
                    for c = 1:ncomp
                        stateInlet.components{c}(model.G.cells.topCell) = components(:,c);
                    end
                end
            end
            stateInlet = model.initStateFunctionContainers(stateInlet);
            if model.thermal
                stateInlet = computeFlashGeothermal(model.parentModel, stateInlet);
            end
            [mobIn, rhoIn, rhoS] = model.parentModel.getProps(stateInlet, 'Mobility', 'Density', 'SurfaceDensity');
            mobIn = applyFunction(@(x) x(model.G.cells.topCell), mobIn);
            rhoIn = applyFunction(@(x) x(model.G.cells.topCell), rhoIn);
            rhoS  = applyFunction(@(x) x(model.G.cells.topCell), rhoS);

            q = -T.*(p - pi);

            [Qph, qph] = deal(cell(1, nph));
            for ph = 1:nph
                % TODO: Add gravity effects
                flag = q > 0;
                mobRho  = flag.*mobIn{ph}.*rhoIn{ph} + ~flag.*mob{ph}.*rho{ph};
                Qph{ph} = mobRho.*q;
                qph{ph} = Qph{ph}./rhoS{ph};
            end
            % TODO: account for compositional flow
            Q = Qph;

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function Qh = computeInletHeatFlux(model, state, W)
            
            [Q, ~, stateInlet] = computeInletFlux(model, state, W);
%             stateInlet.T(model.G.cells.topCell) = vertcat(W.T);
            
            h   = model.parentModel.getProps(state, 'PhaseEnthalpy');
            h = h{1}(model.G.cells.topCell);
            hIn = model.parentModel.getProps(stateInlet, 'PhaseEnthalpy');
            hIn = hIn{1}(model.G.cells.topCell);

            flag = Q{1} > 0;
            Qh = Q{1}.*(flag.*hIn + ~flag.*h);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            [convergence, values, names] = checkConvergence@WrapperModel(model, problem, varargin{:});
            
            dt = problem.dt;
            Q  = abs(value(problem.state.inletFlux{1}));
%             Q = max(Q);
            Q = Q(model.G.cells.wellNo);
            Qh = abs(value(problem.state.inletHeatFlux));
%             Qh = max(Qh);
            Qh = Qh(model.G.cells.wellNo);
            [mass, energy] = model.parentModel.getProps(problem.state, 'ComponentTotalMass', 'TotalThermalEnergy');
            
            scaleMass = max(value(mass{1})./dt, Q);
            
            scaleEnergy = max(value(energy)./dt, Qh);
            
            values(1) = norm(value(problem.equations{1}./scaleMass), inf);
            values(2) = norm(value(problem.equations{2}./scaleEnergy), inf);
            convergence(1:2) = values(1:2) < model.parentModel.nonlinearTolerance;
%             convergence(1:2) = values(1:2) < 1e-2;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            
            [state, report] = updateAfterConvergence@WrapperModel(model, state0, state, dt, drivingForces);
            
            state = model.initStateFunctionContainers(state);
            Qh = model.computeInletHeatFlux(state, drivingForces.W);
            state.effect = Qh;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function n = numWells(model)
            
            n = numel(model.W);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function cells = getWellCells(model)
            
            cells = model.G.cells.global;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function r = getSegmentRadius(model, face)
            
            % Get radius in each cell
            r = vertcat(model.W.r);
            r = r(model.G.cells.order);
            
            if nargin < 2 || ~face, return; end
            
            % Face segment radii computed as mean of adjacent cells
            N = model.G.faces.neighbors; N = N(all(N>0,2), :);
            r = mean(r(N),2);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function l = getSegmentLength(model, face)
            
            % Compute segments within each cell
            l = arrayfun(@(W) sqrt(sum(W.trajectory.vec.^2, 2)), ...
                        model.W, 'UniformOutput', false);
            l = cell2mat(l);
            l = l(model.G.cells.order);
            
            if nargin < 2 || ~face, return; end
            
            % Face segment lengths computed as mean of adjacent cells
            N = model.G.faces.neighbors; N = N(all(N>0,2), :);
            l = sum(l(N).*0.5, 2);

        end
        %-----------------------------------------------------------------%
        
        function names = getComponentNames(model)
            names = model.parentModel.getComponentNames();
        end
        
        function ncomp = getNumberOfComponents(model)
            ncomp = model.parentModel.getComponentNames();
        end
        
        function nph = getNumberOfPhases(model)
            nph = model.parentModel.getNumberOfPhases();
        end
        
        function active = getActivePhases(model)
            active = model.parentModel.getActivePhases();
        end
        

    end
    
end

function rock = extractSubrock(rock, cells)

    n = numel(rock.poro);
    names = fieldnames(rock);
    for name = names'
        v = rock.(name{1});
        if size(v,1) == n
            rock.(name{1}) = v(cells,:);
        end
    end

end