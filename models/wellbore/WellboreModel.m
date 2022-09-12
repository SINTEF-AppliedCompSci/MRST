classdef WellboreModel < WrapperModel
    
    properties
        
        wells
        groups
        water, oil, gas
        thermal = false
        
    end
   
    methods
        
        %-----------------------------------------------------------------%
        function model = WellboreModel(rmodel, trajectories, varargin)
            
            opt = struct('names' , {{}}, ...
                         'WI'    , []  , ...
                         'groups', []);
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
            model.wells = W;
            % Process groups
            if ~isempty(opt.groups)
                model.groups = model.processGroups(opt.groups);
            end

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
        function groups = processGroups(model, groups)
            
            ng = numel(groups);
            wnames = {model.wells.name};
            memberIx = cell(model.numWells, 1);
            for i = 1:ng
                group = groups(i);
                memberIx{i} = find(ismember(wnames, group.members));
            end
            [groups.memberIx] = deal(memberIx{:});
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function G = computeWellGridGeometry(model, G)
            
            % Add inlet cells as copies of the top cell
            G.cells.num = G.cells.num + model.numWells + model.numGroups;
            % Get top cells for each well
            topCells = G.cells.topCell;
            if ~isempty(model.groups)
                % For groups, this we use the top cell of the first member
                % in the group
                gwix = arrayfun(@(group) group.memberIx(1), model.groups);
                topCells = [topCells; topCells(gwix)];
            end
            % Set volumes
            G.cells.volumes = [G.cells.volumes  ; G.cells.volumes(topCells)];
            % Set centroids, with depth equal to refDepth
            refDepth = arrayfun(@(W) W.refDepth, model.wells);
            if ~isempty(model.groups)
               refDepth = [refDepth; refDepth(gwix)]; 
            end
            x = G.cells.centroids(topCells,:); x(:,3) = refDepth;
            G.cells.centroids = [G.cells.centroids; x];
            % Set facePos
            % @@TODO: Add one top face for each well
            nf = diff(G.cells.facePos); nf = nf(topCells);
            facePos = cumsum([G.cells.facePos(end); nf]);
            G.cells.facePos = [G.cells.facePos; facePos(2:end)];
            % Set faces
            faces = G.cells.faces(mcolon(G.cells.facePos(topCells), ...
                                         G.cells.facePos(topCells + 1) - 1),:);
            G.cells.faces = [G.cells.faces; faces];
            % Set inlet cells to be neighbors of the topcells
            topFaces = G.faces.topFace;
            if ~isempty(model.groups)
                topFaces = [topFaces; G.face.num + model.numGroups];
            end
            N = G.faces.neighbors(topFaces,:);
            N = [N; zeros(model.numGroups)];
            N(N == 0) = (G.cells.num-model.numWells-model.numGroups+1:G.cells.num);
            
            G.faces.neighbors(G.faces.topFace,:) = N;
             
            % Set flag indicating if cell is perforation (i.e., connected
            % to a reservoir cell). Otherwise, it is an inlet(/outlet) cell
            G.cells.isPerf = true(G.cells.num,1);
            G.cells.isPerf(end-model.numWells+1:end) = false;
            
            % Compute segment lengths and radii
            [lengthCell, lengthFace] = model.getSegmentLength(G);
            [radiusCell, radiusFace] = model.getSegmentRadius(G);
            
            % Get cell segment lengths and radii and set to G.cells for
            % easy access later
            G.cells.length = lengthCell;
            G.cells.radius = radiusCell;
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
            G.faces.length = lengthFace;
            G.faces.radius = radiusFace;
            
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
            
            iic = model.getInletSegments();
            
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
                state.bhp = bhp(iic);
            end
            % Add surface rates if not already set
            rnames = model.getSurfaceRateNames();
            for rname = rnames
                if ~isfield(state, rname{1})
                    state.(rname{1}) = zeros(model.numWells(), 1);
                end
            end
            
            % Check if we have a thermal model
            if ~(isprop(model.parentModel, 'thermal') && model.parentModel.thermal)
                return
            end
            
            % Add bht if not already set
            if ~isfield(state, 'bht')
                bht       = model.getProp(state, 'temperature');
                state.bht = bht(iic);
            end
            % Add effect if not already set
            if ~isfield(state, 'effect')
                state.effect = zeros(model.numWells(), 1);
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
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            
            % Make wellbore grid and set to model
            GW = model.computeWellGridGeometry(G);
            
            % Make rock structure to represent wellbore
            rock = model.parentModel.rock;
            % Euivalent permeability from laminar pipe flow assumption.
            % This is now obsoloete due to dofs for the segment mass fluxes
            % + wellbore friction loss, but still kept here for reference.
            rock.perm = min(GW.cells.radius.^2/8, 1e-5);
            rock.poro = ones(GW.cells.num, 1);
            fnames = fieldnames(rock)';
            for name = fnames
                v = rock.(name{1});
                if size(v,1) == GW.cells.num - model.numWells
                    rock.(name{1}) = [v; v(GW.cells.topCell,:)];
                end
            end
            
            model.parentModel.rock = rock;
            model.parentModel.G = GW;
            model.G = GW;
                        
            % Compute operators in standard way
            model.parentModel = model.parentModel.setupOperators(GW, rock, varargin{:});
            op = model.parentModel.operators;
            names = {'T', 'Thr', 'Thf'};
            [~, iif] = model.getInletSegments();
            for name = names
                if ~isfield(op, name{1}), continue; end
                fix = isnan(op.(name{1})(iif));
                op.(name{1})(fix) = 0;
                op.([name{1}, '_all'])(op.internalConn) = op.(name{1});
            end
            model.parentModel.operators = op;
            % Set well indices
            WI = vertcat(model.wells.WI);
            WI = WI(model.G.cells.order);
            model.parentModel.operators.WI = WI;
            
            if model.thermal
                WIth = vertcat(model.wells.WIth);
                WIth = WIth(model.G.cells.order);
                model.parentModel.operators.WIth = WIth;
            end
            
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
                qh  = model.getProp(state, 'effect');
                wvars   = [wvars, {bht, qh}];
                wnames  = [wnames, {'bht', 'effect'}];
                worigin = [worigin, repmat({class(model)}, 1, 2)];
            end
            vars   = [vars  , wvars  ];
            names  = [names , wnames ];
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
            
            % Set well controls (bhp, bht, composition) to inlet node
%             state = model.setControls(state, drivingForces);
            
            % Get parent model equations
            [eqs, names, types, state] ...
                = model.parentModel.getModelEquations(state0, state, dt, drivingForces);
            % Mass fluxes, surface rates and control
            [fluxEqs, fluxNames, fluxTypes, state] = model.getMassFluxEquations(state);
            [rateEqs, rateNames, rateTypes, state] = model.getSurfaceRateEquations(state);
            [effEqs , effNames , effTypes , state] = model.getHeatFluxEquations(state);
            [ctrlEqs, ctrlNames, ctrlTypes, state] = model.getControlEquations(state, drivingForces);
            % Replace conservation equations in inlet nodes with equations
            % ensuring inlet pressure/temperature equal to well bhp/bht
            [p, bhp] = model.getProps(state, 'pressure', 'bhp');
            [T, bht] = model.getProps(state, 'T', 'bht');
            iic = model.getInletSegments();
            eqs{1}(iic) = p(iic) - bhp;
            eqs{2}(iic) = T(iic) - bht;
            % Assemble
            eqs   = [eqs  , fluxEqs  , rateEqs  , effEqs  , ctrlEqs  ];
            names = [names, fluxNames, rateNames, effNames, ctrlNames];
            types = [types, fluxTypes, rateTypes, effTypes, ctrlTypes];
            % Set inlet/outlet mass/heat fluxes to state
            [Q, Qh]  = model.getProps(state, 'ComponentTotalFlux', 'HeatFlux');
            [~, iif] = model.getInletSegments();
            Q = applyFunction(@(x) x(iif), Q);
            state.inletFlux = Q;
            state.inletHeatFlux = Qh;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getMassFluxEquations(model, state)
        % Equations realting segment mass fluxes to the pressure by means
        % of a wellbore friction model
            
            [~, roughness] = model.getSegmentRoughness();

            d = model.G.faces.radius.*2;
            l = model.G.faces.length;

            [v, rho, mu, pot, flag] = model.parentModel.getProps(state, ...
                                        'ComponentPhaseFlux'      , ...
                                        'Density'                 , ...
                                        'Viscosity'               , ...
                                        'PhasePotentialDifference', ...
                                        'PhaseUpwindFlag'         );
            rho{1} = model.parentModel.operators.faceAvg(rho{1});

            mu{1} = model.parentModel.operators.faceUpstr(flag{1}, mu{1});
            dp = wellBoreFriction(v{1}, rho{1}, mu{1}, d, l, roughness, 'massRate');

            bhp = model.getProp(state, 'bhp');
            eqs   = (pot{1} - dp)./max(max(value(bhp)), 1);
            eqs   = {eqs};
            names = {'flux'};
            types = {'face'};

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getSurfaceRateEquations(model, state)
        % Equations realting rates into/out of the wellbore inlet/outlet to
        % those computed based on wellbore friction/bhp
            
            % TODO: Implement support for multiphase/multicomponent
        
            % Get properties
            [qs, rhoS, flag] = model.parentModel.getProps(state, ...
                                'ComponentPhaseFlux', ...
                                'SurfaceDensity', ...
                                'PhaseUpwindFlag');                
            % Upstream-wheight density
            rhoS = cellfun(@(flag, rho) ...
                model.parentModel.operators.faceUpstr(flag, rho), ...
                flag, rhoS, 'UniformOutput', false);
            % Get surface rates across inlet face segment
            [~, iif] = model.getInletSegments();
            qs = cellfun(@(qs, rho) qs(iif)./rho(iif), ...
                    qs, rhoS, 'UniformOutput', false);
            
            % Get surface rate dofs
            nph      = model.getNumberOfPhases();
            names    = model.getSurfaceRateNames();
            qsw      = cell(1, nph);
            [qsw{:}] = model.getProps(state, names{:});
            
            % Assemble surface rate equations
            scale = (meter^3/day)^(-1);
            eqs   = cellfun(@(qs, qph) (qs - qph).*scale, qsw, qs, ...
                            'UniformOutput', false);
            types = repmat({'rate'}, 1, nph);

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getHeatFluxEquations(model, state)
        % Equations realting heat flux into/out of the wellbore
        % inlet/outlet (i.e., effect) to 

            % Get properties
            [qh, qhw, h] = model.getProps(state, ...
                                'HeatFlux', ...
                                'effect', ...
                                'enthalpy');
            [~, iif] = model.getInletSegments();
            qh = qh(iif);
            
            scale = (mean(value(h))*meter^3/day).^(-1); %@@ Set reasonable scaling
            eqs   = {(qh - qhw).*scale};
            names = {'effect'};
            types = {'effect'};
        
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getControlEquations(model, state, drivingForces)
        % Equations imposing well control for each well

            % Get driving forces and check that it's non-empty
            wCtrl = drivingForces.W;
            if isempty(wCtrl)
                [eqs, names, types] = deal({});
                return;
            end
        
            bhp = model.getProp(state, 'bhp');
            eqs = model.AutoDiffBackend.convertToAD(zeros(model.numWells,1), bhp);
            ctrlTarget = vertcat(wCtrl.val);
            ctrlType   = {wCtrl.type}';
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
            eqs(is_bhp) = ctrlTarget(is_bhp) - bhp(is_bhp);

            % Set rate control equations
            nph    = model.parentModel.getNumberOfPhases();
            phases = model.parentModel.getPhaseNames();
            rnames = model.getSurfaceRateNames();
            qsw      = cell(1, nph);
            [qsw{:}] = model.getProps(state, rnames{:});
            is_surface_rate = false(model.numWells(), 1);
            wrates = bhp*0;
            % Sum up rate targets
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
            eqs(is_surface_rate) = ctrlTarget(is_surface_rate) - wrates(is_surface_rate);
            
            % Pack in standard eqs/names/types format
            eqs = {eqs}; [names, types] = deal({'control'});
            
            % Add bht equation if applicable
            if model.thermal
                
                % Get logical masks for cell and face inlet segments
                [iic, iif] = model.getInletSegments();
                
                [bht, h, qhw, rhoS, flag] = model.getProps(state, ...
                                    'bht', 'PhaseEnthalpy', 'effect', ...
                                    'SurfaceDensity', 'PhaseUpwindFlag');
                h = applyFunction(@(x) x(iic), h);
                rhoS = applyFunction(@(x) x(iic), rhoS);
                is_inj = flag{1}(iif);
                qh = 0;
                for ph = 1:nph
                    qh = qh + qsw{ph}.*rhoS{ph}.*h{ph};
                end
                eqTh = qh - qhw;
                Tw = vertcat(wCtrl.T);
                eqTh(is_inj) = bht(is_inj) - Tw(is_inj);
                
                eqs   = [eqs, {eqTh}];
                names = [names, 'thermal'];
                types = [types, 'control'];
            end
            
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
            
            pad = zeros(model.numWells,1);
            scaleMass = max(value(mass{1})./dt, [Q; pad]);
            
            scaleEnergy = max(value(energy)./dt, [Qh; pad]);
            
            values(1) = norm(value(problem.equations{1}./scaleMass), inf);
            values(2) = norm(value(problem.equations{2}./scaleEnergy), inf);
            convergence(1:2) = values(1:2) < model.parentModel.nonlinearTolerance;
%             convergence(1:2) = values(1:2) < 1e-2;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function n = numWells(model)
            
            n = numel(model.wells);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function n = numGroups(model)
            
            n = numel(model.groups);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function cells = getWellCells(model)
            
            cells = model.G.cells.global;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [iiCell, iiFace] = getInletSegments(model, G)
            
            if nargin < 2, G = model.G; end
            iiCell = ~G.cells.isPerf;
            iiFace = any(iiCell(model.parentModel.operators.N),2);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [rc, rf] = getSegmentRadius(model, G)
            
            if nargin < 2, G = model.G; end
            
            % Get radius in each cell
            rc = vertcat(model.wells.r);
            rc = rc(G.cells.order);
            rc = [rc; arrayfun(@(W) W.r(1), model.wells)];
            
            % Face segment radii computed as mean of adjacent cells
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            rf = mean(rc(N),2);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [lc, lf] = getSegmentLength(model, G)
            
            if nargin < 2, G = model.G; end
            
            % Compute segments within each cell
            lc = arrayfun(@(W) sqrt(sum(W.trajectory.vec.^2, 2)), ...
                        model.wells, 'UniformOutput', false);
            lc = cell2mat(lc);
            lc = lc(G.cells.order);
            x = G.cells.centroids;
            iic = model.getInletSegments(G);
            l0 = sqrt(sum((x(G.cells.topCell,:) - x(iic,:)).^2, 2));
            lc = [lc; l0];
            
            % Face segment lengths computed as mean of adjacent cells
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            lf = sum(lc(N).*0.5, 2);

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [rc, rf] = getSegmentRoughness(model, G)
            
             if nargin < 2, G = model.G; end
             
            % Get wellbore rougness (currently assumes one value per well)
            r0 = arrayfun(@(W) sum(W.trajectory.roughness, 2), ...
                    model.wells, 'UniformOutput', false);
            rc = cell2mat(r0);
            assert(numel(rc) == model.numWells, ['WellboreModel ', ...
                'currently only supports one roughness value per well']);
            rc = [rc(G.cells.wellNo); cellfun(@(r) r(1), r0)];
            
            % Face segment roughness computed as mean of adjacent cells
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            rf = sum(rc(N).*0.5, 2);

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
            
            names = model.parentModel.getComponentNames();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function ncomp = getNumberOfComponents(model)
            
            ncomp = model.parentModel.getNumberOfComponents();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function nph = getNumberOfPhases(model)
            
            nph = model.parentModel.getNumberOfPhases();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function active = getActivePhases(model)
            
            active = model.parentModel.getActivePhases();
            
        end
        %-----------------------------------------------------------------%
        

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