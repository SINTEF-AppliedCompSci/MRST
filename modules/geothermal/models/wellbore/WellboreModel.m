classdef WellboreModel < WrapperModel
    
    properties
        
        wells
        trajectories
        groups
        allowCtrlSwitching
        shutCtrlType = 'rate';
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
                [wells, G] = processWellboreTrajectories( ...
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
            model.wells = model.processWells(wells);
            model.trajectories = vertcat(wells.trajectory);
            % Process groups
            if ~isempty(opt.groups)
                model.groups = model.processGroups(opt.groups);
            end
            model.allowCtrlSwitching = true(model.numWells() + model.numGroups(),1);

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
            model.operators = model.parentModel.operators;

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function wells = processWells(model, wells)
        % Process wells to see that they are consistently defined
           
            nw = numel(wells);
            for w = 1:nw
                well = wells(w);
                well = model.setMissingLimits(well);
                wells(w) = well;
            end
        
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function groups = processGroups(model, groups)
        % Process groups to see that they are consistently defined
            
            ng = numel(groups);
            wnames = {model.wells.name};
            memberIx = cell(model.numGroups, 1);
            for g = 1:ng
                group = groups(g);
                group = model.setMissingLimits(group);
                groups(g) = group;
                mix = find(ismember(wnames, group.members));
                memberIx{g} = reshape(mix, [], 1);
                
            end
            [groups.memberIx] = deal(memberIx{:});
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function G = computeWellGridGeometry(model, G)
            
            is_fractured = isfield(G.cells, 'hybrid');
            
            % Add inlet cells as copies of the top cell
            topCells = G.cells.topCell;
            [G, wcno] = copyCells(G, G.cells.topCell, 'nnc', repmat(topCells, 1, 2));
            % Set vertical coordinate equal to refDepth
            refDepth = arrayfun(@(W) W.refDepth, model.wells);
            G.cells.centroids(wcno,3) = refDepth;
            G.cells.global = [G.cells.global; G.cells.global(topCells)];
            if is_fractured
                G.cells.hybrid = [G.cells.hybrid; zeros(numel(wcno),1)]; 
            end
            
            gcno = [];
            if ~isempty(model.groups)
                % For groups, we use the top cell of the first member in
                % the group
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                gw  = arrayfun(@(group) group.memberIx, model.groups, 'UniformOutput', false);
                % Make nncs
                nnc = [rldecode(gw1, cellfun(@numel, gw), 1), vertcat(gw{:})];
                nnc = wcno(nnc); if size(nnc,2) == 1; nnc = nnc'; end
                % Copy cells
                [G, gcno] = copyCells(G, wcno(gw1), 'nnc', nnc);
                % Set vertical coord equal to refDepth of first group well
                G.cells.centroids(gcno,3) = refDepth(gw1);
                G.cells.global = [G.cells.global; G.cells.global(gw1)];
                if is_fractured
                    G.cells.hybrid = [G.cells.hybrid; zeros(numel(gcno),1)]; 
                end
            end
            
            % Add heat transmissibility fields for nncs
            [G.nnc.transHr, G.nnc.transHf] = deal(nan(size(G.nnc.trans)));
  
            % Set flag indicating if cell is perforation (i.e., connected
            % to a reservoir cell). Otherwise, it is an inlet(/outlet) cell
            G.cells.type = zeros(G.cells.num, 1);
            G.cells.type(wcno) = 1;
            G.cells.type(gcno) = 2;
            
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
            
            G.connections = struct();
            G.connections.type = zeros(G.faces.num + numel(G.nnc.trans),1);
            G.connections.type((1:model.numWells) + G.faces.num) = 1;
            G.connections.type(G.faces.num+model.numWells+1:end) = 2;
            
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
             % Add bhp if not already set
            if ~isfield(state, 'qt')
                state.qt = zeros(model.numWells() + model.numGroups(), 1);
            end
            % Add surface rates if not already set
            rnames = model.getSurfaceRateNames();
            for rname = rnames
                if ~isfield(state, rname{1})
                    state.(rname{1}) = zeros(model.numWells() + model.numGroups(), 1);
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
                state.effect = zeros(model.numWells() + model.numGroups(), 1);
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@WrapperModel(model);
            forces.groups = [];
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
                case {'qt'}
                    fn = 'qt';
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
            
            % Change flow discretization
            fd = model.parentModel.FlowDiscretization;
            % Replace component phase flux state function
            fd = fd.setStateFunction('ComponentPhaseFlux', WellboreComponentPhaseFlux(model));
            % Replace phase upwind flag state function
            fd = fd.setStateFunction('PhaseUpwindFlag', WellborePhaseUpwindFlag(model));
            % Add surface rate
            fd = fd.setStateFunction('SurfaceRate', WellboreSurfaceRate(model));
            % Add effect
            fd = fd.setStateFunction('Effect', WellboreEffect(model));
            % Set to model
            model.parentModel.FlowDiscretization = fd;
            
            % Change flow property functions
            fp = model.parentModel.FlowPropertyFunctions;
            % Add bottomhole pressure
            fp = fp.setStateFunction('BottomholePressure'   , WellboreBottomholePressure(model));
            % Add bottomhole temperature
            fp = fp.setStateFunction('BottomholeTemperature', WellboreBottomholeTemperature(model));
            % Set Wellbore friction loss
            fp = fp.setStateFunction('FrictionLoss', WellboreFrictionLoss(model));
            
            % Set to model
            model.parentModel.FlowPropertyFunctions = fp;
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            
            % Make wellbore grid and set to model
            nc0 = G.cells.num;
            GW = model.computeWellGridGeometry(G);
            
            % Make rock structure to represent wellbore
            rock = model.parentModel.rock;
            % Equivalent permeability from laminar pipe flow assumption.
            % This is now obsoloete due to dofs for the segment mass fluxes
            % + wellbore friction loss, but still kept here for reference.
            rock.perm = GW.cells.radius.^2/8;
            rock.poro = ones(GW.cells.num, 1);
            fnames = fieldnames(rock)';
            % Add values to rock fields for inlet cells
            cells = GW.cells.topCell;
            if ~isempty(model.groups)
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                cells = [cells; gw1];
            end
            for name = fnames
                v = rock.(name{1});
                if size(v,1) == nc0
                    rock.(name{1}) = [v; v(cells,:)];
                end
            end
            % Set rock and grid to model/parentModel
            model.parentModel.rock         = rock;
            [model.parentModel.G, model.G] = deal(GW);
            
            % Compute operators in standard way
            model.parentModel = model.parentModel.setupOperators(GW, rock, varargin{:});
            op = model.parentModel.operators;
            names = {'T', 'Thr', 'Thf'};
            [~, iif] = model.getInletSegments();
            for name = names
                if ~isfield(op, name{1}), continue; end
                fix = isnan(op.(name{1}));
%                 assert(~any(fix(~iif)), 'Computed non-finite transmissibilites');
                op.(name{1})(fix) = 0;
                op.([name{1}, '_all'])(op.internalConn) = op.(name{1});
            end
            model.parentModel.operators = op;
            
            % Set well indices
            WI = vertcat(model.wells.WI);
            WI = WI(model.G.cells.order);
            model.parentModel.operators.WI = WI;
            
            % Make mapping from interfaces to wells/groups
            nw = model.numWells(); ng = model.numGroups();
            mix = arrayfun(@(group) group.memberIx, model.groups, 'UniformOutput', false);
            ngw = cellfun(@numel, mix);
            ii = [(1:model.numWells)'; rldecode((1:numel(mix))', ngw, 1) + nw];
            jj = [(1:model.numWells)'; (1:sum(ngw))' + nw];
            FS = sparse(ii, jj, 1, nw + ng, nw + sum(ngw));
            
            model.parentModel.operators.fluxSum = @(v) FS*v;
            
            % Make mapping from wells to groups
            ii = [(1:model.numWells)'; rldecode((1:numel(mix))', ngw, 1) + nw];
            jj = [(1:model.numWells)'; vertcat(mix{:})];
            GS = sparse(ii, jj, 1, nw + ng, nw + sum(ngw));
            model.parentModel.operators.groupSum = @(v) GS*v;
            
            if model.thermal
                WIth = vertcat(model.wells.WIth);
                WIth = WIth(model.G.cells.order);
                model.parentModel.operators.WIth = WIth;
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
        % Get primary variables for model. These are the parent model
        % primary variables, in addition to interface mass fluxes
            
            % Get primary variable set of parent model
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
            % Add mass flux
            [vm, qt] = model.getProps(state, 'massFlux', 'qt');
            vars = [vars, {vm, qt}];
            names  = [names, {'massFlux', 'qt'}];
            origin = [origin, class(model), class(model)];

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [model, state] = updateForChangedControls(model, state, forces)
        % Update model and state when controls/drivingForces has changed
        
            model.wells  = forces.W;
            if ~isempty(model.wells)
                model.wells = model.processWells(model.wells);
            end
            
            if ~isfield(forces, 'groups'), return; end
            model.groups = forces.groups;
            if ~isempty(model.groups)
                model.groups = model.processGroups(model.groups);
            end
            
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
            
            state = model.applyLimits(state, drivingForces);
            % Get parent model equations
            [eqs, names, types, state] ...
                = model.parentModel.getModelEquations(state0, state, dt, drivingForces);
            % Mass fluxes
            [fluxEqs, fluxNames, fluxTypes, state] = model.getMassFluxEquations(state);
            % Set controls
            [ctrlEqs, ctrlNames, ctrlTypes, state, eqs] = model.setControls(state, eqs, drivingForces);
            % Assemble
            eqs   = [eqs  , fluxEqs  , ctrlEqs  ];
            names = [names, fluxNames, ctrlNames];
            types = [types, fluxTypes, ctrlTypes];
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getMassFluxEquations(model, state)
        % Equations realting segment mass fluxes to the pressure by means
        % of a wellbore friction model

            
            dp = model.getProp(state, 'FrictionLoss');
            [p, s, rho] = model.parentModel.getProps(state, ...
                'pressure'    , ...
                's'           , ...
                'Density'       ...
            );
            
            nph = model.getNumberOfPhases();
            
            rhoMix = 0;
            if ~iscell(s), s = expandMatrixToCell(s); end
            for ph = 1:nph
                rhoMix = rhoMix + rho{ph}.*s{ph};
            end
            
            rhoMix = model.parentModel.operators.faceAvg(rhoMix);
            
            g   = norm(model.parentModel.gravity);
            dz  = model.parentModel.operators.Grad(model.G.cells.centroids(:,3));
            dpw = model.parentModel.operators.Grad(p);
            
            pot   = dpw - rhoMix.*g.*dz;
            
            eqs   = (pot - dp)./(1*atm);
            eqs   = {eqs};
            names = {'flux'};
            types = {'face'};

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function  [eqs, names, types, state, consEqs] = setControls(model, state, consEqs, drivingForces)
        % Equations imposing well control for each well
        
            % Get control types, targets, and mask for open wells
            [type, target, compi, targetTemp, is_open] = model.getControls(state, drivingForces);
            
            % Convenient well and group indices
            nw = model.numWells();
            ng = model.numGroups();
            [iic, iif] = model.getInletSegments();
        
            % Make mask types
            % No control - well is controlled by group or group is not
            % actively controlling wells
            is_none = strcmpi(type, 'none') & is_open;
            % Well is controlled by its group
            is_group = strcmpi(type, 'group') & is_open;
            % Treat wells controlled by a group as "no control"
            is_none = is_none | is_group;
            % Bottom-hole pressure control
            is_bhp  = strcmpi(type, 'bhp') & is_open;
            % Surface total rates
            is_rate = (strcmp(type, 'rate') | strcmpi(type, 'vrat')) & is_open;
            % Surface oil rates
            is_orat = strcmp(type, 'orat') & is_open;
            % Surface water rates
            is_wrat = strcmp(type, 'wrat') & is_open;
            % Surface gas rates
            is_grat = strcmp(type, 'grat') & is_open;
            % Surface liquid rates (water + oil)
            is_lrat = strcmp(type, 'lrat') & is_open;
            % Reservoir rates (at averaged conditions from previous step)
            is_resv = strcmp(type, 'resv') & is_open;
            % Reservoir rates (at current conditions for each perf.)
            is_volume = strcmp(type, 'volume') & is_open;
            % Check for unsupported well controls
            assert(~any(is_resv | is_volume), 'Not implemented yet');
            
            % Set BHP control equations
            nph = model.parentModel.getNumberOfPhases();
            bhp = model.getProps(state, 'BottomholePressure');
            
            ctrlEqs = model.AutoDiffBackend.convertToAD(zeros(nw + ng, 1), bhp);
            ctrlEqs(is_bhp) = bhp(is_bhp) - target(is_bhp);

            % Set rate control equations
            phases = model.parentModel.getPhaseNames();
            is_surface_rate = false(nw + ng, 1);
            % Find rate targets
            for ph = 1:nph
                switch phases(ph)
                    case 'W'
                        act = is_rate | is_wrat | is_lrat;
                    case 'O'
                        act = is_rate | is_orat | is_lrat;
                    case 'G'
                        act = is_rate | is_grat;
                end
                is_surface_rate(act) = true;
            end
            %@@ TODO: support multiphase flow
            % Convert surface volume rates to mass rates
            [rhoS, Qt] = model.getProps(state, 'SurfaceDensity', 'qt');
            qs = Qt./rhoS{1}(iic);
            % Add in as source in corresponding well cells
            ctrlEqs(is_surface_rate) = target(is_surface_rate) - qs(is_surface_rate);
            
            % Set zero rate for wells that do not have a control
%             ctrlEqs(is_none | ~is_open) = Qt(is_none | ~is_open);
            is_shut = is_none | ~is_open;
            if any(is_shut)
                switch model.shutCtrlType
                    case 'rate'
                        ctrlEqs(is_shut) = Qt(is_shut);
                    case 'bhp'
                        cells = model.operators.N(iif,2);
                        dz = model.G.cells.centroids(cells,3);
                        g = norm(model.parentModel.gravity);
                        rho = model.getProp(state, 'Density');
                        rho = rho{1}(cells);
                        pShut = 1*atm + dz.*g.*rho;
                        ctrlEqs(is_shut) = bhp(is_shut) - pShut;
                end
            end
            
            % Format output
            eqs   = {ctrlEqs  };
            names = {'control'};
            types = {'control'};
            
            % Add source term to mass conservation equation
            consEqs{1}(iic) = consEqs{1}(iic) - Qt;
            
            if ~model.thermal, return; end
            
            % Compute injection/production enthalpy
            [h, bht] = model.getProps(state, 'enthalpy', 'BottomHoleTemperature');
            is_inj = value(Qt) > 0;
            fix = isnan(value(targetTemp));
            if ~isa(targetTemp, 'ADI'), bht = value(bht); end
            targetTemp(fix) = bht(fix);
            hInj = model.parentModel.fluid.hW(bhp, targetTemp);
            h(iic) = is_inj.*hInj + ~is_inj.*h(iic);

            % Add source term to energy conservation equation
            consEqs{2}(iic) = consEqs{2}(iic) - Qt.*h(iic);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [type, target, compi, targetTemp, is_open] = getControls(model, state, drivingForces)
            
            % Get driving forces and check that it's non-empty
            [wCtrl, gCtrl] = deal([]);
            hasWellCtrl  = isfield(drivingForces, 'W') && ~isempty(drivingForces.W);
            hasGroupCtrl = isfield(drivingForces, 'groups') && ~isempty(drivingForces.groups);
            if hasWellCtrl , wCtrl = state.wells;  end
            if hasGroupCtrl, gCtrl = state.groups; end
            if isempty(wCtrl) && isempty(gCtrl), return; end
        
            nw = model.numWells();
            ng = model.numGroups();
            
            nph = model.getNumberOfPhases;
            compi = ones(1, nph)./nph;
            ctrl0 = struct('type', 'none', 'val' , nan, ...
                'compi', compi, 'T', nan, 'status', true);
            if ~hasWellCtrl , wCtrl = repmat(ctrl0, nw, 1); end
            if ~hasGroupCtrl, gCtrl = repmat(ctrl0, ng, 1); end
        
            target = [vertcat(wCtrl.val); vertcat(gCtrl.val)];
            type   = vertcat({wCtrl.type}', {gCtrl.type}');
            compi  = [vertcat(wCtrl.compi); vertcat(gCtrl.compi)];
            
            targetTemp = [];
            if model.thermal
                targetTemp = [vertcat(wCtrl.T); vertcat(gCtrl.T)];
            end
            
            % Wells that are open to flow
            is_open = vertcat(wCtrl.status, gCtrl.status);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function state = applyLimits(model, state, drivingForces)
        % Update solution variables and wellSol based on the well limits.
        % If limits have been reached, this function will attempt to
        % re-initialize the values and change the controls so that the next
        % step keeps within the prescribed ranges.
        
        
            if ~isfield(state, 'wells'), state.wells = model.wells; end
            if ~isfield(state, 'groups'), state.groups = model.groups; end
            % Get driving forces and check that it's non-empty
            [wCtrl, gCtrl] = deal([]);
            hasWellCtrl  = isfield(drivingForces, 'W') && ~isempty(drivingForces.W);
            hasGroupCtrl = isfield(drivingForces, 'groups') && ~isempty(drivingForces.groups);
            if hasWellCtrl , wCtrl = state.wells;  end
            if hasGroupCtrl, gCtrl = state.groups; end
            if isempty(wCtrl) && isempty(gCtrl), return; end
        
            nw = model.numWells();
            ng = model.numGroups();
            
            nph = model.getNumberOfPhases;
            compi = ones(1, nph)./nph;
            ctrl0 = struct('type', 'none', 'val' , nan, ...
                'compi', compi, 'T', nan, 'status', true);
            if ~hasWellCtrl , wCtrl = repmat(ctrl0, nw, 1) ; end
            if ~hasGroupCtrl, gCtrl = repmat(ctrl0, ng, 1); end
        
            target = [vertcat(wCtrl.val); vertcat(gCtrl.val)];
            type   = vertcat({wCtrl.type}', {gCtrl.type}');
            compi  = [vertcat(wCtrl.compi); vertcat(gCtrl.compi)];
            
            targetTemp = [];
            if model.thermal
                targetTemp = [vertcat(wCtrl.T); vertcat(gCtrl.T)];
            end
            
            % Wells that are open to flow
            is_open = vertcat(wCtrl.status, gCtrl.status);
            
            limits = model.getLimits();
            
            % Convenient well and group indices
            nw = model.numWells();
            ng = model.numGroups();
            [iic, iif] = model.getInletSegments();
            
            [bhp, rhoS, Qt] = model.getProps(state, ...
                'BottomholePressure', 'SurfaceDensity', 'qt');
            qs = Qt./rhoS{1}(iic);
            is_inj = vertcat(model.wells.sign) > 0;
            
            lims = [limits.bhp, limits.rate];
            vals = [value(bhp), value(qs)  ];
            d    = vals - lims;
            
            flags = false(nw + ng, 2);
            flags( is_inj,:) = d( is_inj,:) > 0;
            flags(~is_inj,:) = d(~is_inj,:) < 0;
                
            % limits we need to check (all others than w.type):
            
            modes = fieldnames(limits)';
            check = cellfun(@(t) ~strcmpi(t,modes), type, 'UniformOutput', false);
            check = cell2mat(check);
            
            violate       = flags & check;
            has_violation = any(violate, 2);
            
            modes = repmat(modes, nw + ng, 1);
            type(has_violation)   = modes(violate);
            target(has_violation) = lims(violate);
            
            for i = 1:nw
                if has_violation(i)
                    fprintf('Well %s: Control mode changed from %s to %s.\n', ...
                        model.wells(i).name, state.wells(i).type, type{i});
                end
                state.wells(i).type = type{i};
                state.wells(i).val  = target(i);
                state.wells(i).T    = targetTemp(i);
            end
            
            for i = nw + (1:ng)
                state.groups(i-nw).type = type{i};
                state.groups(i-nw).val  = target(i);
                state.groups(i-nw).T    = targetTemp(i);
            end
            
        end

        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function facility = setMissingLimits(model, facility)
            
            limits = facility.lims;
            if isempty(limits), limits = struct(); end
            type = setdiff(facility.type, {'group', 'none'});
            if ~isempty(type) && ~isfield(limits, type)
                limits.(type{:}) = facility.val;
            end
            
            modes = {'bhp', 'rate'};
            modes = setdiff(modes, facility.type);
            missing = modes(~isfield(limits, modes));
            
            sgn = sign(facility.sign);
            if sgn == 0, sgn = 1; end
            val = sgn.*inf;
            for f = missing
                v = val; %if strcmpi(f, 'vrat'), v = -v; end
                limits = setfield(limits, f{:}, v);
            end
            facility.lims = limits;

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, forces)

            parent = strcmpi(problem.types, 'cell');

            % Update parent model variables
            p = problem; p.primaryVariables = p.primaryVariables(parent);
            [state, report] = updateState@WrapperModel(model, state, p, dx, forces);

            
            for i = (nnz(parent) + 1):numel(problem.primaryVariables)
                dxAbsMax = inf;
                p = problem.primaryVariables{i};
%                 if strcmpi(p, 'massFlux'), dxAbsMax = 1e-3*mean(model.G.faces.areas)*litre/second; end
                % Update wellbore model variables
                state = model.updateStateFromIncrement(state, dx, problem, p, inf, dxAbsMax);
            end
            
%             sgn = sign(state.massFlux);
%             v = max(abs(state.massFlux), 1e-3);
%             state.massFlux = v.*sgn;
%             model = model.applyLimits(state);

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            
            [state, report] = model.parentModel.updateAfterConvergence(state0, state, dt, drivingForces);
            [state.bhp, qr] = model.getProps(state, ...
                'BottomholePressure', 'massFlux');

            [s, rhoS, flag] = model.parentModel.getProps(state, ...
                's', 'SurfaceDensity', 'PhaseUpwindFlag');
            
            
            [iic, iif] = model.getInletSegments();
            % Upstream-wheight density
            rhoSf = cellfun(@(flag, rho) ...
                model.parentModel.operators.faceUpstr(flag, rho), ...
                flag, rhoS, 'UniformOutput', false);
            % Set surface rate
            state.qWs = model.parentModel.operators.groupSum(qr(iif)./rhoSf{1}(iif));
            
            if model.thermal
                % Set effect
                state.bht = model.getProp(state, 'BottomholeTemperature');
                qh        = model.parentModel.getProp(state, 'HeatFlux');
                state.effect = model.parentModel.operators.groupSum(qh(iif));
            end
            
            state = rmfield(state, {'wells', 'groups'});

            if ~model.thermal, return; end
            
            [wCtrl, gCtrl] = deal([]);
            hasWellCtrl  = isfield(drivingForces, 'W') && ~isempty(drivingForces.W);
            hasGroupCtrl = isfield(drivingForces, 'groups') && ~isempty(drivingForces.groups);
            if hasWellCtrl , wCtrl = drivingForces.W;      end
            if hasGroupCtrl, gCtrl = drivingForces.groups; end
            
            nw = model.numWells();
            ng = model.numGroups();
            
            ctrl0 = struct('type', 'none', 'val', nan, 'T', nan);
            if ~hasWellCtrl , wCtrl = repmat(ctrl0, nw, 1); end
            if ~hasGroupCtrl, gCtrl = repmat(ctrl0, ng, 1); end
        
            ctrlType = vertcat({wCtrl.type}', {gCtrl.type}');
            is_none  = strcmpi(ctrlType, 'none');
            
            CpW  = model.parentModel.fluid.CpW(state.bhp, state.bht);
            rhoW = model.parentModel.fluid.rhoW(state.bhp, state.bht);
            h = state.effect./(state.qWs.*rhoS{1}(iic));
            T = (h - state.bhp./rhoW)./CpW;
            % Set bottom-hole temperature in inactive groups
            state.bht(is_none) = T(is_none);
            
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
        function cells = getGlobalWellCells(model)
        % Get cell numbering of perforated cells in the reservoir grid
        
            cells = model.G.cells.global(model.G.cells.type == 0);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [iiCell, iiFace] = getInletSegments(model, G)
            
            if nargin < 2, G = model.G; end
            iiCell = G.cells.type ~= 0;
            iiFace = any(iiCell(model.parentModel.operators.N),2);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [rc, rf] = getSegmentRadius(model, G)
            
            if nargin < 2, G = model.G; end
            
            % Get radius in each cell
            rc = vertcat(model.wells.r);
            rc = rc(G.cells.order);
            
            % Add radii for inlet well cells
            rc = [rc; rc(G.cells.topCell)];
            
            if ~isempty(model.groups)
            % Add radii for inlet group cells
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                rc = [rc; rc(G.cells.topCell(gw1))];
            end
            
            % Face segment radii computed as mean of adjacent cells
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            N = [N; G.nnc.cells];
            rf = mean(rc(N),2);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [lc, lf] = getSegmentLength(model, G)
            
            if nargin < 2, G = model.G; end
            
            % Compute lenght of segments within each cell
            lc = arrayfun(@(W) sqrt(sum(W.trajectory.vec.^2, 2)), ...
                        model.wells, 'UniformOutput', false);
            lc = cell2mat(lc);
            lc = lc(G.cells.order);

            % Compute lenght of segments for each inlet cell
            cells = G.cells.topCell;
            if model.numGroups() > 0
                % Get cells that group cells were copied from
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                cells = [cells; gw1];
            end
            lc = [lc; lc(cells)];
            
            % Face segment lengths computed as mean of adjacent cells            
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            N = [N; G.nnc.cells];
            lf = sum(lc(N).*0.5, 2);
%             if model.numGroups > 0
%                 x  = G.cells.centroids;
%                 li = sqrt(sum((x(G.nnc.cells(:,1), :) ...
%                                 - x(G.nnc.cells(:,2), :)).^2, 2)); 
%                 lf = [lf; li];
%             end

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [rc, rf] = getSegmentRoughness(model, G)
            
            if nargin < 2, G = model.G; end
             
            % Get wellbore rougness (currently assumes one value per well)
            r0 = arrayfun(@(trajectory) sum(trajectory.roughness, 2), ...
                    model.trajectories, 'UniformOutput', false);
            rc = cell2mat(r0);
            assert(numel(rc) == model.numWells, ['WellboreModel ', ...
                'currently only supports one roughness value per well']);
            rc = rc(G.cells.wellNo);
            cells = G.cells.topCell;
            if model.numGroups() > 0
                % Get cells that group cells were copied from
                gw1 = arrayfun(@(group) group.memberIx(1), model.groups);
                cells = [cells; gw1];
            end
            rc = [rc; rc(cells)];
            
            % Face segment roughness computed as mean of adjacent cells
            N = G.faces.neighbors; N = N(all(N>0,2), :);
            N = [N; G.nnc.cells];
            rf = sum(rc(N).*0.5, 2);

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function active = getActiveWells(model)
            
            active = vertcat(model.wells.status);
            
        end            
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function limits = getLimits(model)
            
            limits = vertcat(model.wells.lims);
            if ~isempty(model.groups)
                limits = [limits; vertcat(model.groups.lims)];
            end
            names  = fieldnames(limits);
            values = num2cell(cell2mat(struct2cell(limits))', 1);
            limits = cell2struct(values, names, 2);
            
        end            
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function types = getTypes(model)
            types = {model.wells.types};
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        % Conveient parent model calls
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
        % Get names of components in parent model
            
            names = model.parentModel.getComponentNames();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function ncomp = getNumberOfComponents(model)
        % Get number of components in parent model
            
            ncomp = model.parentModel.getNumberOfComponents();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function names = getPhaseNames(model)
        % Get names of phases in parent model
            
            names = model.parentModel.getPhaseNames();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function nph = getNumberOfPhases(model)
        % Get number of phases in parent model
            
            nph = model.parentModel.getNumberOfPhases();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function active = getActivePhases(model)
        % Get active phases in parent model
            
            active = model.parentModel.getActivePhases();
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function rhoS = getSurfaceDensities(model, varargin)
            
            rhoS = model.parentModel.getSurfaceDensities(varargin{:});
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        % Visualization
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function plotWellNetwork(model)
            
            nnc = model.G.nnc.cells(model.numWells + 1:end, :);
            nnc = max(nnc(:)) - nnc + 1;
            names = [{model.wells.name}, {model.groups.name}];
            names = fliplr(names);
            network = graph(nnc(:,1), nnc(:,2), [], names);
            
            plot(network, 'layout', 'layered', 'LineWidth', 2, 'MarkerSize', 8);
            
        end
        %-----------------------------------------------------------------%

    end
    
end

%-------------------------------------------------------------------------%
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
%-------------------------------------------------------------------------%

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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