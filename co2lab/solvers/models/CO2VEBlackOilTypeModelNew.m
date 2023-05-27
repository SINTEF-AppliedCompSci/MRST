classdef CO2VEBlackOilTypeModelNew < ReservoirModel & GenericReservoirModel
    
    properties
    end

    methods
        function model = CO2VEBlackOilTypeModelNew(Gt, rock2D, fluid, varargin)
           
            model = model@ReservoirModel(Gt, varargin{:});
            model.rock = rock2D;
            model.fluid = fluid;
            [model.water, model.oil, model.gas] = deal(true, false, true);
            model.gravity = norm(gravity);
            
            model = model.setupOperators(Gt, rock2D);
            model.Components = {ImmiscibleComponent('brine', 1), ...
                                ImmiscibleComponent('CO2', 2)};
        end
        
% ------------------------------------------------------------------------        

        function model = setupOperators(model, Gt, rock)
        % Compute vertially-integrated transmissibilities if not provided
            rock_tmp      = rock; 
            rock_tmp.perm = rock.perm .* Gt.cells.H; 
            T             = computeTrans(Gt, rock_tmp); 
            cf            = Gt.cells.faces(:, 1); 
            nf            = Gt.faces.num; 
            T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);       
            
            % Computing vertically-integrated pore - volume
            pv = poreVolume(Gt, rock); 
            pv = pv .* Gt.cells.H; 
            
            model.operators = setupOperatorsTPFA(Gt, rock, 'porv', pv, 'trans', T);
            model.operators.T_all = T;
        end
        
% ------------------------------------------------------------------------

        function [eqs, names ,types, state] = ...
                getModelEquations(model, state0, state, dt, drivingForces)
            
            [acc, flux, names, types] = ...
                model.FlowDiscretization.componentConservationEquations(...
                    model, state, state0, dt);
            
            % add sources
            src = model.FacilityModel.getComponentSources(state);
            acc = model.insertSources(acc, src);
            
            % assemble equations
            eqs = cell(1, numel(acc));
            for i = 1:numel(acc)
                eqs{i} = model.operators.AccDiv(acc{i}, flux{i});
            end
            
            % get well equations
            if ~isempty(model.FacilityModel)
                [weqs, wnames, wtypes, state] = ...
                    model.FacilityModel.getModelEquations(state0, state, dt, ...
                                                                  drivingForces);
                % concatenate
                eqs = [eqs, weqs];
                names = [names, wnames];
                types = [types, wtypes];
            end
            
            % add in boundary conditions
            [eqs, state, src] = ...
                model.addBoundaryConditionsAndSources(eqs, names, types, ...
                                                      state, drivingForces);
        end
        
% ------------------------------------------------------------------------        

        function model = validateModel(model, varargin)
            model = validateModel@ReservoirModel(model, varargin{:});
        end

% ------------------------------------------------------------------------        

        function state = validateState(model, state)
            state = validateState@ReservoirModel(model, state);
        end

% ------------------------------------------------------------------------        

        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@ReservoirModel(model, varargin{:});
        
            % FlowPropertyFunctions
            flowprops = model.FlowPropertyFunctions;
            flowprops = flowprops.setStateFunction('CapillaryPressure', ...
                                                   CO2VECapillaryPressure(model));
            flowprops = flowprops.setStateFunction('RelativePermeability', ...
                                                   CO2VERelativePermeability(model));
            model.FlowPropertyFunctions = flowprops;
            
            % PVTPropertyFunctions
            pvtprops = model.PVTPropertyFunctions;
            pvtprops = pvtprops.setStateFunction(...
                'ShrinkageFactors', ...
                ShrinkageFactors(model, 'usePhasePressures', false));
            
            model.PVTPropertyFunctions = pvtprops;
        end

% ------------------------------------------------------------------------        
        
        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                         drivingForces);
            
            sg          = state.s(:,2);
            sg          = min(1, max(0, sg)); %(1-model.fluid.res_water, max(0,sg)); @@
            state.s     = [1-sg, sg];    

            % ensure sGmax between sg and 1
            state.sGmax = max(min(1, state.sGmax), sg);
            
        end

% ------------------------------------------------------------------------        
        
        function [state, report] = ...
                updateAfterConvergence(model, state0, state, dt, drivingForces)
      
            [state, report] = ...
                updateAfterConvergence@ReservoirModel(model, ...
                                                      state0, state, dt, drivingForces);
            
            if model.outputFluxes
                % add internal fluxes
                qWG = model.getProp(state, 'PhaseFlux');
                state = storeFluxes(model, state, qWG{1}, [], qWG{2});
                
                % No need to add boundary fluxes, as this is already done 
                % in the call to 'addBoundaryConditionsAndSources' when 
                % setting up the equations
            end
            
            sGmax0 = model.getProp(state0, 'sGmax');
            sG     = model.getProp(state, 'sg');
            state = model.setProp(state, 'sGmax', max(sG, sGmax0));
            
        end

% ------------------------------------------------------------------------        

        function [fn, index] = getVariableField(model, name, varargin)
            
            switch(lower(name))
              case {'sgmax'}
                index = 1;
                fn = 'sGmax';
              case {'rs'}
                index = 1;
                fn = 'rs';
              otherwise
                [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
            end
        end

% ------------------------------------------------------------------------        
        
       function names = getComponentNames(model)
           names = cellfun(@(c) c.name, model.Components, 'uniformoutput', false);    
           %names = getComponentNames@ReservoirModel(model);
        end

% ------------------------------------------------------------------------        
    
        function [vars, names, origin] = getPrimaryVariables(model, state)
    
            [p, s] = model.getProps(state, 'pressure', 'saturation');
            
            sW = s(:,1); % water saturation.  s(:,2) is CO2 saturation
            
            vars = {p, sW};
            names = {'pressure', 'sW'};
            origin = {class(model), class(model)};
            
            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
        end

% ------------------------------------------------------------------------        
        
        function [eqs, state, src] = ...
                addBoundaryConditionsAndSources(model, eqs, names, types, ...
                                                state, forces)
            
            [p, s, mob, rho] = model.getProps(state, 'PhasePressures', ...
                                                     's'             , ...
                                                     'Mobility'     , ...
                                                     'Density');
            dissolved = {};
            comps = model.Components; 
                        
            [eqs, state, src] = addBoundaryConditionsAndSources@ReservoirModel(...
                model, eqs, names, types, state, p, s, mob, rho, dissolved, comps, forces);
        end

% ------------------------------------------------------------------------        
        
        function ctrl = validateDrivingForces(model, ctrl, index)
            ctrl = validateDrivingForces@ReservoirModel(model, ctrl, index);
        end

% ------------------------------------------------------------------------        
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ReservoirModel(model);
        end

% ------------------------------------------------------------------------        
        
        function state = initStateAD(model, state, vars, names, origin)

            toFacility = strcmp(origin, class(model.FacilityModel));
            
            
            state = model.FacilityModel.initStateAD(state, vars(toFacility), ...
                                                    names(toFacility), ...
                                                    origin(toFacility));

            sw_ix = strcmp(names, 'sW');
            sw = vars{sw_ix};
            sg = 1 - sw;
            state = model.setProp(state, 'saturation', {sw, sg});
            
            remaining = ~toFacility;
            remaining(sw_ix) = false;
            
            state = initStateAD@ReservoirModel(model, state, vars(remaining), ...
                                               names(remaining), ...
                                               origin(remaining));
        end

% ------------------------------------------------------------------------        
            
        function [v_eqs, tolerances, names] = ...
                getConvergenceValues(model, problem, varargin)
            [v_eqs, tolerances, names] = ...
                getConvergenceValues@ReservoirModel(model, problem, varargin{:});
        end

% ------------------------------------------------------------------------        
        
        function [model, state] = updateForChangedControls(model, state, forces)
            [model, state] = updateForChangedControls@ReservoirModel(model, state, forces);
        end

% ------------------------------------------------------------------------        

        function [eq, src] = ...
                addComponentContributions(model, cname, eq, component, src, force)     
            
            cnames = model.getComponentNames();
            ix = strcmpi(cnames, cname);

            qC = src.phaseMass{find(ix)}; % @@ Does not take dissolution into account.  Refine later.
            
            if ~isempty(src.mapping)
                qC = src.mapping*qC;
            end
            
            cells = src.sourceCells;            
            eq(cells) = eq(cells) - qC;
            
        end
    
% ----------------------------------------------------------------------------
   
        function gdz = getGravityGradient(model)
            s  = model.operators;
            Gt = model.G;
            gdz = model.gravity * s.Grad(Gt.cells.z);
        end
    end
    
end



% ---------------------------------------

% function model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, varargin)
% function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
% function [fn, index] = getVariableField(model, name, varargin)
% function [state, report] = updateState(model, state, problem, dx, drivingForces)
% function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
% function model = setupOperators(model, Gt, rock, varargin)
    
    
% function gdz = getGravityGradient(model)
% function rhoS = getSurfaceDensities(model)
% function [problem, state] = getAdjointEquations(model, state0, state, dt, drivingForces, varargin)
