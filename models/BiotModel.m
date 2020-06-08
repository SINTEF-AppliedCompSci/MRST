classdef BiotModel < GenericBlackOilModel

    properties
        
        % Structure that contains the mechanical parameter of the system
        % This structure should contain the fields:
        %
        % E     - Young's modulus (one entry per cell)
        % nu    - Poisson ratio (one entry per cell)
        % el_bc - Structure describing the boundary condition (see VEM_linElast)
        % load  - Structure giving the volumetric forces (see VEM_linElast)
        %
        mech;
        
        % Property container
        MechBiotPropertyFunctions;

    end

    methods
        
        function model = BiotModel(G, rock, fluid, mech_problem, varargin)
            
            model = model@GenericBlackOilModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
            model.OutputStateFunctions = {'Dilatation'};
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end
            % Physical properties of rock and fluid
            model.mech  = mech_problem;

            % Compute stiffness tensor C, if not given
            if ~isfield(model.mech, 'C')
                [model.mech.C, model.mech.invC, model.mech.invCi] = ...
                    Enu2C(model.mech.E, model.mech.nu, model.G);
            end
            
            alpha_scaling = 1;  % default values
            S             = []; % default values
            
            % Add mechanical operators  
            operators = model.operators;
            mechOperators = setupOperatorsVEM(model.G, ...
                                              model.mech.C, ...
                                              model.mech.el_bc, ...
                                              model.mech.load, ...
                                              alpha_scaling, S);

            operators.A             = mechOperators.A;
            operators.rhs           = mechOperators.rhs;
            operators.trace         = mechOperators.trace;
            operators.isdirdofs     = mechOperators.isdirdofs;
            operators.V_dir         = mechOperators.V_dir;
            operators.ovol_div      = mechOperators.ovol_div;
            operators.global_strain = mechOperators.global_strain;
            operators.gradP         = mechOperators.gradP;
            
            model.operators = operators;
            
            % Add mechanical properties
            model.MechBiotPropertyFunctions = MechBiotPropertyFunctions(model);
            
        end
        

        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@GenericBlackOilModel(model);
            containers = [containers, {model.MechBiotPropertyFunctions}];
        end

        function model = setupStateFunctionGroupings(model, varargin) 
            
            model = setupStateFunctionGroupings@GenericBlackOilModel(model, varargin{:});
            
            % Use mechanical coupling in flow equations
            fp = model.FlowPropertyFunctions;
            fp = fp.setStateFunction('BasePoreVolume', BlackOilPoreVolume(model));
            fp = fp.setStateFunction('PoreVolume'    , BiotPoreVolume(model));
            fp = fp.setStateFunction('Dilatation'    , BiotCoupledDilatation(model));
            model.FlowPropertyFunctions = fp;
            
            % Use fluid coupling in mechanical equations
            mp = model.MechBiotPropertyFunctions;
            mp = mp.setStateFunction('BiotGradP', BiotCoupledGradP(model));
            model.MechBiotPropertyFunctions = mp;
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, ...
                                                                       state, ...
                                                                       dt, ...
                                                                       drivingForces, varargin)
            
            % fluid equations
            [fluideqs, flux, fluidnames, fluidtypes] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            
            % Treat source or bc terms
            if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
                [pressures, sat, mob, rho, rs, rv] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'Rs', 'Rv');
                dissolved = model.getDissolutionMatrix(rs, rv);
                fluideqs = model.addBoundaryConditionsAndSources(fluideqs, fluidnames, fluidtypes, state, ...
                                                                 pressures, sat, mob, rho, ...
                                                                 dissolved, {}, ...
                                                                 drivingForces);
            end
            
            % Add aquifer contributions if any.
            if ~isempty(model.AquiferModel)
                fluideqs = addAquifersContribution(model.AquiferModel, fluideqs, fluidnames, state, dt);
            end
            
            % Add sources
            fluideqs = model.insertSources(fluideqs, src);
            
            % Assemble equations
            for i = 1:numel(fluideqs)
                fluideqs{i} = model.operators.AccDiv(fluideqs{i}, flux{i});
            end
            
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
  
            
            % The momentum equations
            
            mecheqs = model.getProp(state, 'MomentumEquations');
            
            eqs = [{mecheqs}, fluideqs, weqs];
            names = [{'momeqs'}, fluidnames, wnames];
            types = [{'nodes'}, fluidtypes, wtypes];
            
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            [vars, names, origin] = getPrimaryVariables@GenericBlackOilModel(model, state);
            xd = model.getProp(state, 'xd');
            vars = [vars, {xd}];
            names = [names, {'xd'}];
            origin = [origin, {class(model)}];
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            
            switch(lower(name))
              case {'xd'}
                % same as 'u' but the degree of freedom where the Dirichlet conditions (fixed
                % displacement) are removed
                fn = 'xd';
                index = ':';
              otherwise
                [fn, index] = getVariableField@GenericBlackOilModel(model, name, varargin{:});
            end

        end

    end

end
