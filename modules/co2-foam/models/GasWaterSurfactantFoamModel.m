classdef GasWaterSurfactantFoamModel < ThreePhaseBlackOilModel
%Model for simulation of CO2 storage with mobility control, modelled as an
%immiscible system with gas, water, and surfactant foam, where the
%surfactant may partition into the gas phase, water phase or both
%
% SYNOPSIS:
%   model = GasWaterSurfactantFoamModel(G, rock, fluid, varargin)
%
% DESCRIPTION: 
%   Fully implicit model for an oil water system with surfactant. All
%   the equations are solved implicitly. A description of the surfactant model
%   that is implemented can be found in the directory ad-eor/docs .
%
% PARAMETERS:
%   G        - Grid
%   rock     - Rock structure
%   fluid    - Fluid structure
%   varargin - optional model parameters
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: equationsGasWaterSurfactantFoam,
%

%{ 
Copyright 2009-2023 SINTEF
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

    properties
        CO2foam
        surfactantmax
        % Use state from previous time step to evaluate mobility reduction
        useLaggedMobMult
        % UNUSED (No LAGGING STRATEGIES NOW)
        laggedExtra
        % Use mimetic method when mapping face values to cell values
        useMimeticMapping
        % disgas is a property of ThreePhaseBlackOilModel that handles gas dissolved in oil.
		% We would like to use something similar for the calculation of CO2 dissolved in water,
		% but it is probably best to make this a special feature of the CO2foam model, to avoid
		% complications. Thus a property to inidcate whether CO2 dissolution should be included.
		% We will use the RS-formulation and see if the ThreePhaseBlackOilModel functions are good.
		% See also the CO2lab module.
		disCO2
        sfwat
        sfgas
    end

    methods

        %-----------------------------------------------------------------%
        function model = GasWaterSurfactantFoamModel(G, rock, fluid, varargin)

            % Call parent constructor and set up operators
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});
            model = model.setupOperators(G, rock, varargin{:});

            % Set active phases
            model.oil   = false;
            model.gas   = true;
            model.water = true;
            
            % Set foam-specific proerties
            model.CO2foam           = true;
            model.useLaggedMobMult  = false; % @@
            model.laggedExtra       = {};
            model.useMimeticMapping = false;
			model.disCO2            = false;
            
            % Set nonlinear tolerance (used to check convergence of scaled
            % foam equation residual) @@ Check that it is indeed the case
            model.nonlinearTolerance = 1e-3;
            
            % Owervrite with user-defined settings
            model = merge_options(model, varargin{:});

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            
            % Call parent setup
            model = setupOperators@ThreePhaseBlackOilModel( ...
                model, model.G, model.rock, varargin{:});
            % Set up operators for computing velocity
            % @@ Figure out where this is used
            model.operators.veloc = computeVelocTPFA( ...
                model.G, model.operators.internalConn);
            model.operators.sqVeloc = computeSqVelocTPFA( ...
                model.G, model.operators.internalConn);
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        % Get model equations
            
            [problem, state] = equationsGasWaterSurfactantFoam( ...
                state0, state, model, dt, drivingForces, varargin{:});
                                                       
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
        % Update state from given Newton increment
        
            [state, report] = updateState@ThreePhaseBlackOilModel( ...
                model, state, problem,  dx, drivingForces);
            % cap the concentration (only if implicit solver for concentration)
            if model.CO2foam
                cs    = model.getProp(state, 'foam');
                state = model.setProp(state, 'foam', max(cs, 0) );
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % Update state after timstep convergence

            [state, report] = updateAfterConvergence@ThreePhaseBlackOilModel( ...
                model, state0, state, dt, drivingForces);
            
              if model.CO2foam
                  % Set maximum foam concentration
                  cs    = model.getProp(state, 'foam');
                  csmax = model.getProp(state, 'fmax');
                  state = model.setProp(state, 'fmax', max(csmax, cs));
              end
              
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            
            if isempty(model.FlowPropertyFunctions)
                model.FlowPropertyFunctions = FoamFlowPropertyFunctions(model);
                % TODO
                % AAG: Add special cases e.g. for no velocity dependence (no shear).
                % @@ This should be filled out or removed
            end
            if isempty(model.FlowDiscretization)
                model.FlowDiscretization = FlowDiscretization(model);
            end
            model = validateModel@ThreePhaseBlackOilModel(model, varargin{:});
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function state = validateState(model, state)
        % Validate state before simulation
            
            % Parent model validates almost everything
            state = validateState@ThreePhaseBlackOilModel(model, state);
            % Check foam properties
            nc = model.G.cells.num;
            model.checkProperty(state, 'foam' , [nc, 1], [1, 2]);
            model.checkProperty(state, 'fmax' , [nc, 1], [1, 2]);
            model.checkProperty(state, 'sfwat', [nc, 1], [1, 2]);
            model.checkProperty(state, 'sfgas', [nc, 1], [1, 2]);
            model.checkProperty(state, 'sfads', [nc, 1], [1, 2]);
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            
            switch(lower(name))
                case {'foam'}
                % Total surfactant concentration (kg/m^3 pore volume)
                    index = 1;
                    fn    = 'cs';
                case {'fmax'}
                % Maximum historical value for cs (@@ add unit)
                    index = ':';
                    fn    = 'csmax';
                case 'qwsft'
                % @@ Add description
                    index = ':';
                    fn    = 'qWSft';
                case {'sfwat'}
                % Surfactant concentration in water (mass fraction)
                    index = 1;
                    fn    = 'cW';
                case {'sfgas'}
                % Surfactant concentration in gas (mass fraction)
                    index = 1;
                    fn    = 'cG';
                case {'sfads'}
                % Adsorbed surfactant concentration (mass fraction)
                    index = 1;
                    fn    = 'cA';
                case {'sfadsmax'}
                % Maximum historical value for cA
                    index = 1;
                    fn    = 'cAmax';
                otherwise
                % Variable is not specific to the model, query parent model
                    [fn, index] = getVariableField@ThreePhaseBlackOilModel(...
                        model, name, varargin{:});
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
            
            names = getComponentNames@ThreePhaseBlackOilModel(model);
            if model.CO2foam
                names{end+1} = 'foam';
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function state = storeSurfData(model, state,  cs, cW, cG, cA, cAmax)
        % @@ Fiugre out if this is really needed
            
            state.SURFACT = double(cs);
            state.SFWAT   = double(cW);
            state.SFGAS   = double(cG);
            state.SFADS   = double(cA);
            %state.SFADSMAX = double(max(cAmax, cs)); % Error comparing variables with different units.
            state.SFADSMAX = double(cAmax);

        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [names, types] = getExtraWellEquationNames(model)
            
            [names, types] = getExtraWellEquationNames@ThreePhaseBlackOilModel(model);
            if model.CO2foam
                 names{end+1} = 'surfactantWells';
                 types{end+1} = 'perf';
            end
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function names = getExtraWellPrimaryVariableNames(model)
            
            names = getExtraWellPrimaryVariableNames@ThreePhaseBlackOilModel(model);
            if model.CO2foam
                  names{end+1} = 'qWSft';
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % PARAMETERS:
        %
        %   model  - (Base class, automatic)
        %
        %   cname  - Name of the component. Must be a property known to the
        %            model itself through `getProp` and `getVariableField`.
        %
        %   eq     - Equation where the source terms are to be added. Should
        %            be one value per cell in the simulation grid (model.G)
        %            so that the src.sourceCells is meaningful.
        %
        %   component - Cell-wise values of the component in question. Used
        %               for outflow source terms only.
        %
        %   src    - Source struct containing fields for fluxes etc. Should
        %            be constructed from force and the current reservoir
        %            state by `computeSourcesAndBoundaryConditionsAD`.
        %
        %   force  - Force struct used to produce src. Should contain the
        %            field defining the component in question, so that the
        %            inflow of the component through the boundary condition
        %            or source terms can accurately by estimated.
        
        % @@ Figure out if this function is needed
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch lower(cname)
              case {'surfactant'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
              otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
            
        end
        %-----------------------------------------------------------------%

        %-----------------------------------------------------------------%
        function [compEqs, compSrc, compNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)

            [compEqs, compSrc, compNames, wellSol] ...
                = getExtraWellContributions@ThreePhaseBlackOilModel( ...
                model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.CO2foam
                % Implementation of surfactant source terms
                assert(model.water && model.gas, ...
                    'Foam model requires both water and gas phases.');
                if well.isInjector()
                    concWell = model.getProp(well.W, 'foam');
                    assert(sum(well.W.compi) == 1, ...
                        'Only one injected phase per well allowed');
                    % Check flowing phase in injection wells
                    % Mass flow rate surfactant
                    cqS = concWell./(1-concWell).*qMass{well.W.compi>0};
                else
                    cqWs = qMass{1}; % Mass flow rate water
                    cqGs = qMass{2}; % Mass flow rate gas
                    cW = packed.components{1}{1};
                    cG = packed.components{1}{2};
                    cqWatSurf = cW./(1-cW).*cqWs;
                    cqGasSurf = cG./(1-cG).*cqGs;
                    cqS = cqWatSurf + cqGasSurf; % Mass flow rate surfactant
                end
                qwsft = packed.extravars{strcmpi(packed.extravars_names, 'qwsft')};
                compEqs{end+1} = qwsft - sum(cqS);
                compSrc{end+1} = cqS;
                compNames{end+1} = 'surfactantWells';
                
            end
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [values, tolerances, names] = getConvergenceValues(model, problem, varargin)
        
            % Get convergence values from parent model (typically CNV/MB
            % for water/gas phases);
            [values, tolerances, names] ...
                = getConvergenceValues@ThreePhaseBlackOilModel( ...
                model, problem, varargin{:});
            
            % Scale foam residual by dt/(pv*rho)
            % Get indices
            cix = strcmpi(names, 'foam');
            eix = strcmpi(problem.equationNames, 'foam');
            
            % Set density equal to sum of adsorbed surfactant, and density
            % of all phases surfactant may be present in
            [SG, SW] = model.surfactantPartitioning();
            rho = model.fluid.rhoRSft;
            if SG, rho = rho + model.fluid.rhoGS; end
            if SW, rho = rho + model.fluid.rhoWS; end
            % Get pore volume and timestep size
            pv = model.operators.pv;
            dt = problem.dt;
            % Scale residual
            res = value(problem.equations{eix});
            res = res.*dt./(pv.*rho);
            % Set value and tolerance
            values(cix)     = norm(res, inf);
            tolerances(cix) = model.nonlinearTolerance;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function [SG, SW] = surfactantPartitioning(model)
        % Determine surfactant partitioning configuration
            
            fluid = model.fluid;
            if isfield(fluid,'Cpart')
                % Check if field is nonempty
                if ~isempty(fluid.Cpart)
                    % Partitioning constant is given
                    Cpart = fluid.Cpart;
                    assert(Cpart>= 0, ...
                        'Partitioning constant must be non-negative.')
                else
                    Cpart = NaN;
                end
            else
                % No partitioning field. Surfactant stays in the phase
                % given by the fluid model, i.e. by 
                % fluid.surfingas = true/false.
                Cpart = NaN; % Not used in this case.
            end
            % Partitioning works as: Cg/Cw=Cpart. Partitioning also needs
            % to be consistent with the adsorption, which is calculated
            % from the concentration in the water phase, except for the
            % case where surfactant is not soluble in water. (This happens
            % if Cpart is not defined, and fluid.surfingas = true.)

            % Resolve various possible configurations
            doPart = ~isnan(Cpart);
            % Conservation of mass of surfactant:
            % First check which phases can contain surfactant
            if ~doPart
                if fluid.surfingas
                    [SG, SW] = deal(true, false);
                else
                    [SG, SW] = deal(false, true); 
                end
            else
                [SG, SW] = deal(Cpart > 0, true);
            end
            
        end
        %-----------------------------------------------------------------%
    
    end
    
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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