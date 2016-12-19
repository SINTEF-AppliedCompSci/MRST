classdef OilWaterPolymerModel < TwoPhaseOilWaterModel
% Oil/water/polymer system
%
%
% SYNOPSIS:
%   model = OilWaterPolymerModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Two phase model with polymer. A description of the polymer model
% that is implemented here can be found in the directory ad-eor/docs .
%
% PARAMETERS:
%   G        - Grid
%   rock     - Rock structure
%   fluid    - Fluid structure
%   varargin - optional parameter
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: ThreePhaseBlackOilPolymerModel
%

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
        % Polymer present
        polymer
        
        % Polymer differene tolerance
        useIncPolymerConvergence
        toleranceIncPolymer
        
        % Add extra output to wellsol/states for polymer quantities
        extraPolymerOutput
        
    end
    
    methods
        function model = OilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            % This is the model parameters for oil/water/polymer
            model.polymer = true;
            
            % Tolerance for the change in polymer concentration
            model.useIncPolymerConvergence = true;
            model.toleranceIncPolymer = 1e-3;
            
            model.extraPolymerOutput = false;
            
            model.wellVarNames = {'qWs', 'qOs', 'qWPoly', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end
        
        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'Polymer', model.G.cells.num, 1);
        end

        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            
            if model.polymer
                % Store the polymer from previous iteration temporarily to
                % use in convergence criteria
                c_prev = model.getProp(state, 'polymer');
            end
            
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);
            
            if model.polymer
                % Limit polymer concentration to [0, fluid.cmax]
                c = model.getProp(state, 'polymer');
                c = min(c, model.fluid.cmax);
                state = model.setProp(state, 'polymer', max(c, 0) );
                state.c_prev = c_prev;
                
                % Shear Thinning Report               
                % We (may) have stored the shear thinning report
                % temporarily in the state structure. We move this over to
                % the report structure instead. The reason for this is that
                % there is no report returned from the equations.
                if isfield(state, 'ShearThinningReport')
                    report.ShearThinning = state.ShearThinningReport;
                    state = rmfield(state, 'ShearThinningReport');
                end
            end
        end
        
        function [convergence, values, names] = checkConvergence(model, ...
                problem, varargin)
            
            if model.useIncPolymerConvergence
                polyEqnInx = find(problem.indexOfEquationName('polymer'));
                if polyEqnInx
                    % The convergence of polymer equation is checked below. In
                    % order to check the remaining equations, the polymer is
                    % removed from the problem before the parent is called.
                    problem_org = problem;
                    problem.equations(polyEqnInx) = [];
                    problem.types(polyEqnInx) = [];
                    problem.equationNames(polyEqnInx) = [];
                    problem.primaryVariables(polyEqnInx) = [];
                end
            end
            
            % Check convergence of all equations except the polymer eqn
            [convergence, values, names] = ...
                checkConvergence@TwoPhaseOilWaterModel(model, ...
                problem, varargin{:});
            
            if model.useIncPolymerConvergence
                if polyEqnInx
                    problem = problem_org;

                    % Compute the convergence norm for polymer
                    polyNorm = Inf;
                    if problem.iterationNo > 1
                        % Check polymer change from previous iteration
                        polyNorm = norm(problem.state.c - ...
                            problem.state.c_prev, Inf) / model.fluid.cmax;
                        convergence = convergence && ...
                            polyNorm < model.toleranceIncPolymer;
                    end

                    % In the printed convergence information (in verbose mode),
                    % we insert the polymer residual after oil and water.
                    nwo = sum([model.water model.oil]);
                    if model.useCNVConvergence
                        inx = 2*nwo + 1; % after water and oil, both CNV and MB
                    else
                        inx = nwo + 1; % after water and oil
                    end

                    % Insert polymer data into values
                    values = [values(1:inx-1)  polyNorm   values(inx:end)];
                    names  = [ names(1:inx-1) {'poly (cell)'}  names(inx:end)];
                end
            end
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, drivingForces);
            if model.polymer
                c     = model.getProp(state, 'polymer');
                cmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cmax, c));
                
                if isfield(state, 'c_prev')
                    % Remove the temporary field used for convergence
                    state = rmfield(state, 'c_prev');
                end
            end
        end

        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'polymer'}
                    % The current polymer concentration in a cell
                    index = 1;
                    fn = 'c';
                case {'polymermax'}
                    % The maximum polymer conctration throughout the
                    % simulation for each cell
                    index = 1;
                    fn = 'cmax';
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(...
                                    model, name);
            end
        end
        
        function scaling = getScalingFactorsCPR(model, problem, names)
            nNames = numel(names);

            scaling = cell(nNames, 1);
            handled = false(nNames, 1);

            for iter = 1:nNames
                name = lower(names{iter});
                switch name
                    case 'polymer'
                        s = 0;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);

                scaling{iter} = s;
                handled(sub) = true;
            end
            if ~all(handled)
                % Get rest of scaling factors
                other = getScalingFactorsCPR@ThreePhaseBlackOilModel(model, problem, names(~handled));
                [scaling{~handled}] = other{:};
            end
        end
        
        
        %------------------------------------------------------------------
        % FUNCTIONS FOR STORING STATE DATA
        %------------------------------------------------------------------
        
        function state = storeFluxes(model, state, vW, vO, vP)
            % Utility function for storing the interface fluxes in the
            % state
            internal = model.operators.internalConn;
            state.flux = zeros(numel(internal), 3);
            state.flux(internal,:) = [double(vW), double(vO), double(vP)];
        end
        
        function state = storeMobilities(model, state, mobW, mobO, mobP) %#ok<INUSL>
            % Utility function for storing the mobilities in the state
            state.mob = [double(mobW), double(mobO), double(mobP)];
        end
        
        function state = storeShearMultiplier(model, state, ...
                shearMult) %#ok<INUSL>
            % Utility function for storing the polymer shear thinning /
            % thickening velocity multiplier in the state. Note that this
            % is the reciprocal of the viscosity multiplier.
            state.shearMult = double(shearMult);
        end
        
        function state = storeEffectiveWaterVisc(model, state, ...
                muWeff) %#ok<INUSL>
            % Utility function for storing the effective water viscosity
            % due to the presence of polymer.
            state.muWeff = double(muWeff);
        end
        
        function state = storeEffectivePolymerVisc(model, state, ...
                muPeff) %#ok<INUSL>
            % Utility function for storing the effective water viscosity
            % due to the presence of polymer.
            state.muPeff = double(muPeff);
        end
        
        function state = storePolymerAdsorption(model, state, ...
                ads) %#ok<INUSL>
            % Utility function for storing the polymer adsorption.
            state.ads = double(ads);
        end
        
        function state = storeRelpermReductionFactor(model, state, ...
                Rk) %#ok<INUSL>
            % Utility function for storing the relative permeability
            % reduction factor due to the presence of polymer.
            state.Rk = double(Rk);
        end
        
    end
end

