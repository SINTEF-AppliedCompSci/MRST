classdef OilWaterPolymerModel < TwoPhaseOilWaterModel
    % Oil/water/polymer system
    % This model is a two phase oil/water model, extended with the polymer
    % phase in addition.
    
    properties
        % Polymer present
        polymer
        
        % Polymer differene tolerence
        tolerancePolymer
        
        % Add extra output to wellsol/states for polymer quantities
        extraPolymerOutput
        
    end
    
    methods
        function model = OilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            % This is the model parameters for oil/water/polymer
            model.polymer = true;
            
            % Tolerance for the change in polymer concentration
            model.tolerancePolymer = 1e-3;
            
            model.extraPolymerOutput = false;
            
            model.wellVarNames = {'qWs', 'qOs', 'qWPoly', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end
        
        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            
            if model.polymer
                % Store the polymer from previous iteration to use in
                % convergence criteria
                c_prev = model.getProp(state, 'polymer');
            end
            
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);
            
            if model.polymer
                c = model.getProp(state, 'polymer');
                c = min(c, model.fluid.cmax);
                state = model.setProp(state, 'polymer', max(c, 0) );
                state.c_prev = c_prev;
            end
        end
        
        function [convergence, values] = checkConvergence(model, ...
                problem, varargin)
            
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
                if model.useCNVConvergence
                    % Hack to print CNV convergence including polymer
                    wasVerbose  = mrstVerbose();
                    mrstVerbose(false);
                end
            end
            
            % Check convergence of all equations except the polymer eqn
            [convergence, values] = ...
                checkConvergence@TwoPhaseOilWaterModel(model, ...
                problem, varargin{:});
            
            if polyEqnInx
                problem = problem_org;
                polyStr = sprintf('-\t');
                if problem.iterationNo > 1
                    % Check polymer change from previous iteration
                    polyNorm = norm(problem.state.c - ...
                        problem.state.c_prev, Inf) / model.fluid.cmax;
                    convergence = convergence && ...
                        polyNorm < model.tolerancePolymer;
                    polyStr = sprintf('%2.2e', polyNorm);
                end
                if model.useCNVConvergence
                    % Hack to print CNV convergence including polymer
                    mrstVerbose(wasVerbose);
                    if mrstVerbose()
                        pinx  = [model.oil, model.water];
                        np    = 2*sum(pinx);
                        if problem.iterationNo == 1
                            ptext = {'CNVO','CNVW','MBO','MBW'};
                            winx  = problem.indexOfType('perf') | ...
                                    problem.indexOfType('well');
                            wtext = problem.equationNames(winx);
                            wtext = strrep(wtext, 'waterWells', 'WATERW');
                            wtext = strrep(wtext, 'oilWells', 'OILW');
                            wtext = strrep(wtext, 'polymerWells', 'POLYW');
                            wtext = strrep(wtext, 'closureWells', 'CLOSE');
                            text  = [ptext([pinx pinx]) {'POLY'} wtext];
                            fprintf('%s\n', sprintf('%s\t\t',text{:}) );
                        end
                        fprintf('%2.2e\t', values(1:np));
                        fprintf('%s\t',    polyStr);
                        fprintf('%2.2e\t', values(np+1:end));
                        fprintf('\n')
                    end
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
                    index = 1;
                    fn = 'c';
                case {'polymermax'}
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
        
        function state = storeFluxes(model, state, vW, vO, vP) %#ok<INUSL>
            % Utility function for storing the interface fluxes in the
            % state
            state.flux = [double(vW), double(vO), double(vP)];
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

