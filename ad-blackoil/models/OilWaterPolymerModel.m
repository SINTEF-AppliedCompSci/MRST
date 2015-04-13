classdef OilWaterPolymerModel < TwoPhaseOilWaterModel
    % Oil/water/polymer system
    % This model is a two phase oil/water model, extended with the polymer
    % phase in addition.
    
    properties
        % Polymer present
        polymer
        
        % Polymer differene tolerence
        tolerancePolymer
        
    end
    
    methods
        function model = OilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            % This is the model parameters for oil/water/polymer
            model.polymer = true;
            
            % Tolerance for the change in polymer concentration
            model.tolerancePolymer = model.toleranceCNV;
            
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
                polyNorm = Inf;
                if problem.iterationNo > 1
                    % Check polymer change from previous iteration
                    polyNorm = norm(problem.state.c - ...
                        problem.state.c_prev, Inf) / model.fluid.cmax;
                    convergence = convergence && ...
                        polyNorm < model.tolerancePolymer;
                end
                if model.useCNVConvergence
                    % Hack to print CNV convergence including polymer
                    mrstVerbose(wasVerbose);
                    if mrstVerbose()
                        inx = [model.oil, model.water];
                        if problem.iterationNo == 1
                            text = {'CNVO','CNVW','MBO','MBW','POLY'};
                            text = text([inx inx model.polymer]);
                            fprintf('%s\n', sprintf('%s\t\t',text{:}) );
                        end
                        fprintf('%2.2e\t', [values polyNorm]);
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

    end
end
