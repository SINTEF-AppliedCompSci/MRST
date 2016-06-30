classdef ThreePhaseBlackOilPolymerModel < ThreePhaseBlackOilModel
    % ThreePhaseBlackOil/polymer system
    % This model is a three phase black oil model, extended with
    % a polymer component.

    properties
        % Polymer present
        polymer
        % Using PLYSHEAR shear model based on water velocity
        usingShear
	    % Using PLYSHLOG shear model based on water velocity
        usingShearLog
	    % Using PLYSHLOG shear model base on water shear rate
        usingShearLogshrate
    end

    methods
        function model = ThreePhaseBlackOilPolymerModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

            % This is the model parameters for oil/water/gas/polymer system
            model.polymer = true;

            if (isfield(fluid,'shrate') && ~isfield(fluid,'plyshlog'))
                error('SHRATE is specified while PLYSHLOG is not specified')
            end
            if (isfield(fluid, 'plyshearMult') && isfield(fluid, 'plyshlog'))
                error('PLYSHLOG and PLYSHEAR are existing together');
            end

            model.usingShear = isfield(fluid, 'plyshearMult');
            model.usingShearLog = (isfield(fluid, 'plyshlog') && ~isfield(fluid, 'shrate'));
            model.usingShearLogshrate = (isfield(fluid, 'plyshlog') && isfield(fluid, 'shrate'));
            model.wellVarNames = {'qWs', 'qOs', 'qGs', 'qWPoly', 'bhp'};
            model = merge_options(model, varargin{:});
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsThreePhaseBlackOilPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end

        % --------------------------------------------------------------------%
        function state = validateState(model, state)
            state = validateState@ThreePhaseBlackOilModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'Polymer', [model.G.cells.num, 1], [1, 2]);
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@ThreePhaseBlackOilModel(model, ...
               state, problem,  dx, drivingForces);

            if model.polymer
                c = model.getProp(state, 'polymer');
                c = min(c, model.fluid.cmax);
                state = model.setProp(state, 'polymer', max(c, 0));
            end
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ThreePhaseBlackOilModel(model, state0, state, dt, drivingForces);
            if model.polymer
                c     = model.getProp(state, 'polymer');
                cmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cmax, c));
            end
        end

        % --------------------------------------------------------------------%
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
                    [fn, index] = getVariableField@ThreePhaseBlackOilModel(...
                                    model, name);
            end
        end

        % --------------------------------------------------------------------%
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
