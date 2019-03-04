classdef TransportOilWaterPolymerModelDG < TransportOilWaterModelDG
   
    properties
        polymer
    end
    
    methods
      function model = TransportOilWaterPolymerModelDG(G, rock, fluid, varargin)
            model = model@TransportOilWaterModelDG(G, rock, fluid, varargin{:});
            model.oil     = true;
            model.water   = true;
            model.polymer = true;
            model.conserveWater = false;
            model.conserveOil   = true;
        end

        % ----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] ...
                = transportEquationOilWaterPolymerDG(state0, state, model, dt, drivingForces, ...
                                  'solveForOil'  , model.conserveOil  , ...
                                  'solveForWater', model.conserveWater, ...
                                  varargin{:}                         );
        end
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch lower(name)
                case {'cdof'}
                    index = 1;
                    fn = 'cdof';
                case {'c', 'polymer', 'polymermax'}
                    c = model.getComponentNames();
                    index = find(strcmpi(c, 'polymer'));
                    if any(strcmpi(name, {'polymer', 'c'}))
                        fn = 'c';
                    else
                        fn = 'cmax';
                    end
                otherwise
                    [fn, index] = getVariableField@TransportOilWaterModelDG(model, name);
            end
        end
        % ----------------------------------------------------------------%
        function names = getComponentNames(model)
            names = getComponentNames@TransportOilWaterModelDG(model);
            if model.polymer
                names{end+1} = 'polymer';
            end
        end
        % ----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            if model.polymer
                % Store the polymer from previous iteration temporarily to
                % use in convergence criteria
                cdof_prev = model.getProp(state, 'polymer');
            end

            [state, report] = updateState@TransportOilWaterModelDG(model, ...
               state, problem, dx, drivingForces);

            if model.polymer
                
                state = model.capProperty(state, 'c', 0, model.fluid.cmax);
%                 % Limit polymer concentration to [0, fluid.cmax]
%                 c = model.getProp(state, 'polymer');
%                 c = min(c, model.fluid.cmax);
%                 state = model.setProp(state, 'polymer', max(c, 0) );
                state.c_prev = cdof_prev;

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
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TransportOilWaterModelDG(model, state0, state, dt, drivingForces);
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

    end
    
end