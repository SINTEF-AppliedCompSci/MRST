classdef OilWaterSurfactantModel < TwoPhaseOilWaterModel
%
%
% SYNOPSIS:
%   model = FullyImplicitOilWaterSurfactantModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Fully implicit model for an oil water system with surfactant. All
% the equations are solved implicitly. A description of the surfactant model
% that is implemented can be found in the directory ad-eor/docs .
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
% SEE ALSO: equationsOilWaterSurfactant, ImplicitExplicitOilWaterSurfactantModel
%



    properties
        surfactant
    end

    methods

        function model = OilWaterSurfactantModel(G, rock, fluid, varargin)

            model = model@TwoPhaseOilWaterModel(G, rock, fluid, varargin{:});
            model = model.setupOperators(G, rock, varargin{:});
            model.surfactant = true;
            model.wellVarNames = {'qWs', 'qOs', 'qWSft', 'bhp'};
            model = merge_options(model, varargin{:});

        end

        function model = setupOperators(model, G, rock, varargin)
            model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
            model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, ...
                                                           varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@OilWaterSurfactantBaseModel(model, state, problem,  dx, ...
                                                              drivingForces);
            % cap the concentration (only if implicit solver for concentration)
            c = model.getProp(state, 'surfactant');
            state = model.setProp(state, 'surfactant', max(c, 0) );

        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, ...
                                                              drivingForces);
            state = updateAdsorption(state0, state, model);

        end
        
        function varargout = evaluateRelPerm(model, sat, varargin)
            error('function evaluateRelPerm is not implemented for surfactant model')
        end

        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            nc = model.G.cells.num;
            model.checkProperty(state, 'Surfactant', [nc, 1], [1, 2]);
            model.checkProperty(state, 'SurfactantMax', [nc, 1], [1, 2]);
        end

        function [fn, index] = getVariableField(model, name)
        % Get the index/name mapping for the model (such as where
        % pressure or water saturation is located in state)
            switch(lower(name))
              case {'ads'} % needed when model.explicitAdsorption
                index = 1;
                fn = 'ads';
              case {'adsmax'} % needed when model.explicitAdsorption
                index = 1;
                fn = 'adsmax';
              case {'surfactant'}
                index = 1;
                fn = 'c';
              case {'surfactantmax'}
                index = 1;
                fn = 'cmax';
              otherwise
                [fn, index] = getVariableField@TwoPhaseOilWaterModel(...
                    model, name);
            end
        end

        function state = storeSurfData(model, state, s, c, Nc, sigma)
            state.SWAT    = double(s);
            state.SURFACT = double(c);
            state.SURFCNM = log(double(Nc))/log(10);
            state.SURFST  = double(sigma);
            % state.SURFADS = double(ads);
        end


    end
end

