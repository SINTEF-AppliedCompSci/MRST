classdef FullyImplicitOilWaterSurfactantModel < OilWaterSurfactantBaseModel
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

    methods

        function model = FullyImplicitOilWaterSurfactantModel(G, rock, fluid, varargin)
            model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
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

    end
end
