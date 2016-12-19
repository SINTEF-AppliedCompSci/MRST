classdef PressureSaturationSurfactantModel < OilWaterSurfactantBaseModel
%
%
% SYNOPSIS:
%   model = PressureSaturationSurfactantModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Model used to solve the pressure-saturation equation implicitly
% for a constant (in time) surfactant concentration.
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
% SEE ALSO: ImplicitExplicitOilWaterSurfactantModel, equationsPressureSaturationForOilWaterSurfactant
%


    methods
        function model = PressureSaturationSurfactantModel(G, rock, fluid, varargin)
            model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsPressureSaturationForOilWaterSurfactant(state0, state, model, ...
                                                              dt, drivingForces, varargin{:});
        end


    end
end
