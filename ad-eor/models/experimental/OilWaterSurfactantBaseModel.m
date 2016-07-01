classdef OilWaterSurfactantBaseModel < TwoPhaseOilWaterModel
%
%
% SYNOPSIS:
%   model = OilWaterSurfactantBaseModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Base class for model with oil, water and surfactant. Not
% instanciated directly but used as a parent class.
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
% SEE ALSO: ImplicitExplicitOilWaterSurfactantModel, FullyImplicitOilWaterSurfactantModel
%


    properties
        surfactant
    end

    methods

        function model = OilWaterSurfactantBaseModel(G, rock, fluid, varargin)

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
