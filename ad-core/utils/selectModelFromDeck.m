function model = selectModelFromDeck(G, rock, fluid, deck, varargin)
%Select simulation model from a ECLIPSE/FrontSim style input deck
%
% SYNOPSIS:
%   model = selectModelFromDeck(G, rock, fluid, deck)
%
% DESCRIPTION:
%   Determine the type of PhysicalModel subclass (if any) most suitable for
%   simulating a given input deck.
%
% REQUIRED PARAMETERS:
%
%   G          - Simulation grid, typically from initEclipseGrid
%
%   rock       - Corresponding rock structure, typically from
%                initEclipseRock.
%
%   fluid      - Fluid model, typically from initDeckADIFluid.
%
%   deck       - Parsed input deck, typically from readEclipseDeck.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%    Are passed onto the model constructor directly.
%
% RETURNS:
%   model      - Subclass of PhysicalModel approprioate for passing along
%                to simulateScheduleAD.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, OilWaterPolymerModel

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

    require ad-blackoil

    rs = deck.RUNSPEC;
    check = @(name) isfield(rs, upper(name)) && rs.(upper(name));

    hasgas  = check('gas');
    hasoil  = check('oil');
    haswat  = check('water');
    haspoly = check('polymer');
    hassurf = check('surfact');
    
    if haspoly
        % Polymer EOR
        assert(~hassurf, 'Polymer model does not support surfactant');
        assert(haswat, 'Polymer model requires water phase to be present');
        if hasgas
            model = ThreePhaseBlackOilPolymerModel(G, rock, fluid, 'inputdata', deck, varargin{:});
        else
            model = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck, varargin{:});
        end
    elseif hassurf
        % Surfactant EOR
        assert(~hasgas, 'Surfactant model does not support gas');
        assert(~haspoly, 'Surfactant model does not support polymer');
        assert(haswat, 'Surfactant model requires water phase to be present');
        model = OilWaterSurfactantModel(G, rock, fluid, 'inputdata', deck, varargin{:});
    elseif hasgas && hasoil && haswat
        % Blackoil three phase
        model = ThreePhaseBlackOilModel(G, rock, fluid, 'inputdata', deck, varargin{:});
    elseif hasoil && haswat
        % Two-phase oil-water
        model = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck, varargin{:});
    elseif haswat && hasgas
        % Two-phase water-gas (Note: Model tailored for CO2 storage, uses
        % CO2lab module).
        require co2lab
        model = twoPhaseGasWaterModel(G, rock, fluid, nan, nan, 'inputdata', deck, varargin{:});
    else
        error('Did not find matching model');
    end
    model.FacilityModel = selectFacilityFromDeck(deck, model);
end
