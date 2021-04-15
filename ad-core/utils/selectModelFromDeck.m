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
% OPTIONAL PARAMETERS:
%    useLegacyModels - Whether or not to construct the original, monolithic
%                      physical models.  Stop-gap solution until all
%                      examples have been ported to the Generic* framework
%                      and only supported for three-phase black-oil
%                      w/polymer (`ThreePhaseBlackOilPolymerModel`).
%
%    Any - All other key/value pair arguments are passed directly on to the
%          model constructor.
%
% RETURNS:
%   model      - Subclass of `PhysicalModel` approprioate for passing along
%                to `simulateScheduleAD`.
%
% SEE ALSO:
%   `GenericBlackOilModel`, `ThreePhaseBlackOilModel`,
%   `TwoPhaseOilWaterModel`, `OilWaterPolymerModel`,
%   `ThreePhaseBlackOilPolymerModel`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('useLegacyModels', false);
    [opt, extra] = merge_options(opt, varargin{:});

    rs = deck.RUNSPEC;
    check = @(name) isfield(rs, upper(name)) && rs.(upper(name));

    hasgas  = check('gas');
    hasoil  = check('oil');
    haswat  = check('water');
    haspoly = check('polymer');
    hassurf = check('surfact');
    hascomp = check('comps');

    arg = [{G, [], fluid}, extra];
    if hascomp
        % Compositional model
        assert(~hassurf, 'Compositional model does not support surfactant');
        assert(~haspoly, 'Compositional model does not support polymer');
        eos = initDeckEOSModel(deck);
        model = GenericNaturalVariablesModel(arg{:}, eos, 'water', haswat);
    elseif haspoly || hassurf
        % Water-based EOR model (surfactant and/or polymer)
        mrstModule add ad-eor
        if ~opt.useLegacyModels
            model = GenericSurfactantPolymerModel(arg{:},...
                'water', haswat, 'oil', hasoil, 'gas', hasgas);
        elseif haspoly && hassurf
            % Surfactant-Polymer EOR
            assert(haswat, 'SurfactantPolymer model requires water phase to be present');
            model = ThreePhaseSurfactantPolymerModel(arg{:});
        elseif haspoly && ~hassurf
            % Polymer EOR
            assert(haswat, 'Polymer model requires water phase to be present');
            if hasgas
                model = ThreePhaseBlackOilPolymerModel(arg{:});
            else
                model = OilWaterPolymerModel(arg{:});
            end
        elseif hassurf && ~haspoly
            % Surfactant EOR
            assert(haswat, 'Surfactant model requires water phase to be present');
            if hasgas
                model = ThreePhaseBlackOilSurfactantModel(arg{:});
            else
                model = OilWaterSurfactantModel(arg{:});
            end
        end
    else
        % Blackoil three phase
        model = GenericBlackOilModel(arg{:}, 'water', haswat, 'oil', hasoil, 'gas', hasgas);
    end
    % Set blackoil specific features
    if isa(model, 'ThreePhaseBlackOilModel') && isfield(deck, 'RUNSPEC')
        if check('BLACKOIL')
            model.disgas = isfield(deck.PROPS, 'PVTO');
            model.vapoil = isfield(deck.PROPS, 'PVTG');
        end
        if isfield(deck.RUNSPEC, 'VAPOIL')
            model.vapoil = deck.RUNSPEC.VAPOIL;
        end
        if isfield(deck.RUNSPEC, 'DISGAS')
            model.disgas = deck.RUNSPEC.DISGAS;
        end
    end
    model.inputdata = deck;
    if ~isempty(rock)
        model = model.setupOperators(G, rock);
    end
    model.rock = rock;
    model.FacilityModel = selectFacilityFromDeck(deck, model);
    model.AquiferModel  = selectAquiferFromDeck(deck, model);
    model.minimumPressure = 0;
end
