classdef BiochemicalFlowDiscretization < FlowDiscretization
    % BiochemicalFlowDiscretization
    % Discretization and state function grouping for bio-chemical flow
    % within a compositional model with microbial growth.

    properties
        % No additional properties needed
    end

    methods
        %-----------------------------------------------------------------%
        function props = BiochemicalFlowDiscretization(model)
            % Constructor: inherit base FlowDiscretization properties
            props = props@FlowDiscretization(model);
            if model.bacteriamodel

                % Set up transmissibility and porosity functions
                if model.dynamicFlowTrans
                    % Dynamic transmissibility computed each nonlinear iteration
                    props = props.setStateFunction('Permeability', BactPermeability(model));
                    props = props.setStateFunction('Transmissibility', ...
                        DynamicFlowTransmissibility(model, 'Permeability'));
                end

                if model.dynamicFlowPv
                    % Dynamic porosity / pore volume
                    props = props.setStateFunction('Porosity', BactPorosity(model));
                    props = props.setStateFunction('PoreVolume', ...
                        DynamicFlowPoreVolume(model, 'Porosity'));
                end

                % Add microbial state functions
                props = props.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                props = props.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
                props = props.setStateFunction('BactConvRate', BactConvertionRate(model));
            end
        end

        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            % Call parent method for standard component conservation
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
        end

        %-----------------------------------------------------------------%
        function [acc, bflux, name, type] = bacteriaConservationEquation(fd, model, state, state0, dt)
            % Computes bacterial mass conservation equation
            bactmass  = model.getProps(state, 'BacterialMass');
            bactmass0 = model.getProps(state0, 'BacterialMass');

            % Accumulation term
            acc = (bactmass - bactmass0) ./ dt;

            % Convert to cell arrays for AD framework
            acc   = {acc};
            bflux = {[]};

            % Output variable names and types
            name = 'bacteria';
            type = 'cell';
        end

        %-----------------------------------------------------------------%
        function c = getPopulationMass(model, state, extra)
            % Compute microbial mass in each phase
            c = component.getPhaseComposition(model, state);

            if nargin < 4
                rho = model.getProps(state, 'Density');
            else
                rho = extra.rho;
            end

            for ph = 1:numel(c)
                if ~isempty(c{ph})
                    c{ph} = rho{ph} .* c{ph};
                end
            end
        end
    end
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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