classdef GrowthBactRateSRC <  StateFunction
    % The bacterial growth rate, given per cell
    properties
    end
    
    methods
        function gp = GrowthBactRateSRC(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'x'}, 'state');
            gp = gp.dependsOn({'s'}, 'state');
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            gp.label = 'Psi_growth';
        end

        function Psigrowth = evaluateOnDomain(prop, model, state)
            Psigrowth = 0;
            namecp = model.ReservoirModel.getComponentNames();
            % Find indices for 'H2' and 'CO2'
            idx_H2 = find(strcmp(namecp, 'H2'));     % Locate 'H2'
            idx_CO2 = find(strcmp(namecp, 'CO2'));   % Locate 'CO2'

            if model.ReservoirModel.bacteriamodel&& model.ReservoirModel.liquidPhase && (~isempty(idx_H2)) && (~isempty(idx_CO2))

                pv = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'PoreVolume');
                rho = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'Density');
                s = model.ReservoirModel.getProps(state, 's');
                nbact = model.ReservoirModel.getProps(state, 'nbact');


                L_ix = model.ReservoirModel.getLiquidIndex();

                if ~iscell(rho)
                    rho = {rho};
                end
                x = model.ReservoirModel.getProps(state, 'x');
                if iscell(x)
                    xH2 = x{idx_H2};     % Mole fraction of H2
                    xCO2 = x{idx_CO2};
                    sL = s{L_ix};
                    rhoL = rho{L_ix};
                else
                    xH2 = x(:, idx_H2);  % Column corresponding to H2 in matrix
                    xCO2 = x(:, idx_CO2);  % Column corresponding to CO2 in matrix
                    sL = s(:,L_ix);
                    rhoL = rho(:,L_ix);
                end
                if iscell(rhoL)
                    Voln = sL{1};
                else
                    Voln = sL;
                end
                
                Voln = max(Voln, 1.0e-8);
                alphaH2 = model.ReservoirModel.alphaH2;
                alphaCO2 = model.ReservoirModel.alphaCO2;
                Psigrowthmax = model.ReservoirModel.Psigrowthmax;
                nbactMax = model.ReservoirModel.nbactMax;
                bact_limit = 1 - (nbact./nbactMax).^0.5;
                % Calculate Psigrowth using H2 and CO2 mole fractions
                Psigrowth = pv.*Psigrowthmax.* (xH2 ./ (alphaH2 + xCO2)) ...
                    .* (xCO2 ./ (alphaCO2 + xH2)).*nbact.*Voln;
            end
        end
    end
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
