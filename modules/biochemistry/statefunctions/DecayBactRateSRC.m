classdef DecayBactRateSRC <  StateFunction
    % The bacterial decay rate, given per cell
    properties           
    end
    
    methods
        function gp = DecayBactRateSRC(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'s'}, 'state');
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            gp.label = 'Psi_decay';
        end

        function Psidecay = evaluateOnDomain(prop, model, state)
            Psidecay = 0;
            bbact = model.ReservoirModel.b_bact;
            nMax = model.ReservoirModel.nbactMax;
            namecp = model.ReservoirModel.getComponentNames();
            pv = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'PoreVolume');
            s = model.ReservoirModel.getProps(state, 's');
            nbact = model.ReservoirModel.getProps(state, 'nbact');
            rho = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'Density');
            L_ix = model.ReservoirModel.getLiquidIndex();
            if ~iscell(rho)
                rho = {rho};
            end
            idx_H2 = find(strcmp(namecp, 'H2'), 1);     % Locate 'H2'
            if iscell(s)
                sL = s{L_ix};
                rhoL = rho{L_ix};
            else
                sL = s(:,L_ix);
                rhoL = rho(:,L_ix);
            end
            if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase && (~isempty(idx_H2))
                if iscell(rhoL)
                    Voln = sL{1};
                else
                    Voln = sL;
                end
                Voln = max(Voln, 1.0e-8);
                Psidecay = pv.*bbact.*nbact.*(nbact.*Voln);
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
