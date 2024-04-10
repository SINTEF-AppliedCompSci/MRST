classdef ConnectionPressureDrop < StateFunction
%State function for perforation pressure gradient using implicit connection
%pressure drop

    properties
        
    end
    
    methods
        
        function gp = ConnectionPressureDrop(varargin)
            
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            gp = gp.dependsOn('bhp', 'state');
            gp.label = 'p_c-p_{bh}-g \Delta z \rho_{w}';
            
        end
        
        function cdp = evaluateOnDomain(prop, model, state)
            
            
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
           
            [rho] = model.ReservoirModel.getProps(state, 'Density');
            rho = cellfun(@(rho) rho(map.cells), rho, 'UniformOutput', false);
            rho = rho{1};
            g = norm(model.ReservoirModel.gravity);
            
            nw = numel(map.W);
            cdp = cell(nw,1);
            
            for i = 1:nw
                dz = diff([0; map.W(i).dZ]);
                rhoW = rho(map.perf2well == i);
                np = numel(dz);
                S = ones(np);
                S = sparse(tril(S));
                cdpW = S*(g.*rhoW.*dz);
                cdp{i} = cdpW;
            end
            cdp = vertcat(cdp{:});
            
        end
        
    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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