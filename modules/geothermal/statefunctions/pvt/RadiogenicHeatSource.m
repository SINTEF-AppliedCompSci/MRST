classdef RadiogenicHeatSource < StateFunction
%State function for heat influx due to radiogenic heat production

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RadiogenicHeatSource(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PoreVolume'});
            rhoQh = model.radiogenicHeatFluxDensity;
            assert(numel(rhoQh) == 1 || numel(rhoQh) == model.G.cells.num, ...
                ['Radiogenic heat flux denisty must be specified ', ...
                 'either as scalar, or one value per grid cell']);
            assert(all(rhoQh > 0), ...
                'Radiogenic heat flux density must be positive');
            gp.label = 'Q_{h,rad}';
        end
        
        %-----------------------------------------------------------------%
        function qh = evaluateOnDomain(prop,model, state)
            
            pv    = prop.getEvaluatedDependencies(state, 'PoreVolume');
            vol   = model.operators.vol - pv;
            rhoQh = model.radiogenicHeatFluxDensity;
            qh    = vol.*rhoQh;
            
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