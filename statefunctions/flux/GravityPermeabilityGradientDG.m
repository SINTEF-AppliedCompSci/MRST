classdef GravityPermeabilityGradientDG < StateFunction
   
    properties
        S
    end
    
    methods
        function gp = GravityPermeabilityGradientDG(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            model = varargin{1};
            dim = model.G.griddim;
            [ii, jj] = blockDiagIndex(dim*ones(dim,1), ones(dim,1));
            gp.S = sparse(ii, jj, 1);
        end
        
        function gRhoKdz = evaluateOnDomain(gp, model, state)
            
            rho = model.getProp(state, 'Density');
            
            g = model.getGravityVector;
            [K, ~, c] = permTensor(model.rock, model.G.griddim);
            gKdz = (K.*g(c))*gp.S;
            gKdz = mat2cell(gKdz, model.G.cells.num, ones(1,model.G.griddim));
            gKdz = SpatialVector(gKdz{:});
            gKdz = gKdz(state.cells,:);
            
            nph = numel(rho);
            gRhoKdz = cell(1,nph);
            for i = 1:nph
                gRhoKdz{i} = rho{i}.*gKdz;
            end
        end
        
    end
        
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
